
# Fish Biomass Pipeline ---------------------------------------------------


# Purpose:
#   Converts underwater visual survey data (per-individual lengths + MaxN
#   counts) into a per-fish biomass estimate, then rolls that up into
#   summary tables by species, OpCode (survey replicate), Site, and Habitat.
#
# Required input files (in the working directory):
#   - ALL_Lengths.csv 
#       One row per *measured* fish. Required columns:
#       OpCode, Site, Replicate, Habitat, Year, Sci, "Length (mm)"
#       Extracted from EventMeasure
#   - ALL_MaxN.csv
#       One row per (OpCode x species) with the MaxN abundance count.
#       Required columns:
#       OpCode, Site, Replicate, Habitat, Year, Family, Genus, Species, Sci, MaxN
#       Extracted from EventMeasure
#   - manual_Sci_to_speccode.csv (created by hand, see Stage 3)
#       Columns: Sci, SpecCode
#   - manual_ab_overrides.csv (created by hand, see Stage 5)
#       Columns (any subset): Sci, a, b, lw_type, alpha, beta
#
# Pipeline stages:
#   1. Read and validate raw inputs
#   2. Expand each (OpCode x species) MaxN count into individual fish rows, using measured lengths first and a median-based imputation ladder (replicate -> site -> habitat -> study) to fill any shortfall
#   3. Map each species (Sci) to a FishBase SpecCode
#   4. Merge in hand-curated SpecCodes for species FishBase couldn't match
#   5. Pull length-weight (a, b) and length-length (alpha, beta) regression parameters from FishBase for each species
#   6. Apply manual overrides for species with missing/unsuitable FishBase parameters
#   7. Compute biomass per fish, with strict validation that every fish has a value (or the script stops and reports what's missing)
#   8. Aggregate biomass/abundance into summary tables
#
# Output files:
#   fish_rows.csv, fish_rows_with_speccode.csv, Sci_to_speccode.csv, All_Fish_Biomass.csv, Manual_Lengths_or_Params_Needed.csv (only if validation fails), biomass_by_species.csv, biomass_by_OpCode.csv, biomass_by_Site.csv, biomass_by_Habitat.csv, species_list.csv, species_by_OpCode_LONG.csv, OpCode_by_species_abundance_biomass_WIDE.csv


# ---- Packages --------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(rfishbase)
library(stringr)
library(tidyr)

# ---- Configuration -----------------------------------------------------
# Length type of your field measurements: "FL" (fork length), "TL" (total length), or "SL" (standard length). Used to (a) prefer FishBase length-weight parameters measured in this type, and (b) decide when a length-length conversion is needed before applying (a, b).
LENGTH_TYPE_INPUT <- "FL"

MM_TO_CM <- 0.1          # field lengths are recorded in mm; L-W equations expect cm
FETCH_PAUSE_SEC <- 0.8   # delay between FishBase API calls, to stay under rate limits
FETCH_RETRIES   <- 2     # retry attempts per species before giving up

`%||%` <- function(x, y) if (!is.null(x)) x else y  # null-coalescing helper



# Stage 1: Read inputs ----------------------------------------------------



lengths <- read_csv("ALL_Lengths.csv", show_col_types = FALSE)
maxn    <- read_csv("ALL_MaxN.csv",    show_col_types = FALSE)

# Use Sci as the species key; lengths column is "Length (mm)"
stopifnot(all(c("OpCode","Site","Replicate","Habitat","Year","Sci","Length (mm)") %in% names(lengths)))
stopifnot(all(c("OpCode","Site","Replicate","Habitat","Year","Family","Genus","Species","Sci","MaxN") %in% names(maxn)))

# Coerce types
maxn <- maxn %>% mutate(MaxN = as.integer(MaxN))
lengths <- lengths %>% mutate(`Length (mm)` = as.numeric(`Length (mm)`))


# Stage 2: Expand MaxN counts into individual fish --------


# ---- 2a) Index the measured lengths by (OpCode, Sci) --------------------
# Keep only measured lengths (non-NA), preserving file order, since Stage 2c takes the first N lengths per group in that order.
lengths_clean <- lengths %>%
  filter(!is.na(`Length (mm)`))

# Named list: key = "OpCode|Sci" -> vector of that group's lengths, in file order
group_lengths_df <- lengths_clean %>%
  group_by(OpCode, Sci) %>%
  summarise(vals = list(`Length (mm)`), .groups = "drop") %>%
  mutate(key = paste(OpCode, Sci, sep = "|"))
length_map <- setNames(group_lengths_df$vals, group_lengths_df$key)

# ---- 2b) Median imputation ladder ----------------------------------------
# When an OpCode has more MaxN individuals than measured lengths, the shortfall is filled with a median length. The median is drawn from the most specific group available, falling back to broader groups: replicate (OpCode) -> site -> habitat -> whole study
# This keeps imputed lengths as locally representative as possible while still producing a value for species measured elsewhere in the study.

replicate_med_df <- lengths_clean %>%
  group_by(OpCode, Sci) %>%
  summarise(med = median(`Length (mm)`), .groups = "drop") %>%
  mutate(key = paste(OpCode, Sci, sep = "|"))
rep_med_vec <- setNames(replicate_med_df$med, replicate_med_df$key)

site_med_df <- lengths_clean %>%
  group_by(Site, Sci) %>%
  summarise(med = median(`Length (mm)`), .groups = "drop") %>%
  mutate(key = paste(Site, Sci, sep = "|"))
site_med_vec <- setNames(site_med_df$med, site_med_df$key)

habitat_med_df <- lengths_clean %>%
  group_by(Habitat, Sci) %>%
  summarise(med = median(`Length (mm)`), .groups = "drop") %>%
  mutate(key = paste(Habitat, Sci, sep = "|"))
hab_med_vec <- setNames(habitat_med_df$med, habitat_med_df$key)

study_med_df <- lengths_clean %>%
  group_by(Sci) %>%
  summarise(med = median(`Length (mm)`), .groups = "drop")
study_med_vec <- setNames(study_med_df$med, study_med_df$Sci)


# A single numeric length in mm, or NA_real_ if the species has no measured lengths anywhere in the study.
get_impute_value <- function(op, site, hab, Sci) {
  k_rep  <- paste(op,   Sci, sep = "|")
  k_site <- paste(site, Sci, sep = "|")
  k_hab  <- paste(hab,  Sci, sep = "|")
  if (!is.na(rep_med_vec[k_rep]))   return(as.numeric(rep_med_vec[k_rep]))
  if (!is.na(site_med_vec[k_site])) return(as.numeric(site_med_vec[k_site]))
  if (!is.na(hab_med_vec[k_hab]))   return(as.numeric(hab_med_vec[k_hab]))
  if (!is.na(study_med_vec[Sci]))   return(as.numeric(study_med_vec[Sci]))
  return(NA_real_)  # species has no lengths anywhere -> leave missing
}

# ---- 2c) Expand each (OpCode x Sci) to MaxN individual-fish rows --------
#    take the first N measured lengths (file order)
#    if there's a shortfall, fill it using the median ladder (2b)
#    if the species has no lengths anywhere, leave those rows NA

maxn_filtered <- maxn %>% filter(!is.na(MaxN), MaxN > 0)

# Expand one MaxN record into MaxN individual-fish rows.
# return A tibble with `MaxN` rows, one per individual, carrying the original grouping columns plus an assigned `Length (mm)` (measured where available, imputed otherwise, NA if never measured anywhere).
expand_one <- function(OpCode, Site, Replicate, Habitat, Year, Family, Genus, Species, Sci, MaxN) {
  key <- paste(OpCode, Sci, sep = "|")
  avail <- length_map[[key]] %||% numeric(0)

  taken <- head(avail, MaxN)
  remaining <- MaxN - length(taken)

  if (remaining > 0) {
    imp <- get_impute_value(OpCode, Site, Habitat, Sci)
    imps <- rep(imp, remaining)
    lens <- c(taken, imps)
  } else {
    lens <- taken
  }

  tibble(
    OpCode    = OpCode,
    Site      = Site,
    Replicate = Replicate,
    Habitat   = Habitat,
    Year      = Year,
    Family    = Family,
    Genus     = Genus,
    Species   = Species,
    Sci       = Sci,
    `Length (mm)` = round(lens, 2)
  )
}

fish_rows <- pmap(
  maxn_filtered[c("OpCode","Site","Replicate","Habitat","Year","Family","Genus","Species","Sci","MaxN")],
  expand_one
) %>% list_rbind()

write_csv(fish_rows, "fish_rows.csv")

# Diagnostic: species where imputation failed for every one of their rows (i.e. that species has no measured length anywhere in the study)
species_with_nas <- fish_rows %>%
  group_by(Sci) %>%
  summarise(
    total = n(),
    n_missing = sum(is.na(`Length (mm)`)),
    .groups = "drop"
  ) %>%
  filter(total == n_missing)  # keep only species where *all* rows are NA

print(species_with_nas)


# Stage 3: Map species (Sci) to FishBase SpecCode -------------------------


# Normalize a species name to "genus species", lowercase and squished.
# Drops any trailing qualifiers (e.g. "cf.", subspecies) beyond the first two whitespace-separated tokens.

norm_Sci <- function(x) x %>% str_squish() %>% str_to_lower() %>% sub("^(\\S+)\\s+(\\S+).*$", "\\1 \\2", .)

# Convert a normalized lowercase binomial ("lutjanus sebae") to the Title Case form ("Lutjanus sebae") that rfishbase expects for name queries.
to_title_binom <- function(x) {
  x <- norm_Sci(x)
  parts <- str_split_fixed(x, "\\s+", 2)
  paste(stringr::str_to_title(parts[, 1]), parts[, 2], sep = " ")
}

# Build the unique list of species to query (binomials only -- excludes unresolved "sp."/"spp." entries, which FishBase can't look up by name)
unique_spp <- fish_rows %>% distinct(Sci)
name_vec <- unique_spp$Sci %>% norm_Sci() %>% unique()
name_vec <- name_vec[grepl("^[a-z]+\\s+[a-z]+$", name_vec) & !grepl("\\bsp\\.?\\b|\\bspp\\.?\\b", name_vec)]
name_vec_title <- to_title_binom(name_vec)

cat("Will query", length(name_vec_title), "species by name, one-by-one.\n")

# Look up one species' SpecCode on FishBase, falling back to SeaLifeBase.


# Fetching species codes from FishBase takes the longest to run -----------


fetch_one <- function(s, pause_sec = FETCH_PAUSE_SEC, retries = FETCH_RETRIES) {
  for (i in 0:retries) {
    Sys.sleep(pause_sec)
    out <- tryCatch(
      rfishbase::species(s, fields = c("SpecCode","Genus","Species"), server = "fishbase"),
      error = function(e) NULL
    )
    if (!is.null(out) && is.data.frame(out) && nrow(out) > 0) return(out)
    # try sealifebase as fallback (invertebrates / species FishBase doesn't cover)
    out <- tryCatch(
      rfishbase::species(s, fields = c("SpecCode","Genus","Species"), server = "sealifebase"),
      error = function(e) NULL
    )
    if (!is.null(out) && is.data.frame(out) && nrow(out) > 0) return(out)
  }
  NULL
}

res_list <- lapply(name_vec_title, fetch_one)  # slow step: one API round-trip per species
ok_idx <- vapply(res_list, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))

cat("Succeeded for", sum(ok_idx), "of", length(ok_idx), "species.\n")

sp <- if (any(ok_idx)) dplyr::bind_rows(res_list[ok_idx]) else tibble()

if (nrow(sp)) {
  # FishBase's `Species` field is inconsistent across builds: sometimes just
  # "sebae", sometimes the full "Lutjanus sebae". Normalize both cases.
  to_valid_norm <- function(genus, species_col) {
    sc <- ifelse(grepl("\\s", species_col), species_col, paste(genus, species_col))
    norm_Sci(sc)
  }

  sp_clean <- sp %>%
    dplyr::mutate(
      Valid_norm = to_valid_norm(Genus, Species),
      SpecCode   = suppressWarnings(as.integer(SpecCode))
    ) %>%
    dplyr::filter(!is.na(Valid_norm), !is.na(SpecCode)) %>%
    dplyr::distinct(Valid_norm, .keep_all = TRUE) %>%
    dplyr::select(Valid_norm, SpecCode)

  # Join: our Sci -> normalized name -> SpecCode
  Sci_map <- unique_spp %>%
    dplyr::mutate(Sci_norm = norm_Sci(Sci)) %>%
    dplyr::left_join(sp_clean, by = c("Sci_norm" = "Valid_norm")) %>%
    dplyr::distinct(Sci, Sci_norm, SpecCode)

  readr::write_csv(Sci_map, "Sci_to_speccode.csv")

  fish_rows_with_speccode <- fish_rows %>%
    dplyr::left_join(dplyr::select(Sci_map, Sci, SpecCode), by = "Sci")

  readr::write_csv(fish_rows_with_speccode, "fish_rows_with_speccode.csv")

  # Report a few misses so they can be inspected / added to the manual file
  misses <- Sci_map %>% dplyr::filter(is.na(SpecCode)) %>% head(20)
  if (nrow(misses)) {
    cat("Examples of unmapped species (first 20):\n")
    print(misses)
  } else {
    cat("All queried species mapped to SpecCode.\n")
  }
} else {
  
}

# Species FishBase couldn't match automatically -- these are the candidates for the manual_Sci_to_speccode.csv file used in Stage 4.
Sci_map <- read_csv("Sci_to_speccode.csv", show_col_types = FALSE)
unmatched <- Sci_map %>%
  filter(is.na(SpecCode)) %>%
  select(Sci, SpecCode)
# To generate a starting template for manual lookups, uncomment:
write_csv(unmatched, "manual_Sci_to_speccode_template.csv")


# Merge hand-curated SpecCodes for unmatched species ----------------------


# manual_Sci_to_speccode.csv is a hand-edited file (created from `unmatched` above) with two columns: Sci, SpecCode. Fill in SpecCode by looking the species up manually on fishbase, then re-run this section.

# 1) Load the FishBase-derived mapping built in Stage 3
Sci_map_fb <- read_csv("fish_rows_with_speccode.csv", show_col_types = FALSE)

# 2) Load the manual file and normalize names
manual <- readr::read_csv("manual_Sci_to_speccode.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    Sci      = stringr::str_squish(Sci),
    SpecCode = as.integer(SpecCode),
    Sci_norm = tolower(Sci)
  ) %>%
  dplyr::filter(!is.na(Sci), !is.na(SpecCode)) %>%
  dplyr::distinct(Sci_norm, .keep_all = TRUE) %>%       # de-dup by normalized key
  dplyr::select(Sci_norm, SpecCode_man = SpecCode)

# 3) Normalize the FishBase map the same way, for joining
Sci_map_fb2 <- Sci_map_fb %>%
  mutate(Sci_norm = str_to_lower(str_squish(Sci))) %>%
  distinct(Sci, Sci_norm, SpecCode)

# 4) Merge without losing existing codes: manual takes precedence where present, otherwise keep whatever FishBase returned
Sci_map_final <- Sci_map_fb2 %>%
  dplyr::left_join(manual, by = "Sci_norm") %>%
  dplyr::mutate(SpecCode = dplyr::coalesce(SpecCode_man, SpecCode)) %>%
  dplyr::select(Sci, SpecCode) %>%
  dplyr::distinct(Sci, .keep_all = TRUE)

write_csv(Sci_map_final, "Sci_to_speccode.csv")

# Re-attach the merged SpecCodes to fish_rows
fish_rows <- read_csv("fish_rows.csv", show_col_types = FALSE)
fish_rows_with_speccode <- fish_rows %>%
  left_join(Sci_map_final, by = "Sci")
write_csv(fish_rows_with_speccode, "fish_rows_with_speccode.csv")

# Coverage check: what fraction of species now have a SpecCode?
Sci_map_final %>% summarise(
  n_species_total = n_distinct(Sci),
  n_with_code     = sum(!is.na(SpecCode)),
  pct_mapped      = round(100 * n_with_code / n_species_total, 1)
) %>% print()

message("Manual SpecCodes merged. Updated Sci_to_speccode.csv and fish_rows_with_speccode.csv.")


# Stage 5: Length-weight (a, b) and length-length (alpha, beta) pa --------


# For each species, pull FishBase's length-weight equation W = a * L^b, and (if the species' equation isn't in our preferred length type) a length-length regression to convert our input length type into that type first.
# Within-type reduction when multiple FishBase records exist for the same (SpecCode, length type): 1 record -> use it; 2 records -> mean; 3+ -> median. Across length types, the species' own LENGTH_TYPE_INPUT is preferred, falling back to the other two types if unavailable.

# Return the first name in `candidates` found in `nm`, or NA.
# Used to match differently-cased/named columns across rfishbase versions.
pick_available <- function(nm, candidates) {
  for (c in candidates) if (c %in% nm) return(c)
  NA_character_
}

# ---- 5a) Pull and normalize length-weight records ------------------------
lw_raw <- tryCatch(rfishbase::length_weight(fields = c("SpecCode","a","b","LType")),
                    error = function(e) rfishbase::length_weight())
nm_lw <- tolower(names(lw_raw)); names(lw_raw) <- nm_lw

sc_col <- pick_available(nm_lw, c("speccode","spec_code","spec.code"))
a_col  <- pick_available(nm_lw, c("a","coefa","intercept"))
b_col  <- pick_available(nm_lw, c("b","coefb","slope"))
lt_col <- pick_available(nm_lw, c("ltype","type","lengthtype","length_type"))

stopifnot(!is.na(sc_col), !is.na(a_col), !is.na(b_col))  # must have SpecCode, a, b

lw <- lw_raw %>%
  dplyr::transmute(
    SpecCode = suppressWarnings(as.integer(.data[[sc_col]])),
    a        = suppressWarnings(as.numeric(.data[[a_col]])),
    b        = suppressWarnings(as.numeric(.data[[b_col]])),
    lw_type  = if (!is.na(lt_col))
      toupper(trimws(as.character(.data[[lt_col]])))
    else NA_character_
  ) %>%
  dplyr::filter(!is.na(SpecCode), !is.na(a), !is.na(b)) %>%
  dplyr::mutate(lw_type = dplyr::if_else(is.na(lw_type), "TL", lw_type)) %>%  # default to TL if missing
  dplyr::filter(lw_type %in% c("FL","TL","SL"))

# ---- 5b) Reduce multiple records per (SpecCode, lw_type) to one ----------
lw_reduced_by_type <- lw %>%
  dplyr::group_by(SpecCode, lw_type) %>%
  dplyr::summarise(
    n_rec   = dplyr::n(),
    a_mean  = mean(a, na.rm = TRUE),
    a_med   = stats::median(a, na.rm = TRUE),
    a_first = dplyr::first(a),
    b_mean  = mean(b, na.rm = TRUE),
    b_med   = stats::median(b, na.rm = TRUE),
    b_first = dplyr::first(b),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    a = dplyr::case_when(n_rec == 1 ~ a_first, n_rec == 2 ~ a_mean, TRUE ~ a_med),
    b = dplyr::case_when(n_rec == 1 ~ b_first, n_rec == 2 ~ b_mean, TRUE ~ b_med)
  ) %>%
  dplyr::select(SpecCode, lw_type, a, b)

# ---- 5c) Choose one length type per SpecCode, preferring our input type --
type_pref <- c(LENGTH_TYPE_INPUT, setdiff(c("FL","TL","SL"), LENGTH_TYPE_INPUT))
lw_by_speccode <- lw_reduced_by_type %>%
  dplyr::mutate(pref_order = match(lw_type, type_pref)) %>%
  dplyr::arrange(SpecCode, pref_order) %>%
  dplyr::group_by(SpecCode) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(SpecCode, a, b, lw_type)

# ---- 5d) Map SpecCode -> Sci, collapsing multiple SpecCodes per Sci ------
stopifnot(all(c("Sci","SpecCode") %in% names(Sci_map)))

lw_by_sci_all_types <- Sci_map %>%
  dplyr::filter(!is.na(SpecCode)) %>%
  dplyr::left_join(lw_by_speccode, by = "SpecCode") %>%
  dplyr::filter(!is.na(a), !is.na(b), lw_type %in% c("FL","TL","SL")) %>%
  dplyr::group_by(Sci, lw_type) %>%
  dplyr::summarise(
    a = stats::median(a, na.rm = TRUE),  # across SpecCodes mapped to the same Sci & type
    b = stats::median(b, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 5e) Pick the preferred length type per species -----------------------
lw_by_type <- lw_by_sci_all_types %>%
  dplyr::mutate(pref_order = match(lw_type, type_pref)) %>%
  dplyr::arrange(Sci, pref_order) %>%
  dplyr::group_by(Sci) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sci, a, b, lw_type)

# ---- 5f) Length-length regressions (for converting our input type into whatever type a species' length-weight equation uses)
# Pull and normalize FishBase's length-length regression table.
# Reduces multiple records per (SpecCode, from, to) pair to their median.
# A data frame of SpecCode/from/to/alpha/beta, or NULL if the table is empty or FishBase's response is missing an expected column.
get_llr <- function() {
  tbl <- tryCatch(rfishbase::length_length(), error = function(e) NULL)
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  nm <- tolower(names(tbl)); names(tbl) <- nm

  sc_col   <- pick_available(nm, c("speccode","spec_code","spec.code"))
  from_col <- pick_available(nm, c("type1","ltype1","length1type","length_type1"))
  to_col   <- pick_available(nm, c("type2","ltype2","length2type","length_type2"))
  a_col    <- pick_available(nm, c("a","intercept"))
  b_col    <- pick_available(nm, c("b","slope"))

  if (any(is.na(c(sc_col, from_col, to_col, a_col, b_col)))) return(NULL)

  tbl %>%
    dplyr::transmute(
      SpecCode = suppressWarnings(as.integer(.data[[sc_col]])),
      from     = toupper(trimws(as.character(.data[[from_col]]))),
      to       = toupper(trimws(as.character(.data[[to_col]]))),
      alpha    = suppressWarnings(as.numeric(.data[[a_col]])),
      beta     = suppressWarnings(as.numeric(.data[[b_col]]))
    ) %>%
    dplyr::filter(!is.na(SpecCode), !is.na(alpha), !is.na(beta),
                  from %in% c("FL","TL","SL"), to %in% c("FL","TL","SL")) %>%
    dplyr::group_by(SpecCode, from, to) %>%
    dplyr::summarise(alpha = stats::median(alpha, na.rm = TRUE),
                      beta  = stats::median(beta,  na.rm = TRUE), .groups = "drop")
}

llr_map <- NULL
llr_raw <- get_llr()
if (!is.null(llr_raw)) {
  llr_map <- Sci_map %>%
    dplyr::filter(!is.na(SpecCode)) %>%
    dplyr::left_join(llr_raw, by = "SpecCode") %>%
    dplyr::group_by(Sci, from, to) %>%
    dplyr::summarise(alpha = stats::median(alpha, na.rm = TRUE),
                      beta  = stats::median(beta,  na.rm = TRUE), .groups = "drop")
}

# ---- 5g) Join (a, b) and LLR onto fish_rows -------------------------------
fish_rows <- fish_rows %>%
  dplyr::select(-dplyr::any_of(c("a","b","lw_type","alpha","beta","length_cm_for_lw"))) %>%
  dplyr::left_join(lw_by_type, by = "Sci") %>%
  dplyr::mutate(input_type = LENGTH_TYPE_INPUT)

if (!is.null(llr_map)) {
  fish_rows <- fish_rows %>%
    dplyr::left_join(llr_map, by = c("Sci" = "Sci", "input_type" = "from", "lw_type" = "to"))
} else {
  fish_rows <- fish_rows %>% dplyr::mutate(alpha = NA_real_, beta = NA_real_)
}

# Species without a usable LLR fall through to Stage 6/7's manual overrides or the missing-data report.


# Stage 6: Manual overrides for length-weight / length-length para --------


# manual_ab_overrides.csv is a hand-edited file for species where FishBase has no usable regression, or where its regression should be overridden.
# Expected columns (any subset may be present): Sci, a, b, lw_type, alpha, beta
#   - a, b:         length-weight coefficients (W = a * L^b) in `lw_type` units
#   - alpha, beta:  length-length regression converting LENGTH_TYPE_INPUT -> lw_type


has_col <- function(df, nm) nm %in% names(df)  # does column `nm` exist in `df`?

manual_raw <- read_csv("manual_ab_overrides.csv")
manual_ab <- manual_raw %>%
  dplyr::mutate(
    Sci = str_squish(Sci),

    a_man     = if (has_col(., "a"))     readr::parse_number(str_remove(as.character(a),     "^'\\s*")) else NA_real_,
    b_man     = if (has_col(., "b"))     readr::parse_number(str_remove(as.character(b),     "^'\\s*")) else NA_real_,
    alpha_man = if (has_col(., "alpha")) readr::parse_number(str_remove(as.character(alpha), "^'\\s*")) else NA_real_,
    beta_man  = if (has_col(., "beta"))  readr::parse_number(str_remove(as.character(beta),  "^'\\s*")) else NA_real_,

    lw_type_man = if (has_col(., "lw_type"))
      toupper(str_squish(as.character(lw_type)))
    else NA_character_
  ) %>%
  dplyr::filter(!is.na(Sci) & Sci != "") %>%
  dplyr::mutate(Sci_norm = tolower(Sci)) %>%    # normalize key so the join below is one-to-one
  dplyr::group_by(Sci_norm) %>%
  dplyr::slice_tail(n = 1) %>%                  # if duplicates exist, keep the last entry
  dplyr::ungroup() %>%
  dplyr::select(Sci, a_man, b_man, lw_type_man, alpha_man, beta_man)

# ---- 6a) Collapse any remaining duplicate Sci entries ---------------------
dupes <- manual_ab %>% count(Sci, name = "n") %>% filter(n > 1)
if (nrow(dupes) > 0) {
  message("Found ", nrow(dupes), " species with duplicate overrides; collapsing to one per Sci.")

  first_non_na <- function(x) { i <- which(!is.na(x))[1]; if (is.na(i)) NA else x[i] }
  mode_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]  # most common non-NA value
  }

  manual_ab <- manual_ab %>%
    group_by(Sci) %>%
    summarise(
      a_man       = first_non_na(a_man),
      b_man       = first_non_na(b_man),
      lw_type_man = mode_non_na(lw_type_man),
      alpha_man   = first_non_na(alpha_man),
      beta_man    = first_non_na(beta_man),
      .groups = "drop"
    )
}

# ---- 6b) Apply overrides ---------------------------------------------------
prefer_manual <- TRUE  # TRUE = manual values always win; FALSE = only fill gaps FishBase left

fish_rows <- fish_rows %>%
  dplyr::left_join(manual_ab, by = "Sci") %>%
  dplyr::mutate(
    a = case_when(
      prefer_manual & !is.na(a_man) ~ a_man,
      !prefer_manual & is.na(a) & !is.na(a_man) ~ a_man,
      TRUE ~ a
    ),
    b = case_when(
      prefer_manual & !is.na(b_man) ~ b_man,
      !prefer_manual & is.na(b) & !is.na(b_man) ~ b_man,
      TRUE ~ b
    ),
    lw_type = case_when(
      prefer_manual & !is.na(lw_type_man) ~ lw_type_man,
      !prefer_manual & is.na(lw_type) & !is.na(lw_type_man) ~ lw_type_man,
      TRUE ~ lw_type
    ),
    alpha = case_when(
      prefer_manual & !is.na(alpha_man) ~ alpha_man,
      !prefer_manual & is.na(alpha) & !is.na(alpha_man) ~ alpha_man,
      TRUE ~ alpha
    ),
    beta = case_when(
      prefer_manual & !is.na(beta_man) ~ beta_man,
      !prefer_manual & is.na(beta) & !is.na(beta_man) ~ beta_man,
      TRUE ~ beta
    )
  ) %>%
  dplyr::select(-dplyr::ends_with("_man"))


# Stage 7: Biomass per fish, with strict validation -----------------------


# 0) length_mm must exist
if (!"length_mm" %in% names(fish_rows)) {
  if ("Length (mm)" %in% names(fish_rows)) {
    fish_rows <- fish_rows %>% mutate(length_mm = round(as.numeric(`Length (mm)`), 2))
  } else stop("Need 'length_mm' or 'Length (mm)' in fish_rows for biomass.")
}

# 1) Compute biomass: W = a * L^b, where L is in the length type the species' (a, b) equation expects (converting via alpha/beta if needed)
fish_rows <- fish_rows %>%
  mutate(
    input_type       = LENGTH_TYPE_INPUT,
    length_cm_input  = length_mm * MM_TO_CM,
    length_cm_for_lw = case_when(
      is.na(lw_type) ~ NA_real_,                                       # no target type known
      lw_type == input_type ~ length_cm_input,                          # already in the right type
      !is.na(alpha) & !is.na(beta) ~ alpha + beta * length_cm_input,    # convert input_type -> lw_type
      TRUE ~ NA_real_                                                   # type mismatch, no LLR available
    ),
    biomass_g = if_else(!is.na(length_cm_for_lw) & !is.na(a) & !is.na(b),
                         a * (length_cm_for_lw ^ b), NA_real_)
  )

# 2) Strict validation: every fish must end up with a biomass value. If any don't, write out *why* for each one and stop 
na_cases <- fish_rows %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      is.na(length_mm) ~ "No length_mm (imputation ladder failed)",
      is.na(a) | is.na(b) ~ "Missing a/b (FishBase + manual)",
      lw_type != input_type & (is.na(alpha) | is.na(beta)) ~ "Type mismatch: FL->(TL/SL) LLR missing",
      is.na(lw_type) ~ "Missing lw_type",
      is.na(length_cm_for_lw) ~ "Converted length missing (LLR or type)",
      TRUE ~ "Other/compound"
    )
  ) %>%
  dplyr::filter(is.na(biomass_g)) %>%
  dplyr::select(
    OpCode, Site, Replicate, Habitat, Year,
    Family, Genus, Species, Sci,
    length_mm, input_type, lw_type, a, b, alpha, beta,
    length_cm_for_lw, biomass_g, reason
  ) %>%
  dplyr::arrange(OpCode, Site, Replicate, Sci)

if (nrow(na_cases) > 0) {
  out_file <- "Manual_Lengths_or_Params_Needed.csv"
  readr::write_csv(na_cases, out_file, na = "NA")
  message("Biomass missing for ", nrow(na_cases), " rows. See: ", normalizePath(out_file))
  stop("Biomass must exist for every fish row. Patch missing lengths / a,b / LLR in manual_ab_overrides.csv and re-run.")
}

readr::write_csv(fish_rows, "All_Fish_Biomass.csv")

# Quick sanity check: top 10 individual fish by biomass
top_individuals <- fish_rows %>%
  arrange(desc(biomass_g)) %>%
  slice_head(n = 10) %>%
  dplyr::select(biomass_g)
print(top_individuals)


# Stage 8: Summary tables -------------------------------------------------

# ---- 8a) Biomass by species -----------------------------------------------
biomass_by_species <- fish_rows %>%
  filter(!is.na(biomass_g)) %>%
  group_by(Sci) %>%
  summarise(
    n_individuals    = n(),
    total_biomass_g  = sum(biomass_g, na.rm = TRUE),
    mean_biomass_g   = mean(biomass_g, na.rm = TRUE),
    median_biomass_g = median(biomass_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_biomass_g))
readr::write_csv(biomass_by_species, "biomass_by_species.csv")

# ---- 8b) Biomass by OpCode (survey replicate) -----------------------------
stopifnot("biomass_g" %in% names(fish_rows))

biomass_by_OpCode <- fish_rows %>%
  group_by(OpCode) %>%
  summarise(
    total_biomass_g = sum(biomass_g, na.rm = TRUE),
    n_individuals   = sum(!is.na(biomass_g)),
    n_rows          = n(),
    n_missing       = sum(is.na(biomass_g)),
    .groups = "drop"
  ) %>%
  arrange(OpCode)

biomass_by_OpCode_export <- biomass_by_OpCode %>%
  mutate(
    total_biomass_g  = round(total_biomass_g, 3),
    total_biomass_kg = round(total_biomass_g / 1000, 3)
  )
write_csv(biomass_by_OpCode_export, "biomass_by_OpCode.csv", na = "")

# ---- 8c) Biomass by Site ---------------------------------------------------
stopifnot("biomass_g" %in% names(fish_rows))

biomass_by_Site <- fish_rows %>%
  group_by(Site) %>%
  summarise(
    total_biomass_g = sum(biomass_g, na.rm = TRUE),
    n_individuals   = sum(!is.na(biomass_g)),
    n_rows          = n(),
    n_missing       = sum(is.na(biomass_g)),
    .groups = "drop"
  ) %>%
  arrange(Site)

biomass_by_Site_export <- biomass_by_Site %>%
  mutate(
    total_biomass_g  = round(total_biomass_g, 3),
    total_biomass_kg = round(total_biomass_g / 1000, 3)
  )
write_csv(biomass_by_Site_export, "biomass_by_Site.csv", na = "")

# ---- 8d) Biomass by Habitat -------------------------------------------------
stopifnot("biomass_g" %in% names(fish_rows))

biomass_by_Habitat <- fish_rows %>%
  group_by(Habitat) %>%
  summarise(
    total_biomass_g = sum(biomass_g, na.rm = TRUE),
    n_individuals   = sum(!is.na(biomass_g)),
    n_rows          = n(),
    n_missing       = sum(is.na(biomass_g)),
    .groups = "drop"
  ) %>%
  arrange(Habitat)

biomass_by_Habitat_export <- biomass_by_Habitat %>%
  mutate(
    total_biomass_g  = round(total_biomass_g, 3),
    total_biomass_kg = round(total_biomass_g / 1000, 3)
  )
write_csv(biomass_by_Habitat_export, "biomass_by_Habitat.csv", na = "")

# ---- 8e) Species list (Family, Genus, Species, Sci) ------------------------
species_list <- fish_rows %>%
  distinct(Family, Genus, Species, Sci) %>%
  arrange(Family, Genus, Species, Sci)
write_csv(species_list, "species_list.csv")

# ---- 8f) Biomass and abundance per species x OpCode -----------------------
fish <- readr::read_csv("All_Fish_Biomass.csv", show_col_types = FALSE)

req <- c("OpCode","Site","Replicate","Habitat","Sci","biomass_g")
stopifnot(all(req %in% names(fish)))

# Environmental lookup: one row per OpCode
env_lookup <- fish %>%
  dplyr::distinct(OpCode, Site, Replicate, Habitat)

# Aggregate to species x OpCode
#   abundance = number of individual rows (expanded from MaxN in Stage 2)
#   biomass_g = summed biomass over those individuals
species_opc_long <- fish %>%
  dplyr::group_by(OpCode, Sci) %>%
  dplyr::summarise(
    abundance = dplyr::n(),
    biomass_g = sum(biomass_g, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(species_opc_long, "species_by_OpCode_LONG.csv")

# Wide format: two columns per species, <Species>_abundance and <Species>_biomass_g. Column names will contain spaces where Sci does; that's fine for CSV/Excel.
species_opc_wide_only <- species_opc_long %>%
  tidyr::pivot_wider(
    names_from  = Sci,
    values_from = c(abundance, biomass_g),
    names_glue  = "{Sci}_{.value}",
    values_fill = list(abundance = 0L, biomass_g = 0)
  )

species_opc_wide <- env_lookup %>%
  dplyr::left_join(species_opc_wide_only, by = "OpCode") %>%
  dplyr::arrange(OpCode)

readr::write_csv(species_opc_wide, "OpCode_by_species_abundance_biomass_WIDE.csv")

