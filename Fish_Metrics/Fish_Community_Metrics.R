# This script builds one row of community metrics per OpCode (survey replicate), using every habitat. 
# Elasmobranchs (sharks and rays) are excluded so the metrics reflect teleosts only. 
# For each OpCode it calculates abundance, biomass, species count, and Shannon diversity for all teleosts, targeted species only, and non-targeted species only, using the fisheries status lookup from SP_TRAITS.csv. 
# Site, Replicate, Habitat, and Year are carried through from the biomass data so the output can be filtered by habitat afterwards.

# Inputs : biomass.csv, SP_TRAITS.csv
# Output : community_metrics_fish.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr); library(purrr)
})

# Paths
biomass_path <- "input" # User's biomass csv
traits_path  <- "SP_TRAITS.csv" # User's species traits (i.e. Family, Genus, Species, Sci)


# 1. Load data ------------------------------------------------------------

fish_all <- read_csv(biomass_path, show_col_types = FALSE)
traits_raw <- read_csv(traits_path, show_col_types = FALSE)

fish <- fish_all  # no habitat filter - keep reef, Halls, HMAS Ledge together

message("Total fish observations: ", nrow(fish_all))

# 2. Exclude elasmobranchs (sharks and rays) ------------------------------

elasmobranch_families <- c(
  "Carcharhinidae",    # requiem sharks
  "Hemiscylliidae",    # bamboo sharks
  "Orectolobidae",     # wobbegongs
  "Aetobatidae",        # eagle rays
  "Dasyatidae"          # stingrays
)

elasmo_before <- fish %>%
  filter(Family %in% elasmobranch_families)

fish <- fish %>%
  filter(!Family %in% elasmobranch_families)

message("Removed ", nrow(elasmo_before), " elasmobranch observations - ",
        nrow(fish), " teleost observations remain across ",
        length(unique(fish$OpCode)), " OpCodes")

# 3. Helper functions -----------------------------------------------------

# Normalize a species name to "genus species", lowercase and squished, dropping
# any trailing qualifiers beyond the first two whitespace-separated tokens.
norm_sci <- function(x) {
  x %>%
    str_squish() %>%
    str_to_lower() %>%
    sub("^(\\S+)\\s+(\\S+).*$", "\\1 \\2", .)
}

# Shannon diversity index (natural log) for a vector of counts.
shannon_H <- function(counts) {
  counts <- counts[counts > 0]
  if (!length(counts)) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p))
}

# Build abundance, biomass, species count, and Shannon diversity per OpCode for whatever subset of fish rows is passed in.
build_metrics <- function(dfin) {
  sp_tab <- dfin %>%
    group_by(OpCode, Sci) %>%
    summarise(
      abundance = n(),
      biomass   = sum(biomass_g, na.rm = TRUE),
      .groups = "drop"
    )

  sp_tab %>%
    group_by(OpCode) %>%
    summarise(
      total_abundance   = sum(abundance),
      total_biomass     = sum(biomass, na.rm = TRUE),
      species_count     = n(),
      shannon_diversity = shannon_H(abundance),
      .groups = "drop"
    ) %>%
    mutate(
      total_biomass     = round(total_biomass, 3),
      shannon_diversity = round(shannon_diversity, 6)
    )
}

# Rename a metrics block's columns with a prefix (e.g. "all", "targeted").
rename_block <- function(tbl, tag) {
  tbl %>%
    dplyr::select(OpCode, total_abundance, total_biomass, shannon_diversity, species_count) %>%
    rename(
      !!paste0(tag, "_abundance")     := total_abundance,
      !!paste0(tag, "_biomass")       := total_biomass,
      !!paste0(tag, "_shannon")       := shannon_diversity,
      !!paste0(tag, "_species_count") := species_count
    )
}


# 4. Prepare fisheries status (targeted / non-targeted) from SP_TR --------

status_col <- c("Fisheries status", "Fisheries_status", "fisheries_status", "FisheriesStatus")
status_col <- status_col[status_col %in% names(traits_raw)]
if (!length(status_col)) stop("Could not find a fisheries status column in SP_TRAITS.csv")
status_col <- status_col[1]

traits_clean <- traits_raw %>%
  mutate(
    Sci_norm = norm_sci(sci),
    FisheriesStatus = toupper(str_squish(.data[[status_col]]))
  ) %>%
  mutate(
    FisheriesStatus = case_when(
      FisheriesStatus %in% c("T", "TARGETED") ~ "T",
      FisheriesStatus %in% c("NT", "NON-TARGETED", "NON TARGETED", "NONTARGETED", "NON TARGET") ~ "NT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Sci_norm)) %>%
  group_by(Sci_norm) %>%
  summarise(
    FisheriesStatus = {
      s <- FisheriesStatus[!is.na(FisheriesStatus)]
      if (!length(s)) NA_character_ else names(sort(table(s), decreasing = TRUE))[1]
    },
    .groups = "drop"
  )

fish2 <- fish %>%
  mutate(Sci_norm = norm_sci(Sci)) %>%
  left_join(traits_clean, by = "Sci_norm") %>%
  dplyr::select(-Sci_norm)

miss_stat <- fish2 %>% filter(is.na(FisheriesStatus)) %>% distinct(Sci) %>% arrange(Sci)
if (nrow(miss_stat)) {
  readr::write_csv(miss_stat, "species_missing_fisheries_status_all_habitats.csv")
  message(nrow(miss_stat), " species missing fisheries status - see species_missing_fisheries_status_all_habitats.csv")
}


# 5. Build metrics blocks (all, targeted, non-targeted) -------------------

m_all   <- build_metrics(fish2)
opc_all <- rename_block(m_all, "all")

df_T  <- fish2 %>% filter(FisheriesStatus == "T")
df_NT <- fish2 %>% filter(FisheriesStatus == "NT")

message("Fish counts - all: ", nrow(fish2), " | targeted: ", nrow(df_T), " | non-targeted: ", nrow(df_NT))

opc_T  <- if (nrow(df_T))  rename_block(build_metrics(df_T),  "targeted")    else NULL
opc_NT <- if (nrow(df_NT)) rename_block(build_metrics(df_NT), "nontargeted") else NULL


# 6. Consolidate to a single OpCode table ---------------------------------

blocks <- list(opc_all, opc_T, opc_NT)
blocks <- blocks[!vapply(blocks, is.null, logical(1))]
opcode_metrics <- reduce(blocks, ~ full_join(.x, .y, by = "OpCode")) %>%
  arrange(OpCode)


# 7. Add site metadata (Site, Replicate, Habitat, Year) from the f --------

site_info <- fish %>%
  group_by(OpCode) %>%
  slice(1) %>%
  dplyr::select(OpCode, Site, Replicate, Habitat, Year) %>%
  ungroup()

opcode_metrics <- opcode_metrics %>%
  left_join(site_info, by = "OpCode")


# 8. Write output (one row per OpCode) ------------------------------------

final_order <- c(
  "OpCode", "Site", "Replicate", "Habitat", "Year",
  "all_biomass", "all_abundance", "all_shannon", "all_species_count",
  "targeted_biomass", "targeted_abundance", "targeted_shannon", "targeted_species_count",
  "nontargeted_biomass", "nontargeted_abundance", "nontargeted_shannon", "nontargeted_species_count"
)

for (nm in final_order) {
  if (!nm %in% names(opcode_metrics)) opcode_metrics[[nm]] <- NA_real_
}

opcode_metrics <- opcode_metrics %>% dplyr::select(all_of(final_order))

write_csv(opcode_metrics, "community_metrics_fish.csv", na = "")

message("Saved community_metrics_fish.csv - ", nrow(opcode_metrics), " OpCodes")

# 9. Family-level abundance and biomass for targeted families (carangidae,
#    sparidae, lutjanidae, serranidae, lethrinidae), restricted to the OpCodes
#    validated in all_final.csv. Reads the biomass file fresh rather than
#    reusing opcode_metrics above, since it needs individual-fish rows.
#    Output: family_abundances.csv, family_biomass.csv

all_fish <- read.csv("All_Fish_Biomass.csv")
final_dat <- read.csv("all_final.csv")

valid_opcodes <- unique(final_dat$opcode)

all_fish <- all_fish %>%
  mutate(opcode = tolower(OpCode)) %>%
  filter(opcode %in% valid_opcodes)

if (nrow(all_fish) == 0) {
  stop("No rows remain after opcode filtering - check that opcodes match between files.")
}

target_families <- c(
  "carangidae",
  "sparidae",
  "lutjanidae",
  "serranidae",
  "lethrinidae"
)

all_fish <- all_fish %>%
  mutate(family_lower = tolower(Family)) %>%
  filter(family_lower %in% target_families)

if (nrow(all_fish) == 0) {
  stop("No rows remain after family filtering - check target_families against the Family column.")
}

message("Family-filtered rows: ", nrow(all_fish), " across families: ",
        paste(unique(all_fish$family_lower), collapse = ", "))

# 9a. Abundance per opcode x family, wide format
family_abundance_wide <- all_fish %>%
  group_by(opcode, family_lower) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = family_lower,
    values_from = abundance,
    values_fill = 0
  )

for (fam in target_families) {
  if (!fam %in% names(family_abundance_wide)) {
    family_abundance_wide[[fam]] <- 0
  }
}

family_abundance_wide <- family_abundance_wide %>%
  dplyr::select(opcode, dplyr::all_of(target_families))

write.csv(family_abundance_wide, "family_abundances.csv", row.names = FALSE)
message("Saved family_abundances.csv - ", nrow(family_abundance_wide), " opcodes")

# 9b. Biomass per opcode x family, wide format
family_biomass_wide <- all_fish %>%
  group_by(opcode, family_lower) %>%
  summarise(biomass = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = family_lower,
    values_from = biomass,
    values_fill = 0
  )

for (fam in target_families) {
  if (!fam %in% names(family_biomass_wide)) {
    family_biomass_wide[[fam]] <- 0
  }
}

family_biomass_wide <- family_biomass_wide %>%
  dplyr::select(opcode, dplyr::all_of(target_families))

write.csv(family_biomass_wide, "family_biomass.csv", row.names = FALSE)
message("Saved family_biomass.csv - ", nrow(family_biomass_wide), " opcodes")
