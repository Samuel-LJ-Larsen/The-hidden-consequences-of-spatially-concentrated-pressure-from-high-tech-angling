

# Fish Community dbRDA — Replicate Level

# Tests whether fish assemblage composition differs between habitat types at the replicate level, and identifies which environmental predictors and species drive that structure. Uses distance-based redundancy analysis (dbRDA/capscale) constrained by habitat + environmental predictors, cross-checked with PERMANOVA, and visualised as a dbRDA plot.

# Steps:
  #   1. Load master csv, filter to habitats of interest, check required columns
  #   2. Identify environmental predictor columns present in the data
  #   3. Identify species abundance columns (genus_species_abundance)
  #   4. Build replicate x species matrix, zero-fill NAs, drop all-zero rows
  #   5. Transform abundances (sqrt) and filter out rare taxa
  #   6. Join environmental predictors, drop rows with missing values
  #   7. Check dispersion homogeneity (betadisper), fit dbRDA (capscale), test significance, check VIFs, cross-check with PERMANOVA
  #   8. Test each predictor individually for its own R^2/significance
  #   9. Extract site/species ordination scores; select species to plot via envfit
  #  10. Select significant environmental predictor arrows via envfit
  #  11. Plot dbRDA biplot with habitat ellipses, species and predictor arrows

# Outputs: replicate_predictor_stats.csv, species_envfit_results_replicate.csv, fish_dbrda_replicate_sqrt.pdf/.svg


# Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(dplyr)
})


# 0. File and options -----------------------------------------------------


csv_path <- "input" # Input master csv file

habitats_keep     <- c("indurated_reef", "rocky_reef") # Input habitat's relevant to user
min_reps_present  <- 3        # rarity filter: taxa present in >= this many replicates
transform_method  <- "sqrt"   

# envfit-driven species selection (by r^2 only; for replicate plot)
r2_cutoff          <- 0.29
max_species_labels <- 80
plot_predictors    <- TRUE
set.seed(1)



# Helpers for taxon labels ------------------------------------------------


pretty_taxon <- function(x) {
  y <- sub("_abundance$", "", x, ignore.case = TRUE)
  y <- gsub("[_.]+", " ", y)
  y <- trimws(gsub("\\s+", " ", y))
  vapply(strsplit(y, " "), function(tokens) {
    tokens <- tolower(tokens)
    if (length(tokens) > 0) {
      tokens[1] <- paste0(toupper(substring(tokens[1], 1, 1)), substring(tokens[1], 2))
    }
    paste(tokens, collapse = " ")
  }, character(1))
}

.cap1 <- function(s) paste0(toupper(substr(s,1,1)), tolower(substring(s,2)))

pretty_scientific <- function(x) {
  y <- sub("_abundance$", "", x, ignore.case = TRUE)
  y <- gsub("[_.]+", " ", y)
  y <- trimws(gsub("\\s+", " ", y))
  vapply(strsplit(y, " "), function(tok) {
    if (length(tok) == 0) return("")
    tok <- tolower(tok)
    tok[1] <- .cap1(tok[1])
    paste(tok, collapse = " ")
  }, character(1))
}


# 1. Load and basic checks ------------------------------------------------

dat_raw <- read.csv(csv_path, check.names = FALSE)

# ID columns
id_cols <- c("opcode","site","replicate","habitat")
stopifnot(all(id_cols %in% names(dat_raw)))

# Keep only habitats of interest
dat_raw <- dat_raw %>% dplyr::filter(habitat %in% habitats_keep)
if (nrow(dat_raw) == 0) stop("No rows left after filtering to habitats: ",
                             paste(habitats_keep, collapse=", "))

# Habitat factor + colours
dat_raw$habitat <- factor(dat_raw$habitat,
                          levels = c("indurated_reef","rocky_reef"))
cols <- c("indurated_reef" = "#d55e00", "rocky_reef" = "#0069a2")


# 2. Environmental predictors (replicate level) ---------------------------


want_preds <- c(
  "coral_cover",
  "coral_shannon",
  "distance_to_estuary",
  "depth",
  "slope_mean_500m",
  "New_Log_RS",
  "prox_index_rocky_reef_500m_log",
  "relief_max_500m"
)

env_pred_cols <- intersect(want_preds, names(dat_raw))  # only those present
if (!length(env_pred_cols)) {
  stop("None of the requested predictors found in the data:\n",
       paste(want_preds, collapse = ", "))
}
cat("Using predictors (replicate-level):\n")
print(env_pred_cols)


# 3. Species columns (genus_species_abundance)  ---------------------------


species_regex <- "^[A-Za-z][A-Za-z\\.]*_[A-Za-z\\.]+_abundance$"
non_taxa_prefixes <- c("^all_", "^targeted_", "^nontargeted_", "^gt\\d+_",
                       "^lt\\d+_", "^site_", "^opc_", "^total_", "^family_",
                       "^richness_", "^diversity_")

cand_species <- names(dat_raw)[grepl(species_regex, names(dat_raw))]
cand_species <- setdiff(
  cand_species,
  grep(paste(non_taxa_prefixes, collapse="|"),
       cand_species, value = TRUE, ignore.case = TRUE)
)

num_cols     <- names(dat_raw)[sapply(dat_raw, is.numeric)]
species_cols <- intersect(cand_species, num_cols)
species_cols <- setdiff(species_cols, c(id_cols, env_pred_cols))

if (!length(species_cols)) {
  message("No species columns matched. Diagnostics:")
  message(" - Regex candidates: ", length(cand_species))
  if (length(cand_species))
    message("   e.g.: ", paste(utils::head(cand_species, 10), collapse=", "))
  not_numeric <- setdiff(cand_species, num_cols)
  if (length(not_numeric)) {
    message(" - Matched but not numeric: ",
            paste(utils::head(not_numeric, 10), collapse=", "))
  }
  stop("No usable species columns found. Check names and types.")
}

cat("Detected species columns (n=", length(species_cols), ")\n", sep = "")
print(utils::head(species_cols, 10))


# 4. Build replicate tables --------------------------------------------------


rep_species <- dat_raw[, c(id_cols, species_cols)]
rep_env     <- dat_raw[, c(id_cols, env_pred_cols), drop = FALSE]

# Species matrix (replicate-level)
sp_rep <- rep_species[, species_cols, drop = FALSE]

# Replace NA with 0 (counts/abundance)
sp_rep[is.na(sp_rep)] <- 0

# Assert non-negative
neg_cols <- names(which(colSums(sp_rep < 0, na.rm = TRUE) > 0))
if (length(neg_cols)) {
  stop("Species matrix must be non-negative. Negatives in: ",
       paste(neg_cols, collapse=", "))
}

# Drop all-zero rows BEFORE transform and keep metadata aligned
keep_rep_rows <- rowSums(sp_rep) > 0
if (!all(keep_rep_rows))
  message("Dropping ", sum(!keep_rep_rows),
          " replicate rows with all-zero species.")

sp_rep  <- sp_rep [keep_rep_rows, , drop = FALSE]
dat_rep <- rep_species[keep_rep_rows, c("opcode","site","replicate","habitat")]


# 5. Transform species (sqrt) ---------------------------------------------


if (identical(transform_method, "hellinger")) {
  sp_rep_t <- vegan::decostand(sp_rep, method = "hellinger")
} else if (identical(transform_method, "sqrt")) {
  # Simple sqrt of abundance
  sp_rep_t <- sqrt(sp_rep)
} else if (is.null(transform_method)) {
  sp_rep_t <- as.matrix(sp_rep)
} else {
  stop("Unknown transform_method: ", transform_method,
       " (use 'hellinger', 'sqrt', or NULL)")
}

# Rarity filter at replicate level
keep_taxa <- colSums(sp_rep_t > 0, na.rm = TRUE) >= min_reps_present
if (!any(keep_taxa))
  stop("All taxa filtered by 'min_reps_present' at replicate level.")
sp_rep_f <- sp_rep_t[, keep_taxa, drop = FALSE]

cat("Replicates kept:", nrow(dat_rep),
    "  Taxa kept:", ncol(sp_rep_f), "\n")


# 6. Replicate-level site/env table (predictors) --------------------------

# Bring predictors into the replicate data 
rep_env_only <- dat_raw %>%
  dplyr::select(dplyr::any_of(c(id_cols, env_pred_cols))) %>%
  dplyr::distinct()

# Align rows with species matrix after zero-row filter
dat_rep_env <- dplyr::left_join(dat_rep, rep_env_only,
                                by = c("opcode","site","replicate","habitat"))

# Drop replicates with missing predictors
cc <- stats::complete.cases(dat_rep_env[, env_pred_cols, drop = FALSE])
if (!all(cc))
  message("Dropping ", sum(!cc),
          " replicate rows with NA in predictors for replicate-level dbRDA.")

dat_rep_env <- dat_rep_env[cc, , drop = FALSE]
sp_rep_f_cc <- sp_rep_f[cc, , drop = FALSE]

# 7. Replicate-level dispersion + dbRDA + PERMANOVA -----------------------

# (A) Dispersion (homogeneity) by habitat
D_rep  <- vegdist(sp_rep_f_cc, method = "bray")
bd_rep <- betadisper(D_rep, dat_rep_env$habitat)
cat("betadisper (replicate level) ANOVA:\n")
print(anova(bd_rep))
cat("\nPermutation test for dispersion (replicate level):\n")
print(permutest(bd_rep, permutations = 999))

# (B) dbRDA formula (replicate level)
form_rep <- as.formula(
  paste("sp_rep_f_cc ~ habitat +", paste(env_pred_cols, collapse = " + "))
)

dbr_rep <- capscale(form_rep, data = dat_rep_env, distance = "bray")

cat("\nReplicate-level dbRDA tests:\n")
print(anova(dbr_rep, permutations = 999))                 # overall
print(anova(dbr_rep, by = "terms", permutations = 999))   # by term
print(anova(dbr_rep, by = "axis", permutations = 999))    # constrained axes
cat("Adj. R^2 (replicate-level):",
    RsquareAdj(dbr_rep)$adj.r.squared, "\n\n")

# Collinearity check among constraints
suppressWarnings(print(vegan::vif.cca(dbr_rep)))

# PERMANOVA (same fixed set, replicate level)
cat("REPLICATE-LEVEL PERMANOVA (adonis2):\n")
print(adonis2(form_rep, data = dat_rep_env,
              permutations = 999, method = "bray"))

# 8. Predictor statistics table -------------------------------------------

terms_anova <- anova(dbr_rep, by = "terms", permutations = 999)
terms_df    <- as.data.frame(terms_anova)
adj_r2_full <- RsquareAdj(dbr_rep)$adj.r.squared

predictor_stats <- data.frame()

for (pred in c("habitat", env_pred_cols)) {
  form_single <- as.formula(paste("sp_rep_f_cc ~", pred))
  mod_single <- tryCatch(
    capscale(form_single, data = dat_rep_env, distance = "bray"),
    error = function(e) NULL
  )
  if (!is.null(mod_single)) {
    adj_r2_single <- RsquareAdj(mod_single)$adj.r.squared
    anova_single  <- anova(mod_single, permutations = 999)
    predictor_stats <- rbind(predictor_stats, data.frame(
      Predictor   = pred,
      Adj_R2      = adj_r2_single,
      df          = anova_single$Df[1],
      F           = anova_single$F[1],
      p           = anova_single$`Pr(>F)`[1],
      Significant = ifelse(anova_single$`Pr(>F)`[1] < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
}

predictor_stats <- predictor_stats[order(predictor_stats$p), ]


cat("Replicate level predictor statistics (sorted by significance)\n")
print(predictor_stats, row.names = FALSE)

write.csv(predictor_stats,
          "replicate_predictor_stats.csv",
          row.names = FALSE)
cat("\nSaved: replicate_predictor_stats.csv\n")
cat("Full model Adjusted R²:", round(adj_r2_full, 4), "\n")


# 9. Scores + species ENFIT -----------------------------------------------

eig_con_rep <- vegan::eigenvals(dbr_rep, model = "constrained")
if (length(eig_con_rep) && sum(eig_con_rep) > 0) {
  x_lab_rep <- sprintf("dbRDA1 (%.1f%% of constrained)",
                       100 * eig_con_rep[1] / sum(eig_con_rep))
  y_lab_rep <- if (length(eig_con_rep) >= 2)
    sprintf("dbRDA2 (%.1f%% of constrained)",
            100 * eig_con_rep[2] / sum(eig_con_rep)) else "dbRDA2"
} else {
  x_lab_rep <- "dbRDA1"; y_lab_rep <- "dbRDA2"
}

site_sc <- as.data.frame(scores(dbr_rep, display = "sites", choices = 1:2))
site_sc$habitat <- dat_rep_env$habitat
colnames(site_sc)[1:2] <- c("dbRDA1","dbRDA2")

sp_sc <- as.data.frame(scores(dbr_rep, display = "species", choices = 1:2))
if (nrow(sp_sc) > 0 && ncol(sp_sc) >= 2) {
  colnames(sp_sc)[1:2] <- c("dbRDA1","dbRDA2")
  sp_sc$Species_raw <- rownames(sp_sc)
  sp_sc$Species     <- pretty_taxon(sp_sc$Species_raw)
  sp_sc$arrow_len   <- sqrt(sp_sc$dbRDA1^2 + sp_sc$dbRDA2^2)
  have_sp <- TRUE
} else {
  sp_sc <- data.frame(
    dbRDA1 = numeric(0),
    dbRDA2 = numeric(0),
    Species_raw = character(0),
    Species = character(0),
    arrow_len = numeric(0),
    stringsAsFactors = FALSE
  )
  have_sp <- FALSE
}

# ----- envfit species selection -----
if (have_sp) {
  fit <- vegan::envfit(dbr_rep, sp_rep_f_cc, permutations = 9999)
  
  vec_raw <- tryCatch(rownames(fit$vectors$arrows), error = function(e) NULL)
  r_vec   <- tryCatch(as.numeric(fit$vectors$r),     error = function(e) numeric(0))
  p_vec   <- tryCatch(as.numeric(fit$vectors$pvals), error = function(e) numeric(0))
  
  if (!is.null(vec_raw) && length(r_vec) == length(vec_raw) && length(r_vec) > 0) {
    r2 <- r_vec
    stats_tbl <- data.frame(
      Species_raw = vec_raw,
      Species     = pretty_taxon(vec_raw),
      r2 = r2,
      p_value = if (length(p_vec) == length(vec_raw)) p_vec else rep(NA, length(vec_raw)),
      keep_r2 = r2 >= r2_cutoff,
      row.names = NULL
    )
    cat("\nENVFIT species r^2 and p-values (replicate ordination):\n")
    print(stats_tbl[order(-stats_tbl$r2),
                    c("Species","r2","p_value","keep_r2")])
    
    sig_raw <- stats_tbl$Species_raw[stats_tbl$keep_r2]
    sp_lab <- dplyr::semi_join(sp_sc,
                               data.frame(Species_raw = sig_raw),
                               by = "Species_raw")
    
    if (nrow(sp_lab) == 0 && max_species_labels > 0) {
      sp_lab <- sp_sc %>%
        dplyr::arrange(dplyr::desc(arrow_len)) %>%
        dplyr::slice_head(n = max_species_labels)
    }
    
    write.csv(stats_tbl[order(-stats_tbl$r2), ],
              "species_envfit_results_replicate.csv",
              row.names = FALSE)
    cat("\n✓ Saved: species_envfit_results_replicate.csv\n")
  } else {
    sp_lab <- sp_sc[0, ]
  }
} else {
  sp_lab <- sp_sc[0, ]
}


# X/Y limits including ellipses + species ---------------------------------

# Build a temporary plot to obtain full ellipse coordinates
p_tmp <- ggplot() +
  stat_ellipse(data = site_sc,
               aes(dbRDA1, dbRDA2, fill = habitat),
               type = "norm", geom = "polygon",
               alpha = 0.18, color = NA)

gb   <- ggplot_build(p_tmp)

# gather ellipse coordinates from all groups
ellx <- unlist(lapply(gb$data, function(d) d$x))
elly <- unlist(lapply(gb$data, function(d) d$y))

# combine with sites + species arrows + origin
x_all <- c(site_sc$dbRDA1,
           if (nrow(sp_lab)) sp_lab$dbRDA1 * 0.8 else numeric(0),
           ellx,
           0)
y_all <- c(site_sc$dbRDA2,
           if (nrow(sp_lab)) sp_lab$dbRDA2 * 0.8 else numeric(0),
           elly,
           0)

x_all <- x_all[is.finite(x_all)]
y_all <- y_all[is.finite(y_all)]
rng_x <- range(x_all); rng_y <- range(y_all)
padx  <- diff(rng_x) * 0.08
pady  <- diff(rng_y) * 0.08
xlim  <- c(rng_x[1] - padx, rng_x[2] + padx)
ylim  <- c(rng_y[1] - pady, rng_y[2] + pady)



# 10. ENVFIT for predictors -----------------------------------------------

pred_sig_alpha <- 0.05
env_sig <- data.frame()

env_mat <- as.data.frame(dat_rep_env[, env_pred_cols, drop = FALSE])
env_mat[] <- lapply(env_mat, function(col_) suppressWarnings(as.numeric(col_)))

fit_env <- vegan::envfit(dbr_rep, env_mat, permutations = 9999)

env_ar_raw <- tryCatch(fit_env$vectors$arrows, error = function(e) NULL)
env_p      <- tryCatch(as.numeric(fit_env$vectors$pvals), error = function(e) NULL)

env_ar_df <- if (!is.null(env_ar_raw))
  as.data.frame(env_ar_raw, stringsAsFactors = FALSE, check.names = FALSE) else NULL

if (!is.null(env_ar_df)) {
  num_cols <- which(vapply(env_ar_df, is.numeric, logical(1)))
  if (length(num_cols) >= 2) {
    env_ar_df$ax1 <- as.numeric(env_ar_df[[num_cols[1]]])
    env_ar_df$ax2 <- as.numeric(env_ar_df[[num_cols[2]]])
  } else env_ar_df <- NULL
}

if (!is.null(env_ar_df) && !is.null(env_p) && nrow(env_ar_df) == length(env_p)) {
  env_ar_df$var <- rownames(env_ar_df)
  env_ar_df$p   <- env_p
  env_ar_df$sig <- env_ar_df$p <= pred_sig_alpha
  
  env_sig <- env_ar_df[env_ar_df$sig, c("ax1","ax2","var","p","sig"),
                       drop = FALSE]
  
  if (nrow(env_sig) > 0) {
    # Scale arrows to current xlim/ylim (which already include ellipses)
    max_ax1 <- max(abs(env_sig$ax1), na.rm = TRUE); if (!is.finite(max_ax1) || max_ax1 == 0) max_ax1 <- 1
    max_ax2 <- max(abs(env_sig$ax2), na.rm = TRUE); if (!is.finite(max_ax2) || max_ax2 == 0) max_ax2 <- 1
    span_x <- diff(xlim) / 2
    span_y <- diff(ylim) / 2
    scaleF <- 0.9 * min(span_x / max_ax1, span_y / max_ax2)
    
    env_sig$x0 <- 0;  env_sig$y0 <- 0
    env_sig$x1 <- env_sig$ax1 * scaleF
    env_sig$y1 <- env_sig$ax2 * scaleF
    env_sig$label <- gsub("_", " ", env_sig$var)
  } else {
    env_sig <- data.frame()
  }
} else {
  env_sig <- data.frame()
}

# Optional label overrides
norm_key <- function(s) { s <- tolower(gsub("[. ]", "_", s)); trimws(s) }
predictor_label_map <- c(
  "coral_cover"                  = "Coral cover",
  "coral_shannon"                = "Coral evenness",
  "distance_to_estuary"          = "Distance to estuary",
  "depth"                        = "Depth",
  "prox_index_rocky_reef_500m_log" = "Rocky reef proximity index",
  "relief_max_500m"             = "Relief",
  "slope_mean_500m"             = "Slope",
  "New_Log_RS" = "Fishing pressure index"
)
if (nrow(env_sig) > 0) {
  keys <- norm_key(env_sig$var)
  hit  <- keys %in% names(predictor_label_map)
  env_sig$label[hit] <- unname(predictor_label_map[keys[hit]])
}


# 11. Plot ----------------------------------------------------------------


p <- ggplot() +
  stat_ellipse(data = site_sc,
               aes(dbRDA1, dbRDA2, fill = habitat),
               type = "norm", geom = "polygon",
               alpha = 0.18, color = NA) +
  geom_point(data = site_sc,
             aes(dbRDA1, dbRDA2, color = habitat),
             size = 2, alpha = 0.9) +
  scale_fill_manual(values = cols, name = "Habitat",
                    labels = c("Indurated reef","Rocky reef")) +
  scale_color_manual(values = cols, name = "Habitat",
                     labels = c("Indurated reef","Rocky reef"))

if (nrow(sp_lab) > 0) {
  p <- p +
    geom_segment(data = sp_lab,
                 aes(x = 0, y = 0,
                     xend = dbRDA1 * 0.8,
                     yend = dbRDA2 * 0.8),
                 arrow = arrow(length = unit(0.20, "cm")),
                 linewidth = 0.25, alpha = 0.6) +
    geom_text_repel(data = sp_lab,
                    aes(dbRDA1, dbRDA2, label = Species),
                    size = 2.6,
                    max.overlaps = Inf,
                    force = 3,
                    force_pull = 0.5,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0,
                    segment.size = 0.2,
                    max.iter = 10000)
}

if (plot_predictors && nrow(env_sig)) {
  p <- p +
    geom_segment(data = env_sig,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 arrow = arrow(length = unit(0.20, "cm")),
                 linewidth = 0.4, alpha = 0.8, colour = "grey25") +
    geom_text_repel(data = env_sig,
                    aes(x = x1, y = y1, label = label),
                    size = 2.8, min.segment.length = 0,
                    box.padding = 0.2, point.padding = 0.1)
}

p <- p +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = xlim, expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = ylim, expand = expansion(mult = 0.02)) +
  labs(title = "",
       x = x_lab_rep, y = y_lab_rep) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p)
ggsave("fish_dbrda_replicate_sqrt.pdf", plot = p,
       width = 10, height = 8, device = "pdf")

ggsave("fish_dbrda_replicate_sqrt.svg", plot = p,
       width = 10, height = 8, device = "svg")


cat("\nReplicate-level dbRDA — axis test:\n")
print(anova(dbr_rep, by = "axis", permutations = 999))
cat("\nReplicate-level VIFs:\n")
suppressWarnings(print(vegan::vif.cca(dbr_rep)))
