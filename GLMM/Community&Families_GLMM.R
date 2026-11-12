
# GLMMTMB Analysis: Targeted vs Non-targeted Fish Metrics by Habitat

# Purpose:
#   Fits a mixed-effects model (GLMM via glmmTMB) for each of 8 response variables (targeted/non-targeted x abundance/biomass/richness/Shannon), testing whether the response differs between two habitat types ("rocky_reef" vs "indurated_reef")
#   For each response it runs diagnostics (Anova, DHARMa residual checks), extracts pairwise habitat contrasts (emmeans), and produces a bar-plot panel with compact-letter-display (CLD) significance groupings.

# Required input file (in the working directory):
#   - Master csv
#       One row per replicate/site. Must contain:
#         - a habitat column with values including "rocky_reef" and "indurated_reef".
#         - a site column.
#         - the 8 response columns named in `family_config` below (targeted_abundance, nontargeted_abundance, targeted_biomass, nontargeted_biomass, targeted_richness, nontargeted_richness, targeted_shannon, nontargeted_shannon)

# Pipeline stages:
#   1. Load data, clean blank/unnamed columns, detect habitat/site column names, and filter to just the two habitats of interest.
#   2. Manually configured error-family choices per response variable (`family_config`) - these are assumed to already be decided from prior diagnostic exploration (not performed in this script).
#   3. Helper functions:
#        - standardise_contrast_cols(): normalizes emmeans contrast-table column names across the different possible model/link types.
#        - fit_with_family(): builds the response transform (if any) and fits the glmmTMB model for a given family choice.
#   4. fit_one(): for each response variable, subsets the data, fits the model, prints Anova + DHARMa diagnostics, extracts pairwise habitat contrast via emmeans, and returns a tidy one-row-per-contrast table.
#   5. Runs fit_one() over all 8 response variables and writes the combined pairwise-contrast results to CSV.
#   6. Plotting: create_formatted_plot() re-fits each response's model, computes estimated marginal means + compact letter display, and builds a styled bar chart with error bars and CLD significance letters. All 8 plots are assembled into a 4x2 panel and saved as SVG.

# This is all then repeated for the biomass of key fisheries targeted families. 

# Output files:
#   glmmtmb_pairwise_targeted_nontargeted.csv  (pairwise contrast results)
#   glmmtmb_targeted_nontarget_panel.png       (8-panel bar chart figure)


# ----- Load Required Packages -----
library(dplyr)          # Data manipulation
library(tidyr)          # Data tidying
library(stringr)        # String operations
install.packages("glmmTMB", type = "source")       # Generalized linear mixed models
library(DHARMa)         # Residual diagnostics for GLMMs
library(emmeans)        # Estimated marginal means
library(car)            # ANOVA for models
library(purrr)          # Functional programming tools
library(ggplot2)        # Plotting
library(multcomp)       # Multiple comparisons
library(multcompView)   # Compact letter display
library(patchwork)      # Combining plots


# ----- Set Analysis Configuration -----
DH_NSIM      <- 250     # Number of simulations for DHARMa diagnostics
TMB_CTRL     <- glmmTMBControl(optCtrl = list(iter.max = 2e4, eval.max = 2e4))
options(glmmTMB.verbose = FALSE)



# Stage 1: Load & prep data -----------------------------------------------

# ----- Load Data -----
dat <- read.csv("input") # Input user's master csv file


# ----- Clean Column Names -----
# Some exported CSVs (e.g. from Excel) carry blank/NA header names for
# spacer or leftover columns. This block finds truly-empty columns (all NA or all-whitespace) under a blank header and drops them, then gives any remaining blank-named columns a synthetic "unnamed_<index>" name so later
# code never has to deal with "" as a column name.
blank_idx <- which(names(dat) == "" | is.na(names(dat)))
if (length(blank_idx)) {
  empty_cols <- blank_idx[sapply(blank_idx, function(i) {
    v <- dat[[i]]
    all(is.na(v) | (is.character(v) && trimws(v) == ""))
  })]
  if (length(empty_cols)) dat <- dat[, -empty_cols, drop = FALSE]
  remaining_idx <- which(names(dat) == "" | is.na(names(dat)))
  if (length(remaining_idx)) names(dat)[remaining_idx] <- paste0("unnamed_", remaining_idx)
}

# ----- Identify Core Columns and Filter Data -----
# Habitat and Site column names vary between exports (e.g. "habitat" vs "habitat.new"), so this picks the first matching alias present in `dat`. stopifnot() guards against an ambiguous/missing column rather than silently proceeding with the wrong one.
hab_col  <- intersect(c("habitat","habitat.new","habitat_new"), names(dat))[1]
site_col <- intersect(c("site","site.name","site_name"),       names(dat))[1]
stopifnot(length(hab_col) == 1, length(site_col) == 1)

# Restrict to just the two habitats being compared, and set an explicit
# factor level order (rocky_reef first) so it becomes the reference level
# in model output and the baseline for emmeans contrasts below.
dat <- dat %>%
  filter(.data[[hab_col]] %in% c("rocky_reef","indurated_reef")) %>%
  mutate(
    !!hab_col  := factor(.data[[hab_col]], levels = c("rocky_reef","indurated_reef")),
    !!site_col := factor(.data[[site_col]])
  )

message("Data after habitat filter:")
print(table(dat[[hab_col]]))

# Stage 2: Manual family configuration (edit these to whatever you want)  --------

family_config <- list(
  # Abundance (counts)
  "targeted_abundance"   = "nb2",
  "nontargeted_abundance"  = "nb2",

  # Biomass (positive, skewed)
  "targeted_biomass"     = "gaussian_log1p",
  "nontargeted_biomass"    = "gaussian_log1p",

   # Richness (counts)
  "targeted_richness"    = "nb2",
  "nontargeted_richness"   = "nb2",
  
  #  Shannon (continuous, approx normal)
   "targeted_shannon"     = "gaussian",
   "nontargeted_shannon"    = "gaussian"
  
  
)

# ----- Get Response Variables -----
# The names of family_config double as the authoritative list of response
# columns to analyze (Stage 4) and plot (Stage 6).
analysis_cols <- names(family_config)
message("Analyzing ", length(analysis_cols), " response variables with manual family choices")



# Stage 3: Helper functions for model fitting and result tidying ----------


# ----- Helper: Standardize Contrast Column Names -----
# Normalize emmeans pairwise-contrast output to a common column schema.
# 'emmeans' contrasts output uses different column names depending on the model family/link (e.g. "ratio" for log-link count models, "odds.ratio" for logit-link models, "emmean" for identity-link Gaussian models, "z.ratio" vs "t.ratio" depending on whether df are estimated). This renames whichever variant is present into a single consistent set of columns (estimate, t.ratio, SE, df, contrast, p.value) so results from every response/family combination can be row-bound into one table.

standardise_contrast_cols <- function(tbl) {
  tbl <- as.data.frame(tbl)
  
  if (!("estimate" %in% names(tbl))) {
    if ("ratio" %in% names(tbl))           tbl <- dplyr::rename(tbl, estimate = ratio)
    else if ("odds.ratio" %in% names(tbl)) tbl <- dplyr::rename(tbl, estimate = odds.ratio)
    else if ("response" %in% names(tbl))   tbl <- dplyr::rename(tbl, estimate = response)
    else if ("emmean" %in% names(tbl))     tbl <- dplyr::rename(tbl, estimate = emmean)
    else tbl$estimate <- NA_real_
  }
  
  if (!("t.ratio" %in% names(tbl))) {
    if ("z.ratio" %in% names(tbl))   tbl <- dplyr::rename(tbl, t.ratio = z.ratio)
    else if ("stat" %in% names(tbl)) tbl <- dplyr::rename(tbl, t.ratio = stat)
    else tbl$t.ratio <- NA_real_
  }
  
  if (!("SE" %in% names(tbl)))       tbl$SE      <- NA_real_
  if (!("df" %in% names(tbl)))       tbl$df      <- NA_real_
  if (!("contrast" %in% names(tbl))) tbl$contrast <- "contrast"
  if (!("p.value" %in% names(tbl)))  tbl$p.value <- NA_real_
  
  tbl
}

# ----- Helper: Fit Model with Specified Family -----
# Build the response transform (if any) and fit a glmmTMB model.

# Family-choice reference:
#   - "beta"/"beta_disp": rescale Y into (0,1) (handles exact 0/1 via the standard n-observation correction), fit Beta regression (with or without a habitat-dependent dispersion formula).
#   - "lognormal": log(Y), Gaussian.
#   - "gaussian_log1p": log1p(Y), Gaussian (safe for Y containing zeros).
#   - "gaussian_logit": qlogis(Y), Gaussian.
#   - "poisson"/"nb2"/"zi_poisson"/"zi_nb2": count models (negative
#     binomial "nb2" variants allow overdispersion; "zi_" adds a zero-
#     inflation intercept term).
#   - "gamma"/"tweedie"/"tweedie_disp"/"zi_tweedie": continuous positive
#     (possibly zero-inflated / heteroscedastic) response models.
#   - "gaussian": untransformed Gaussian.

fit_with_family <- function(df, family_choice) {
  
  if (family_choice %in% c("beta", "beta_disp")) {
    n <- nrow(df)
    if (any(df$Y %in% c(0,1))) {
      df$Y_trans <- (df$Y * (n - 1) + 0.5) / n
    } else {
      df$Y_trans <- pmin(pmax(df$Y, .Machine$double.eps), 1 - .Machine$double.eps)
    }
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "lognormal") {
    df$Y_trans <- log(df$Y)
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "gaussian_log1p") {
    df$Y_trans <- log1p(df$Y)
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "gaussian_logit") {
    df$Y_trans <- qlogis(df$Y)
    formula_lhs <- "Y_trans"
    
  } else {
    formula_lhs <- "Y"
  }
  
  fml <- as.formula(paste(formula_lhs, "~ habitat + (1 | site)"))
  
  mod <- switch(
    family_choice,
    "poisson"        = glmmTMB(fml, data = df, family = poisson(link = "log"),  control = TMB_CTRL),
    "nb2"            = glmmTMB(fml, data = df, family = nbinom2(link = "log"),  control = TMB_CTRL),
    "zi_poisson"     = glmmTMB(fml, data = df, family = poisson(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "zi_nb2"         = glmmTMB(fml, data = df, family = nbinom2(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "gamma"          = glmmTMB(fml, data = df, family = Gamma(link = "log"),    control = TMB_CTRL),
    "tweedie"        = glmmTMB(fml, data = df, family = tweedie(link = "log"),  control = TMB_CTRL),
    "tweedie_disp"   = glmmTMB(fml, data = df, family = tweedie(link = "log"),  dispformula = ~ habitat, control = TMB_CTRL),
    "zi_tweedie"     = glmmTMB(fml, data = df, family = tweedie(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "beta"           = glmmTMB(fml, data = df, family = beta_family(link = "logit"), control = TMB_CTRL),
    "beta_disp"      = glmmTMB(fml, data = df, family = beta_family(link = "logit"), dispformula = ~ habitat, control = TMB_CTRL),
    "gaussian"       = glmmTMB(fml, data = df, family = gaussian(link = "identity"), control = TMB_CTRL),
    "lognormal"      = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    "gaussian_log1p" = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    "gaussian_logit" = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    stop("Unknown family: ", family_choice)
  )
  
  mod
}



# Stage 4: Fit + test each response variable ------------------------------


# ----- Main Function: Fit and Test One Model -----
# Fit the configured model for one response column and extract results.
# For a single response variable (`ycol`): subsets `dat` to the relevant columns and renames them to the generic `Y`/`habitat`/`site` used by `fit_with_family()`, drops non-finite/missing rows, fits the model, prints a Type II Anova (car::Anova) and DHARMa residual diagnostics(uniformity, dispersion, outlier tests) to the console for manual inspection, then extracts the pairwise habitat contrast from emmeans and returns it as a tidy one-row data frame.

fit_one <- function(ycol) {
  message("\n=== Fitting model for: ", ycol, " ===")
  
  family_choice <- family_config[[ycol]]
  if (is.null(family_choice)) {
    warning("No family specified for ", ycol, ". Skipping.")
    return(NULL)
  }
  message("Using family: ", family_choice)
  
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(hab_col, site_col, ycol))) %>%
    dplyr::rename(habitat = !!hab_col, site = !!site_col, Y = !!ycol) %>%
    dplyr::mutate(Y = suppressWarnings(as.numeric(Y))) %>%
    dplyr::filter(is.finite(Y), !is.na(habitat), !is.na(site))
  
  if (nrow(df) == 0 || dplyr::n_distinct(df$habitat) < 2L) {
    warning("Skipping ", ycol, ": insufficient rows or only one habitat level.")
    return(NULL)
  }
  
  mod <- try(fit_with_family(df, family_choice), silent = TRUE)
  if (inherits(mod, "try-error")) {
    warning("Model failed for ", ycol, " with family ", family_choice)
    return(NULL)
  }
  
  # Type II Anova on the fixed effect (habitat) -- printed for manual
  # review, not captured into the returned table.
  an <- try(car::Anova(mod), silent = TRUE)
  if (!inherits(an, "try-error")) {
    message("\nType II Anova:")
    print(an)
  }
  
  # DHARMa simulated-residual diagnostics -- also printed for manual
  # review only (uniformity/goodness-of-fit, dispersion, outlier tests).
  dh <- try(DHARMa::simulateResiduals(mod, n = DH_NSIM, plot = FALSE), silent = TRUE)
  if (!inherits(dh, "try-error")) {
    message("\nDHARMa tests:")
    print(DHARMa::testUniformity(dh))
    print(DHARMa::testDispersion(dh))
    print(DHARMa::testOutliers(dh))
  } else {
    message("DHARMa diagnostics failed")
  }
  
  emm_type <- if (family_choice %in% c("lognormal", "gaussian_log1p")) "link" else "response"
  pairs_obj <- try(emmeans::emmeans(mod, specs = pairwise ~ habitat, type = emm_type)$contrasts,
                   silent = TRUE)
  if (inherits(pairs_obj, "try-error") || is.null(pairs_obj)) {
    warning("No contrasts available for ", ycol)
    return(NULL)
  }
  
  pairs_obj %>%
    as.data.frame() %>%
    standardise_contrast_cols() %>%
    dplyr::mutate(response = ycol, family = family_choice) %>%
    dplyr::select(response, family, contrast, estimate, SE, df, t.ratio, p.value)
}

# Stage 5: Run all models and export the pairwise-contrast table ----------

# ----- Run All Models -----

message("Starting final analysis")


# Runs fit_one() over every configured response column and row-binds the
# (possibly NULL-dropping) results into one combined table.
results <- purrr::map_dfr(analysis_cols, fit_one)

# Restore the original family_config ordering (abundance/biomass/richness/
# Shannon x targeted/non-targeted) rather than whatever order map_dfr
# produced, since responses that fail silently drop out of `results`.
results <- results %>%
  dplyr::relocate(response, family, contrast) %>%
  dplyr::arrange(factor(response, levels = analysis_cols))


message("Results summary")
print(results)

write.csv(results, "glmmtmb_pairwise_targeted_nontargeted.csv", row.names = FALSE)
message('\nSaved: glmmtmb_pairwise_targeted_nontargeted.csv')


# Stage 6: Plotting ----------------------------------------------------------------

# Each response variable's model is re-fit (rather than reusing the Stage 4
# `results`/model objects) specifically to derive estimated marginal means.

# ----- Custom Labels for Subplot Titles -----

custom_labels <- list(
  "targeted_abundance"   = "Targeted",
  "nontargeted_abundance"  = "Non-targeted",
  "targeted_biomass"     = "Targeted",
  "nontargeted_biomass"    = "Non-targeted",
  "targeted_richness"    = "Targeted",
  "nontargeted_richness"   = "Non-targeted",
  "targeted_shannon"     = "Targeted",
  "nontargeted_shannon"    = "Non-targeted"
)

# ----- Function: Create Formatted Plot for One Response -----
# Fit one response's model and render it as a styled habitat bar chart.
# Re-fits the response's glmmTMB model (same family/formula as Stage 4), computes habitat-level estimated marginal means via emmeans, and derives a compact letter display marking which habitats differ significantly. Falls back to plain emmeans summary output (with blank CLD letters) if `cld()` fails. The resulting table is normalized to a `response`/`SE` column pair, then plotted as a bar + errorbar chart with letters above each bar.


create_formatted_plot <- function(ycol, show_x_label = FALSE, show_y_label = FALSE,
                                  y_label_text = "", position = c(1, 1), total_rows = 2) {
  
  message("\nCreating formatted plot for: ", ycol)
  
  family_choice <- family_config[[ycol]]
  if (is.null(family_choice)) return(NULL)
  
  # Same subsetting as fit_one(), plus habitat-label formatting ("rocky_reef" -> "Rocky reef") for display on the plot axis/legend.
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(hab_col, site_col, ycol))) %>%
    dplyr::rename(habitat = !!hab_col, site = !!site_col, Y = !!ycol) %>%
    dplyr::mutate(
      Y = suppressWarnings(as.numeric(Y)),
      habitat = gsub("_", " ", as.character(habitat)),
      habitat = dplyr::case_when(
        tolower(habitat) == "rocky reef"     ~ "Rocky reef",
        tolower(habitat) == "indurated reef" ~ "Indurated reef",
        TRUE ~ habitat
      ),
      habitat = factor(habitat, levels = c("Rocky reef", "Indurated reef"))
    ) %>%
    dplyr::filter(is.finite(Y), !is.na(habitat), !is.na(site))
  
  if (nrow(df) == 0 || dplyr::n_distinct(df$habitat) < 2) {
    message(" Insufficient data")
    return(NULL)
  }
  
  mod <- try(fit_with_family(df, family_choice), silent = TRUE)
  if (inherits(mod, "try-error")) {
    message(" Model failed")
    return(NULL)
  }
  
  emm_type <- if (family_choice %in% c("lognormal", "gaussian_log1p")) "link" else "response"
  emm <- try(emmeans::emmeans(mod, specs = ~ habitat, type = emm_type), silent = TRUE)
  if (inherits(emm, "try-error") || is.null(emm)) {
    message("emmeans failed")
    return(NULL)
  }
  
  # Compact letter display: habitats sharing a letter are not significantly different (Tukey-adjusted pairwise comparisons). If cld() itself fails (can happen with some family/link combinations), fall back to a plain emmeans summary table with empty grouping letters so the plot can still be built, just without significance annotations.
  cld_tbl <- try(multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame(),
                 silent = TRUE)
  if (inherits(cld_tbl, "try-error")) {
    cld_tbl <- as.data.frame(summary(emm))
    cld_tbl$.group <- ""
  }
  
  # Normalize whichever point-estimate column name this family/link produced (emmean / rate / prob / other numeric) to a common `response` column so the plotting code below doesn't need to branch on family.
  if ("emmean" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = emmean)
  } else if ("rate" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = rate)
  } else if ("prob" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = prob)
  } else if (!"response" %in% names(cld_tbl)) {
    numeric_cols <- names(cld_tbl)[sapply(cld_tbl, is.numeric)]
    exclude <- c("SE", "df", "asymp.LCL", "asymp.UCL", "lower.CL", "upper.CL")
    candidate <- setdiff(numeric_cols, exclude)[1]
    if (!is.na(candidate)) {
      cld_tbl <- dplyr::rename(cld_tbl, response = !!candidate)
    }
  }
  
  # If cld()/summary() didn't return an SE column directly, back it out from whichever confidence-interval columns are present (assuming a ~95% CI, i.e. +/- 1.96 SE); otherwise default to 0 so errorbars still render (as zero-width) rather than erroring.
  if (!"SE" %in% names(cld_tbl)) {
    if (all(c("lower.CL","upper.CL") %in% names(cld_tbl))) {
      cld_tbl$SE <- (cld_tbl$upper.CL - cld_tbl$lower.CL) / (2 * 1.96)
    } else if (all(c("asymp.LCL","asymp.UCL") %in% names(cld_tbl))) {
      cld_tbl$SE <- (cld_tbl$asymp.UCL - cld_tbl$asymp.LCL) / (2 * 1.96)
    } else {
      cld_tbl$SE <- 0
    }
  }
  

  # y_lab_pos places the significance letter just above each error bar.
  cld_tbl <- cld_tbl %>%
    dplyr::mutate(
      ymin = pmax(response - SE, 0),
      ymax = response + SE,
      y_lab_pos = ymax * 1.05
    )
  
  plot_title <- if (ycol %in% names(custom_labels)) custom_labels[[ycol]] else ycol
  
  hab_colors <- c("Rocky reef" = "#0069a2", "Indurated reef" = "#d55e00")
  
  p <- ggplot2::ggplot(cld_tbl, ggplot2::aes(habitat, response, fill = habitat)) +
    ggplot2::geom_col(width = 0.5, color = NA) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                           width = 0.1, linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .group, y = y_lab_pos),
                       vjust = 0, color = "black", size = 4) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    ggplot2::scale_fill_manual(values = hab_colors, guide = "none") +
    ggplot2::labs(
      title = plot_title,
      y = if (show_y_label) y_label_text else NULL,
      x = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        vjust = 1,
        size = 15,
        margin = margin(b = 5)
      ),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1.7),
      axis.title.y = if (show_y_label) {
        ggplot2::element_text(
          size = 14,
          face = "bold",
          vjust = 0.5,
          angle = 90,
          margin = margin(r = 10)
        )
      } else ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8, color = "black"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = if (show_x_label) {
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size  = 11,
          color = "black"
        )
      } else ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank()
    )
  
  return(p)
}

# ----- Create Panel Layout -----
message("Creating panel plot")


# Order: rows = Abundance, Biomass, Richness, Shannon
#        columns = Targeted, Non-targeted
# Fixed display order for the panel grid (independent of family_config's list order, though in this case they happen to match).
layout_order <- c(
  "targeted_abundance",  "nontargeted_abundance",
  "targeted_biomass",    "nontargeted_biomass",
  "targeted_richness",   "nontargeted_richness",
  "targeted_shannon",    "nontargeted_shannon"
)

# Defensive filter in case family_config is ever edited to drop a response that layout_order still references.
layout_order <- layout_order[layout_order %in% names(family_config)]
n_plots <- length(layout_order)
n_cols  <- 2
n_rows  <- ceiling(n_plots / n_cols)

formatted_plots <- list()

# Build each panel, computing per-panel label visibility (only bottom row gets x-axis tick labels, only left column gets a y-axis title) and the appropriate y-axis title text based on which metric the column name refers to.
for (i in seq_along(layout_order)) {
  ycol <- layout_order[i]
  row <- ceiling(i / n_cols)
  col <- ((i - 1) %% n_cols) + 1
  
  show_x_label <- (row == n_rows)
  show_y_label <- (col == 1)
  
  y_label <- if (grepl("abundance", ycol)) {
    "Abundance"
  } else if (grepl("biomass", ycol)) {
    "Biomass (kg)"
  } else if (grepl("richness", ycol)) {
    "Richness"
  } else if (grepl("shannon", ycol)) {
    "Shannon diversity"
  } else {
    ""
  }
  
  formatted_plots[[i]] <- create_formatted_plot(
    ycol,
    show_x_label = show_x_label,
    show_y_label = show_y_label,
    y_label_text = y_label,
    position     = c(row, col),
    total_rows   = n_rows
  )
}

# ----- Stage 7: Combine and Save Panel Plot -----

if (length(formatted_plots) > 0) {
  panel_plot <- patchwork::wrap_plots(formatted_plots, ncol = n_cols)
  print(panel_plot)
  
  plot_height <- n_rows * 3.5
  plot_width  <- n_cols * 3.5
  
  ggplot2::ggsave(
    "glmmtmb_targeted_nontarget_panel.png",
    panel_plot,
    width  = plot_width,
    height = plot_height,
    dpi    = 300
  )
  
  message("\n Custom panel plot saved: glmmtmb_targeted_nontarget_panel.svg")
  message("  Dimensions: ", plot_width, " x ", plot_height, " inches")
  message("  Layout: ", n_rows, " rows x ", n_cols, " columns")
} else {
  message("No plots generated")
}


# GLMMTMB family biomass
# Carangidae, Lethrinidae, Lutjanidae, Serranidae, Sparidae


# Load & prep data -----------------------------------------------

# Family Biomass configuration

family_config <- list(
  carangidae_biomass   = "gaussian_log1p",
  lethrinidae_biomass  = "gaussian_log1p",
  lutjanidae_biomass   = "gaussian_log1p",
  serranidae_biomass   = "gaussian_log1p",
  sparidae_biomass     = "gaussian_log1p"
)

analysis_cols <- names(family_config)
message("\nAnalyzing ", length(analysis_cols), " family biomass variables")


# Helper functions --------------------------------------------------------

standardise_contrast_cols <- function(tbl) {
  tbl <- as.data.frame(tbl)
  
  if (!("estimate" %in% names(tbl))) {
    if ("ratio" %in% names(tbl))           tbl <- dplyr::rename(tbl, estimate = ratio)
    else if ("odds.ratio" %in% names(tbl)) tbl <- dplyr::rename(tbl, estimate = odds.ratio)
    else if ("response" %in% names(tbl))   tbl <- dplyr::rename(tbl, estimate = response)
    else if ("emmean" %in% names(tbl))     tbl <- dplyr::rename(tbl, estimate = emmean)
    else tbl$estimate <- NA_real_
  }
  
  if (!("t.ratio" %in% names(tbl))) {
    if ("z.ratio" %in% names(tbl))   tbl <- dplyr::rename(tbl, t.ratio = z.ratio)
    else if ("stat" %in% names(tbl)) tbl <- dplyr::rename(tbl, t.ratio = stat)
    else tbl$t.ratio <- NA_real_
  }
  
  if (!("SE" %in% names(tbl)))       tbl$SE      <- NA_real_
  if (!("df" %in% names(tbl)))       tbl$df      <- NA_real_
  if (!("contrast" %in% names(tbl))) tbl$contrast <- "contrast"
  if (!("p.value" %in% names(tbl)))  tbl$p.value <- NA_real_
  
  tbl
}

fit_with_family <- function(df, family_choice) {
  
  if (family_choice %in% c("beta", "beta_disp")) {
    n <- nrow(df)
    if (any(df$Y %in% c(0,1))) {
      df$Y_trans <- (df$Y * (n - 1) + 0.5) / n
    } else {
      df$Y_trans <- pmin(pmax(df$Y, .Machine$double.eps), 1 - .Machine$double.eps)
    }
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "lognormal") {
    df$Y_trans <- log(df$Y)
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "gaussian_log1p") {
    df$Y_trans <- log1p(df$Y)       # log(1 + biomass)
    formula_lhs <- "Y_trans"
    
  } else if (family_choice == "gaussian_logit") {
    df$Y_trans <- qlogis(df$Y)
    formula_lhs <- "Y_trans"
    
  } else {
    formula_lhs <- "Y"
  }
  
  fml <- as.formula(paste(formula_lhs, "~ habitat + (1 | site)"))
  
  mod <- switch(
    family_choice,
    "poisson"        = glmmTMB(fml, data = df, family = poisson(link = "log"),  control = TMB_CTRL),
    "nb2"            = glmmTMB(fml, data = df, family = nbinom2(link = "log"),  control = TMB_CTRL),
    "zi_poisson"     = glmmTMB(fml, data = df, family = poisson(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "zi_nb2"         = glmmTMB(fml, data = df, family = nbinom2(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "gamma"          = glmmTMB(fml, data = df, family = Gamma(link = "log"),    control = TMB_CTRL),
    "tweedie"        = glmmTMB(fml, data = df, family = tweedie(link = "log"),  control = TMB_CTRL),
    "tweedie_disp"   = glmmTMB(fml, data = df, family = tweedie(link = "log"),  dispformula = ~ habitat, control = TMB_CTRL),
    "zi_tweedie"     = glmmTMB(fml, data = df, family = tweedie(link = "log"),  ziformula = ~1, control = TMB_CTRL),
    "beta"           = glmmTMB(fml, data = df, family = beta_family(link = "logit"), control = TMB_CTRL),
    "beta_disp"      = glmmTMB(fml, data = df, family = beta_family(link = "logit"), dispformula = ~ habitat, control = TMB_CTRL),
    "gaussian"       = glmmTMB(fml, data = df, family = gaussian(link = "identity"), control = TMB_CTRL),
    "lognormal"      = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    "gaussian_log1p" = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    "gaussian_logit" = glmmTMB(fml, data = df, family = gaussian(),            control = TMB_CTRL),
    stop("Unknown family: ", family_choice)
  )
  
  mod
}

# Transformed-response families stay on "link" scale
get_emm_type <- function(family_choice) {
  if (family_choice %in% c("lognormal", "gaussian_log1p", "gaussian_logit")) {
    return("link")
  }
  return("response")
}


# Fit + test each response variable ---------

fit_one <- function(ycol) {
  message("\n=== Fitting model for: ", ycol, " ===")
  
  family_choice <- family_config[[ycol]]
  if (is.null(family_choice)) {
    warning("No family specified for ", ycol)
    return(NULL)
  }
  
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(hab_col, site_col, ycol))) %>%
    dplyr::rename(habitat = !!hab_col, site = !!site_col, Y = !!ycol) %>%
    dplyr::mutate(Y = suppressWarnings(as.numeric(Y))) %>%
    dplyr::filter(is.finite(Y), !is.na(habitat), !is.na(site))
  
  if (nrow(df) == 0 || dplyr::n_distinct(df$habitat) < 2L) {
    warning("Skipping ", ycol, ": insufficient rows or only one habitat level.")
    return(NULL)
  }
  
  mod <- try(fit_with_family(df, family_choice), silent = TRUE)
  if (inherits(mod, "try-error")) {
    warning("Model failed for ", ycol, " with family ", family_choice)
    return(NULL)
  }
  
  an <- try(car::Anova(mod), silent = TRUE)
  if (!inherits(an, "try-error")) {
    message("\nType II Anova:")
    print(an)
  }
  
  dh <- try(DHARMa::simulateResiduals(mod, n = DH_NSIM, plot = FALSE), silent = TRUE)
  if (!inherits(dh, "try-error")) {
    message("\nDHARMa tests:")
    print(DHARMa::testUniformity(dh))
    print(DHARMa::testDispersion(dh))
    print(DHARMa::testOutliers(dh))
  }
  
  emm_type  <- get_emm_type(family_choice)
  pairs_obj <- try(
    emmeans::emmeans(mod, specs = pairwise ~ habitat, type = emm_type)$contrasts,
    silent = TRUE
  )
  if (inherits(pairs_obj, "try-error") || is.null(pairs_obj)) {
    warning("No contrasts available for ", ycol)
    return(NULL)
  }
  
  pairs_obj %>%
    as.data.frame() %>%
    standardise_contrast_cols() %>%
    dplyr::mutate(response = ycol, family = family_choice) %>%
    dplyr::select(response, family, contrast, estimate, SE, df, t.ratio, p.value)
}

message("Starting family biomass analysis")

results_family <- purrr::map_dfr(analysis_cols, fit_one) %>%
  dplyr::relocate(response, family, contrast) %>%
  dplyr::arrange(factor(response, levels = analysis_cols))

message("Summary")
print(results_family)

write.csv(results_family, "glmmtmb_pairwise_family_biomass.csv", row.names = FALSE)
message('\nSaved: glmmtmb_pairwise_family_biomass.csv')

# Panel plot

custom_labels <- list(
  carangidae_biomass  = "Carangidae",
  lethrinidae_biomass = "Lethrinidae",
  lutjanidae_biomass  = "Lutjanidae",
  serranidae_biomass  = "Serranidae",
  sparidae_biomass    = "Sparidae"
)

create_formatted_plot <- function(ycol, show_x_label = FALSE, show_y_label = FALSE,
                                  y_label_text = "") {
  message("\n=== Creating formatted plot for: ", ycol, " ===")
  
  family_choice <- family_config[[ycol]]
  if (is.null(family_choice)) return(NULL)
  
  df <- dat %>%
    dplyr::select(dplyr::all_of(c(hab_col, site_col, ycol))) %>%
    dplyr::rename(habitat = !!hab_col, site = !!site_col, Y = !!ycol) %>%
    dplyr::mutate(
      Y = suppressWarnings(as.numeric(Y)),
      habitat = gsub("_", " ", as.character(habitat)),
      habitat = dplyr::case_when(
        tolower(habitat) == "rocky reef"     ~ "Rocky reef",
        tolower(habitat) == "indurated reef" ~ "Indurated reef",
        TRUE ~ habitat
      ),
      habitat = factor(habitat, levels = c("Rocky reef", "Indurated reef"))
    ) %>%
    dplyr::filter(is.finite(Y), !is.na(habitat), !is.na(site))
  
  if (nrow(df) == 0 || dplyr::n_distinct(df$habitat) < 2) {
    message(" Insufficient data")
    return(NULL)
  }
  
  mod <- try(fit_with_family(df, family_choice), silent = TRUE)
  if (inherits(mod, "try-error")) return(NULL)
  
  emm_type <- get_emm_type(family_choice)
  emm <- try(emmeans::emmeans(mod, specs = ~ habitat, type = emm_type), silent = TRUE)
  if (inherits(emm, "try-error") || is.null(emm)) return(NULL)
  
  cld_tbl <- try(multcomp::cld(emm, Letters = letters, adjust = "tukey") %>% as.data.frame(),
                 silent = TRUE)
  if (inherits(cld_tbl, "try-error")) {
    cld_tbl <- as.data.frame(summary(emm))
    cld_tbl$.group <- ""
  }
  
  if ("emmean" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = emmean)
  } else if ("rate" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = rate)
  } else if ("prob" %in% names(cld_tbl)) {
    cld_tbl <- dplyr::rename(cld_tbl, response = prob)
  } else if (!"response" %in% names(cld_tbl)) {
    numeric_cols <- names(cld_tbl)[sapply(cld_tbl, is.numeric)]
    exclude      <- c("SE","df","asymp.LCL","asymp.UCL","lower.CL","upper.CL")
    candidate    <- setdiff(numeric_cols, exclude)[1]
    if (!is.na(candidate)) {
      cld_tbl <- dplyr::rename(cld_tbl, response = !!candidate)
    }
  }
  
  if (!"SE" %in% names(cld_tbl)) {
    if (all(c("lower.CL","upper.CL") %in% names(cld_tbl))) {
      cld_tbl$SE <- (cld_tbl$upper.CL - cld_tbl$lower.CL) / (2 * 1.96)
    } else if (all(c("asymp.LCL","asymp.UCL") %in% names(cld_tbl))) {
      cld_tbl$SE <- (cld_tbl$asymp.UCL - cld_tbl$asymp.LCL) / (2 * 1.96)
    } else {
      cld_tbl$SE <- 0
    }
  }
  
  cld_tbl <- cld_tbl %>%
    dplyr::mutate(
      ymin      = pmax(response - SE, 0),
      ymax      = response + SE,
      y_lab_pos = ymax * 1.05
    )
  
  plot_title <- custom_labels[[ycol]]
  
  hab_colors <- c("Rocky reef" = "#0069a2", "Indurated reef" = "#d55e00")
  
  ggplot(cld_tbl, aes(habitat, response, fill = habitat)) +
    geom_col(width = 0.5, color = NA) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.1, linewidth = 0.8) +
    geom_text(aes(label = .group, y = y_lab_pos),
              vjust = 0, color = "black", size = 4) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_fill_manual(values = hab_colors, guide = "none") +
    labs(
      title = plot_title,
      y     = if (show_y_label) y_label_text else NULL,
      x     = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        face = "bold", hjust = 0.5, vjust = 1,
        size = 15, margin = margin(b = 5)
      ),
      panel.border = element_rect(color = "black", linewidth = 1.7),
      axis.title.y = if (show_y_label) {
        element_text(size = 14, face = "bold", vjust = 0.5,
                     angle = 90, margin = margin(r = 10))
      } else element_blank(),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x  = if (show_x_label) {
        element_text(angle = 45, hjust = 1, vjust = 1,
                     size = 11, color = "black")
      } else element_blank(),
      axis.ticks   = element_blank()
    )
}

# Build 1 × 5 panel
layout_order <- analysis_cols
n_plots <- length(layout_order)
n_cols  <- 5
n_rows  <- ceiling(n_plots / n_cols)

formatted_plots_family <- vector("list", n_plots)

for (i in seq_along(layout_order)) {
  ycol <- layout_order[i]
  row  <- ceiling(i / n_cols)
  col  <- ((i - 1) %% n_cols) + 1
  
  show_x_label <- (row == n_rows)
  show_y_label <- (col == 1)
  
  formatted_plots_family[[i]] <- create_formatted_plot(
    ycol,
    show_x_label = show_x_label,
    show_y_label = show_y_label,
    y_label_text = "Biomass (kg)"
  )
}

if (length(formatted_plots_family) > 0) {
  family_panel <- patchwork::wrap_plots(formatted_plots_family, ncol = n_cols)
  print(family_panel)
  
  plot_height <- n_rows * 3.5
  plot_width  <- n_cols * 3.5
  
  ggsave(
    "glmmtmb_family_biomass_panel.pdf",
    family_panel,
    width  = plot_width,
    height = plot_height,
    dpi    = 300
  )
  
  message("\n✓ Family BIOMASS panel plot saved: glmmtmb_family_biomass_panel.svg")
  message("  Dimensions: ", plot_width, " x ", plot_height, " inches")
  message("  Layout: ", n_rows, " rows x ", n_cols, " columns")
} else {
  message("No family Biomass plots generated")
}

