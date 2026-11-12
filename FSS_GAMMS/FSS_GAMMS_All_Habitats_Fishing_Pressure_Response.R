
# LARSEN ET AL. 2026 | FSS-GAMM: WHAT PREDICTS FISHING PRESSURE?
# Standalone complementary analysis; does NOT change existing fish-response models.
# Response: New_Log_RS (transformed fishing-pressure index).
# Candidate ecological predictors: five 2D/3D geomorphological/seascape metrics.
# Baseline: spatial smooths + habitat random intercept, as in existing workflow.
# Model ranking: ML AICc; selected model refit by REML for diagnostics/plots.
# All candidate combinations have 0, 1 or 2 ecological predictors.
# Requires: FSSgam, mgcv, MuMIn, ggplot2, cowplot; DHARMa for diagnostics.


# ---- 0. EDIT THESE SETTINGS --------------------------------------------------
PROJECT_DIR <- "input" # Where user wants project directory
DATA_FILE <- file.path(PROJECT_DIR, "input") # Add master csv data file

# Dedicated folder; does not overwrite combined_reef_gamms_REVISED analyses.
OUT_DIR <- file.path(PROJECT_DIR, "Fishing_Pressure_GAMMS",
                     "Outputs")

RESPONSE <- "New_Log_RS"
CONTINUOUS_PREDICTORS <- c(
  "distance_to_estuary",
  "depth",
  "relief_max_500m",
  "slope_mean_500m",
  "prox_index_rocky_reef_500m_log"
)

# Keep habitat as a RANDOM intercept, matching your latest working analysis.
INCLUDE_HABITAT_RANDOM_EFFECT <- TRUE
INCLUDE_SITE_RANDOM_EFFECT <- FALSE
SITE_COLUMN <- "site"     # Only needed if INCLUDE_SITE_RANDOM_EFFECT = TRUE.

# Gaussian is a defensible starting family for a transformed continuous index.
# Set to "beta" only if EVERY observed value is strictly between 0 and 1;
# beta and Gaussian are compared on exactly the same response scale if used.
RESPONSE_FAMILY <- "gaussian"  # "gaussian" or "beta"

MAX_PREDICTORS <- 2L
SMOOTH_K <- 3L
SPATIAL_K <- 3L
CORRELATION_CUTOFF <- 0.70
MIN_ROWS <- 25L
SEED <- 123
RUN_DHARMA <- TRUE
RUN_SPATIAL_TEST <- TRUE
N_SIMULATIONS <- 1000L
CI_LEVEL <- 0.95

PREDICTOR_HEADINGS <- c(
  distance_to_estuary = "Distance to estuary (m)",
  depth = "Depth (m)",
  relief_max_500m = "Maximum relief (500 m)",
  slope_mean_500m = "Mean slope (500 m)",
  prox_index_rocky_reef_500m_log = "Rocky reef proximity index"
)

# ---- 1. PACKAGES / DATA / AUDIT ---------------------------------------------
NEEDED <- c("FSSgam", "mgcv", "MuMIn", "ggplot2", "cowplot")
missing_pkgs <- NEEDED[!vapply(NEEDED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages: ", paste(missing_pkgs, collapse = ", "))
}
if (RUN_DHARMA && !requireNamespace("DHARMa", quietly = TRUE)) {
  stop("Install DHARMa or set RUN_DHARMA <- FALSE")
}
if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)
if (MAX_PREDICTORS < 1L || MAX_PREDICTORS > length(CONTINUOUS_PREDICTORS)) {
  stop("Choose MAX_PREDICTORS between 1 and the number of candidate predictors")
}
set.seed(SEED)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "models"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "diagnostics"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)

pretty_predictor <- function(x) {
  if (x %in% names(PREDICTOR_HEADINGS)) unname(PREDICTOR_HEADINGS[[x]])
  else gsub("_", " ", x)
}
write_table <- function(x, path) utils::write.csv(x, path, row.names = FALSE)

raw <- utils::read.csv(DATA_FILE, check.names = TRUE, stringsAsFactors = FALSE)
required <- c(RESPONSE, "longitude", "latitude", CONTINUOUS_PREDICTORS,
              if (INCLUDE_HABITAT_RANDOM_EFFECT) "habitat",
              if (INCLUDE_SITE_RANDOM_EFFECT) SITE_COLUMN)
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) stop("Required columns missing: ", paste(missing_cols, collapse = ", "))

for (v in unique(c(RESPONSE, "longitude", "latitude", CONTINUOUS_PREDICTORS))) {
  raw[[v]] <- suppressWarnings(as.numeric(raw[[v]]))
  raw[[v]][!is.finite(raw[[v]])] <- NA_real_
}
if (INCLUDE_HABITAT_RANDOM_EFFECT) raw$habitat <- factor(raw$habitat)
if (INCLUDE_SITE_RANDOM_EFFECT) raw[[SITE_COLUMN]] <- factor(raw[[SITE_COLUMN]])

# Use the same complete-case observations for every candidate model.
# Candidate predictors with excessive missingness or insufficient variation
# are dropped BEFORE complete-case selection, consistent with the fish script.
predictors <- CONTINUOUS_PREDICTORS[
  vapply(CONTINUOUS_PREDICTORS, function(x) {
    mean(is.na(raw[[x]])) <= 0.25 &&
      length(unique(stats::na.omit(raw[[x]]))) >= 4L
  }, logical(1))
]
if (!length(predictors)) stop("No eligible continuous predictors")
excluded <- setdiff(CONTINUOUS_PREDICTORS, predictors)
if (length(excluded)) warning("Excluded unavailable/low-variation predictors: ",
                             paste(excluded, collapse = ", "))

cols <- unique(c(RESPONSE, "longitude", "latitude", predictors,
                 if (INCLUDE_HABITAT_RANDOM_EFFECT) "habitat",
                 if (INCLUDE_SITE_RANDOM_EFFECT) SITE_COLUMN))
# Audit-only site names must not determine the complete-case sample size.
complete_rows <- which(stats::complete.cases(raw[, cols, drop = FALSE]))
dd <- droplevels(raw[complete_rows, cols, drop = FALSE])
rownames(dd) <- NULL
if (nrow(dd) < MIN_ROWS) stop("Only ", nrow(dd), " complete cases; required ", MIN_ROWS)
if (length(unique(dd[[RESPONSE]])) < 3L) stop("Fishing pressure response has too few unique values")
if (length(unique(dd$longitude)) < 4L || length(unique(dd$latitude)) < 4L) {
  stop("Insufficient coordinate variation for spatial smooths")
}
if (INCLUDE_HABITAT_RANDOM_EFFECT && nlevels(dd$habitat) < 2L) {
  stop("Habitat random intercept needs at least two observed levels")
}
if (INCLUDE_SITE_RANDOM_EFFECT && nlevels(dd[[SITE_COLUMN]]) < 3L) {
  stop("Site random effect needs at least three observed levels")
}
predictors <- predictors[vapply(predictors, function(x) length(unique(dd[[x]])) >= 4L,
                                logical(1))]
if (!length(predictors)) stop("No varying candidate predictors after complete-case selection")

# Alias the response, so the formula works unchanged in FSSgam.
dd$response <- dd[[RESPONSE]]
if (RESPONSE_FAMILY == "gaussian") {
  family_object <- stats::gaussian()
} else if (RESPONSE_FAMILY == "beta") {
  if (any(dd$response <= 0 | dd$response >= 1)) {
    stop("Beta regression requires all fishing index values strictly inside (0, 1).")
  }
  family_object <- mgcv::betar(link = "logit")
} else {
  stop("RESPONSE_FAMILY must be 'gaussian' or 'beta'")
}

# Check possible pseudoreplication: is the exposure index identical at all
# replicate positions for a given reef/site? If so, use site-level analysis
# or an explicit site correlation structure rather than claiming independent n.
site_audit <- NULL
if (SITE_COLUMN %in% names(raw)) {
  site_names <- raw[[SITE_COLUMN]][complete_rows]
  site_names <- ifelse(is.na(site_names), "<missing site>", as.character(site_names))
  spl <- split(dd$response, site_names)
  site_audit <- data.frame(
    site = names(spl),
    n_rows = vapply(spl, length, integer(1)),
    n_unique_pressure = vapply(spl, function(x) length(unique(x)), integer(1)),
    stringsAsFactors = FALSE
  )
  write_table(site_audit, file.path(OUT_DIR, "SITE_REPLICATION_AUDIT.csv"))
  repeated_constant <- subset(site_audit, n_rows > 1 & n_unique_pressure == 1)
  if (nrow(repeated_constant)) warning(
    nrow(repeated_constant), " sites have identical pressure in repeated records. ",
    "Review independence/pseudoreplication before interpreting p-values.")
}

write_table(data.frame(
  response = RESPONSE, distribution = RESPONSE_FAMILY,
  n_initial = nrow(raw), n_complete = nrow(dd),
  min_response = min(dd$response), max_response = max(dd$response),
  zeros = sum(dd$response == 0),
  predictors = paste(predictors, collapse = "; "),
  habitat_random = INCLUDE_HABITAT_RANDOM_EFFECT,
  site_random = INCLUDE_SITE_RANDOM_EFFECT,
  stringsAsFactors = FALSE
), file.path(OUT_DIR, "INPUT_AND_MODEL_AUDIT.csv"))
write_table(data.frame(predictor = CONTINUOUS_PREDICTORS,
                       heading = vapply(CONTINUOUS_PREDICTORS, pretty_predictor, character(1)),
                       eligible = CONTINUOUS_PREDICTORS %in% predictors),
            file.path(OUT_DIR, "PREDICTOR_HEADINGS.csv"))

if (length(predictors) > 1L) {
  cor_mat <- stats::cor(dd[, predictors, drop = FALSE], use = "complete.obs")
  utils::write.csv(cor_mat, file.path(OUT_DIR, "PREDICTOR_CORRELATIONS.csv"))
}

# ---- 2. SPATIAL + RANDOM INTERCEPT BASELINE ---------------------------------
base_terms <- c(sprintf("s(longitude, k=%d)", SPATIAL_K),
                sprintf("s(latitude, k=%d)", SPATIAL_K))
if (INCLUDE_HABITAT_RANDOM_EFFECT) base_terms <- c(base_terms, "s(habitat, bs='re')")
if (INCLUDE_SITE_RANDOM_EFFECT) {
  base_terms <- c(base_terms, sprintf("s(%s, bs='re')", SITE_COLUMN))
}
BASELINE_TERMS <- paste(base_terms, collapse = " + ")
null_formula <- stats::as.formula(paste("response ~", BASELINE_TERMS))
null_fit <- mgcv::gam(null_formula, family = family_object, data = dd,
                      method = "ML", na.action = stats::na.fail)
if (!isTRUE(null_fit$converged)) stop("Spatial/random-intercept baseline failed to converge")
# FSSgam internally updates this call; preserve the concrete family object.
null_fit$call$family <- family_object
saveRDS(null_fit, file.path(OUT_DIR, "models", "SPATIAL_HABITAT_BASELINE_ML.rds"))

# ---- 3. GENERATE AND FIT FULL SUBSETS ---------------------------------------
generate_fun <- if ("generate_model_set" %in% getNamespaceExports("FSSgam")) {
  getExportedValue("FSSgam", "generate_model_set")
} else getExportedValue("FSSgam", "generate.model.set")
fit_fun <- if ("fit_model_set" %in% getNamespaceExports("FSSgam")) {
  getExportedValue("FSSgam", "fit_model_set")
} else getExportedValue("FSSgam", "fit.model.set")

model_set <- generate_fun(
  use.dat = dd,
  test.fit = null_fit,
  pred.vars.cont = predictors,
  pred.vars.fact = NA,
  null.terms = BASELINE_TERMS,
  max.predictors = min(MAX_PREDICTORS, length(predictors)),
  cov.cutoff = CORRELATION_CUTOFF,
  k = SMOOTH_K,
  factor.smooth.interactions = FALSE,
  factor.factor.interactions = FALSE,
  smooth.smooth.interactions = FALSE
)

n_candidates <- length(model_set$mod.formula)
if (!n_candidates) stop("FSSgam did not generate any candidate formulas")
message("Fitting ", n_candidates, " FSS-GAMM candidates for fishing pressure")
fit_results <- fit_fun(model_set, save.model.fits = TRUE, parallel = FALSE,
                       max.models = max(200L, n_candidates + 10L), VI.mods = "all")
if (!is.null(fit_results$mod.data.out)) {
  write_table(fit_results$mod.data.out, file.path(OUT_DIR, "FSSgam_NATIVE_MODEL_TABLE.csv"))
}
if (!is.null(fit_results$variable.importance)) {
  saveRDS(fit_results$variable.importance,
          file.path(OUT_DIR, "FSSgam_NATIVE_IMPORTANCE.rds"))
}
if (length(fit_results$failed.models)) {
  writeLines(capture.output(str(fit_results$failed.models)),
             file.path(OUT_DIR, "FAILED_CANDIDATE_MODELS.txt"))
}

fits <- c(list(SPATIAL_BASELINE = null_fit), fit_results$success.models)
fits <- fits[vapply(fits, inherits, logical(1), what = "gam")]
if (!length(fits)) stop("No usable fitted candidate models")

# Guard against accidental loss of the mandatory habitat random effect.
if (INCLUDE_HABITAT_RANDOM_EFFECT) {
  re_present <- vapply(fits, function(m) {
    any(vapply(m$smooth, function(s) {
      inherits(s, "random.effect") && ("habitat" %in% s$term)
    }, logical(1)))
  }, logical(1))
  if (!all(re_present)) stop("Some candidate fits lost the habitat random intercept")
}

# AICc recomputed from the actual fitted ML models; same observations throughout.
candidate_table <- do.call(rbind, lapply(seq_along(fits), function(i) {
  m <- fits[[i]]
  data.frame(model_id = names(fits)[i],
             AICc = tryCatch(as.numeric(MuMIn::AICc(m)), error = function(e) NA_real_),
             AIC = as.numeric(stats::AIC(m)),
             n = stats::nobs(m),
             n_parameters = as.numeric(attr(stats::logLik(m), "df")),
             converged = isTRUE(m$converged),
             formula = paste(deparse(stats::formula(m)), collapse = " "),
             deviance_explained_ML = summary(m)$dev.expl,
             stringsAsFactors = FALSE)
}))
candidate_table <- subset(candidate_table, converged & is.finite(AICc))
if (!nrow(candidate_table)) stop("No converged models with valid AICc")
if (length(unique(candidate_table$n)) != 1L || unique(candidate_table$n) != nrow(dd)) {
  stop("Candidate models used inconsistent observation counts")
}
# Keep one model per identical formula (the null can also be in FSSgam output).
candidate_table <- candidate_table[!duplicated(candidate_table$formula), , drop = FALSE]
candidate_table <- candidate_table[order(candidate_table$AICc,
                                         candidate_table$n_parameters), , drop = FALSE]
candidate_table$delta_AICc <- candidate_table$AICc - min(candidate_table$AICc)
w <- exp(-0.5 * candidate_table$delta_AICc)
candidate_table$weight_AICc <- w / sum(w)
candidate_table$within_delta2 <- candidate_table$delta_AICc <= 2
candidate_table$is_best <- seq_len(nrow(candidate_table)) == 1L
options_within_2 <- which(candidate_table$within_delta2)
parsimonious_index <- options_within_2[which.min(candidate_table$n_parameters[options_within_2])]
candidate_table$parsimonious_within_2 <- seq_len(nrow(candidate_table)) == parsimonious_index
write_table(candidate_table, file.path(OUT_DIR, "ALL_CANDIDATES_RANKED_AICc.csv"))
write_table(subset(candidate_table, within_delta2),
            file.path(OUT_DIR, "TOP_MODELS_DELTA_AICc_LE_2.csv"))

# ---- 4. VARIABLE IMPORTANCE / BEST-FIT SUMMARY -----------------------------
# Importance = sum of normalised AICc weights for all admissible models
# that include each predictor (including baseline model in denominator).
incidence <- vapply(seq_len(nrow(candidate_table)), function(i) {
  m <- fits[[as.character(candidate_table$model_id[i])]]
  predictors %in% all.vars(stats::formula(m))
}, logical(length(predictors)))
incidence <- matrix(incidence, nrow = length(predictors),
                    ncol = nrow(candidate_table))
rownames(incidence) <- predictors
importance <- as.numeric(incidence %*% candidate_table$weight_AICc)
best_model_id <- as.character(candidate_table$model_id[1])
best_ml <- fits[[best_model_id]]
best_predictors <- intersect(predictors, all.vars(stats::formula(best_ml)))
importance_table <- data.frame(
  predictor = predictors,
  predictor_heading = vapply(predictors, pretty_predictor, character(1)),
  AICc_importance = importance,
  n_candidate_models_including = rowSums(incidence),
  n_candidate_models_total = nrow(candidate_table),
  in_best_model = predictors %in% best_predictors,
  stringsAsFactors = FALSE
)
importance_table <- importance_table[order(-importance_table$AICc_importance), ]
importance_table$rank <- seq_len(nrow(importance_table))
write_table(importance_table, file.path(OUT_DIR, "VARIABLE_IMPORTANCE_AICc_LABELLED.csv"))
write_table(as.data.frame(as.list(stats::setNames(importance_table$AICc_importance,
                                                  importance_table$predictor_heading)),
                          check.names = FALSE),
            file.path(OUT_DIR, "VARIABLE_IMPORTANCE_WIDE_LABELLED.csv"))

saveRDS(best_ml, file.path(OUT_DIR, "models", "BEST_ML.rds"))
best_reml <- tryCatch(stats::update(best_ml, method = "REML"), error = function(e) {
  warning("REML refit failed, using ML best fit: ", conditionMessage(e))
  best_ml
})
if (!isTRUE(best_reml$converged)) {
  warning("REML refit did not converge; using ML fit for plots and diagnostics")
  best_reml <- best_ml
}
saveRDS(best_reml, file.path(OUT_DIR, "models", "BEST_REML.rds"))
writeLines(c(
  "Best selected ML model for fishing-pressure response:",
  candidate_table$formula[1],
  paste("Family:", RESPONSE_FAMILY),
  paste("n:", nrow(dd)),
  paste("AICc:", signif(candidate_table$AICc[1], 7)),
  paste("Other model(s) within delta <= 2:", sum(candidate_table$within_delta2) - 1L),
  "", capture.output(summary(best_reml))
), file.path(OUT_DIR, "BEST_MODEL_SUMMARY.txt"))

# ---- 5. MODEL DIAGNOSTICS ---------------------------------------------------
tryCatch(utils::write.csv(mgcv::concurvity(best_reml, full = TRUE),
                          file.path(OUT_DIR, "diagnostics", "CONCURVITY.csv")),
         error = function(e) warning("Concurvity failed: ", conditionMessage(e)))
tryCatch(utils::write.csv(mgcv::k.check(best_reml),
                          file.path(OUT_DIR, "diagnostics", "K_CHECK.csv")),
         error = function(e) warning("k-check failed: ", conditionMessage(e)))

png(file.path(OUT_DIR, "diagnostics", "GAM_CHECK.png"),
    width = 1350, height = 1000, res = 150)
try(mgcv::gam.check(best_reml), silent = TRUE)
grDevices::dev.off()

if (RUN_DHARMA) {
  set.seed(SEED)
  sim <- tryCatch(DHARMa::simulateResiduals(best_reml, n = N_SIMULATIONS),
                  error = function(e) {warning(conditionMessage(e)); NULL})
  if (!is.null(sim)) {
    png(file.path(OUT_DIR, "diagnostics", "DHARMA_OVERVIEW.png"),
        width = 1300, height = 950, res = 150)
    plot(sim)
    grDevices::dev.off()
    fmt_test <- function(f) tryCatch(capture.output(print(f(sim))),
                                     error = function(e) paste("ERROR:", conditionMessage(e)))
    results <- c("Uniformity", fmt_test(DHARMa::testUniformity),
                 "", "Dispersion", fmt_test(DHARMa::testDispersion),
                 "", "Outliers", fmt_test(DHARMa::testOutliers),
                 "", "Zero inflation", fmt_test(DHARMa::testZeroInflation))
    if (RUN_SPATIAL_TEST && requireNamespace("sf", quietly = TRUE)) {
      coords <- best_reml$model[, c("longitude", "latitude"), drop = FALSE]
      if (!anyDuplicated(paste(coords$longitude, coords$latitude)) &&
          all(stats::complete.cases(coords))) {
        spatial_test <- tryCatch({
          pts <- sf::st_as_sf(coords, coords = c("longitude", "latitude"), crs = 4326)
          xy <- sf::st_coordinates(sf::st_transform(pts, 32756))
          capture.output(print(DHARMa::testSpatialAutocorrelation(
            sim, x = xy[, 1], y = xy[, 2])))
        }, error = function(e) paste("ERROR:", conditionMessage(e)))
        results <- c(results, "", "Spatial autocorrelation (exploratory)", spatial_test)
      } else {
        results <- c(results, "", "Spatial Moran's I skipped: duplicate/missing coordinates")
      }
    }
    writeLines(results, file.path(OUT_DIR, "diagnostics", "DHARMA_TESTS.txt"))
  }
}

# ---- 6. PLOTS: ONE FOR EACH SELECTED EFFECT; SE WITH MEAN -------------------
# type='iterms', iterms.type=1 is mgcv's seWithMean=TRUE counterpart for
# smooth effects. The curve is the centred partial term; uncertainty is widened
# to include estimated overall-mean uncertainty. Confidence bands are
# approximate pointwise 95% intervals on the linear predictor scale.
plot_data <- best_reml$model
make_ref <- function(model, n = 150L) {
  dd0 <- model$model
  nd <- dd0[rep(1L, n), , drop = FALSE]
  for (nm in names(dd0)) {
    if (is.factor(dd0[[nm]])) {
      mode_level <- names(sort(table(dd0[[nm]]), decreasing = TRUE))[1L]
      nd[[nm]] <- factor(rep(mode_level, n), levels = levels(dd0[[nm]]))
    } else if (is.numeric(dd0[[nm]])) {
      nd[[nm]] <- rep(stats::median(dd0[[nm]], na.rm = TRUE), n)
    }
  }
  nd
}

smooth_plots <- list()
for (sm in best_reml$smooth) {
  if (inherits(sm, "random.effect")) next
  if (length(sm$term) != 1L) next
  v <- sm$term[[1L]]
  if (!(v %in% predictors) || !is.numeric(plot_data[[v]])) next
  nd <- make_ref(best_reml)
  nd[[v]] <- seq(min(plot_data[[v]], na.rm = TRUE),
                 max(plot_data[[v]], na.rm = TRUE), length.out = nrow(nd))
  pr <- stats::predict(best_reml, newdata = nd, type = "iterms",
                       terms = sm$label, se.fit = TRUE, iterms.type = 1)
  if (!(sm$label %in% colnames(pr$fit))) next
  z <- data.frame(x = nd[[v]], fit = as.numeric(pr$fit[, sm$label]),
                  se = as.numeric(pr$se.fit[, sm$label]))
  delta <- stats::qnorm(1 - (1 - CI_LEVEL)/2) * z$se
  z$lower <- z$fit - delta
  z$upper <- z$fit + delta
  p <- ggplot2::ggplot(z, ggplot2::aes(x, fit)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#ad1735", alpha = 0.22) +
    ggplot2::geom_line(colour = "#ad1735", linewidth = 0.85) +
    ggplot2::geom_rug(data = data.frame(x = plot_data[[v]]),
                      ggplot2::aes(x = x), inherit.aes = FALSE,
                      colour = "#ad1735", alpha = 0.17) +
    ggplot2::labs(x = pretty_predictor(v), y = "Partial effect on fishing pressure (link scale)") +
    ggplot2::theme_classic(base_size = 11)
  smooth_plots[[v]] <- p
  ggplot2::ggsave(file.path(OUT_DIR, "figures", paste0("PARTIAL_EFFECT_", v, ".png")),
                  p, width = 5.5, height = 4.2, dpi = 300)
}
if (length(smooth_plots)) {
  combined <- cowplot::plot_grid(plotlist = smooth_plots, ncol = min(2L, length(smooth_plots)),
                                 align = "hv", axis = "tblr")
  ggplot2::ggsave(file.path(OUT_DIR, "figures", "FISHING_PRESSURE_SELECTED_PARTIAL_EFFECTS.pdf"),
                  combined, width = 10.5, height = 5.1 * ceiling(length(smooth_plots)/2))
  ggplot2::ggsave(file.path(OUT_DIR, "figures", "FISHING_PRESSURE_SELECTED_PARTIAL_EFFECTS.png"),
                  combined, width = 10.5, height = 5.1 * ceiling(length(smooth_plots)/2), dpi = 300)
} else {
  writeLines("Spatial/habitat-only baseline selected; no ecological partial effects.",
             file.path(OUT_DIR, "figures", "NO_ECOLOGICAL_PARTIAL_EFFECTS.txt"))
}

# Observed vs fitted values (Gaussian on transformed index: both in index units).
observed_fitted <- data.frame(observed = as.numeric(best_reml$y),
                              fitted = as.numeric(stats::fitted(best_reml)))
p_of <- ggplot2::ggplot(observed_fitted, ggplot2::aes(observed, fitted)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey65") +
  ggplot2::geom_point(alpha = .65, colour = "#166785") +
  ggplot2::labs(x = "Observed fishing pressure index",
                y = "Fitted fishing pressure index",
                subtitle = sprintf("Family: %s | n = %d | deviance explained = %.1f%%",
                                   RESPONSE_FAMILY, nrow(dd),
                                   100 * summary(best_reml)$dev.expl)) +
  ggplot2::theme_classic(base_size = 11)
ggplot2::ggsave(file.path(OUT_DIR, "figures", "OBSERVED_VS_FITTED.pdf"),
                p_of, width = 6.2, height = 5.2)

# Show raw relationships for ALL FIVE predictors, not just the selected two.
# These are descriptive scatterplots, not conditional GAM partial effects.
raw_plots <- lapply(predictors, function(v) {
  ggplot2::ggplot(dd, ggplot2::aes(x = .data[[v]], y = response)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.7, colour = "#166785") +
    ggplot2::labs(x = pretty_predictor(v), y = "Fishing pressure index") +
    ggplot2::theme_classic(base_size = 10)
})
raw_fig <- cowplot::plot_grid(plotlist = raw_plots, ncol = 2)
ggplot2::ggsave(file.path(OUT_DIR, "figures", "RAW_RELATIONSHIPS_ALL_FIVE_PREDICTORS.pdf"),
                raw_fig, width = 11, height = 4.5 * ceiling(length(raw_plots)/2))

writeLines(capture.output(utils::sessionInfo()),
           file.path(OUT_DIR, "R_SESSION_INFO.txt"))
message("\nDONE: fishing-pressure response FSS-GAMM saved in: ", OUT_DIR)
message("Best model: ", candidate_table$formula[1])
message("Candidate set: ", nrow(candidate_table), " converged unique models")
message("See AICc ranking, labelled importance, diagnostics and figures.")
