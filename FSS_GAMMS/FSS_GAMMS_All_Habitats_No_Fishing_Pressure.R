
# LARSEN ET AL. - FULL-SUBSETS GAM ANALYSIS. FISHING PRESSURE REMOVED
# Rocky + indurated reefs | spatial baseline | habitat RANDOM effect
# Model selection: ML + AICc. Selected model presentation: REML.
# Exports: fitted models, candidate tables, diagnostics, combined figures.
# Run in RStudio with Source. Designed for FSSgam >= 1.2.0.

# ---- 0. USER SETTINGS --------------------------------------------------------
PROJECT_DIR <- "input" # Where user wants project directory
DATA_FILE <- file.path(PROJECT_DIR, "input") # Add master csv data file
OUT_DIR <- file.path(PROJECT_DIR, "Combined_Reefs_GAMMS")

RUN_MODELS <- TRUE  # FALSE = rebuild combined plots from saved best models
RUN_DHARMA <- TRUE
RUN_SPATIAL_TEST <- TRUE
INCLUDE_SITE_RANDOM_EFFECT <- FALSE # inspect grouping first; TRUE adds s(site,bs="re")
SEED <- 123
MAX_PREDICTORS <- 2L  # ecological candidate predictors only; habitat random baseline
MIN_ROWS <- 25L
SMOOTH_K <- 3L
CORRELATION_CUTOFF <- 0.7  # maximum correlation allowed in same candidate model
FISHING_PREDICTOR <- "New_Log_RS"  

# Logged responses: Gaussian; OPTIONAL Gamma only if every value > 0.
# Prefer analysing ORIGINAL biomass using Tweedie when those columns exist.
TEST_GAMMA_ON_LOGGED <- FALSE
# Switch to TRUE to replace each *_log response with the corresponding raw
# response if present. Raw biomass: Tweedie (+Gamma if strictly positive);
# raw abundance: Poisson + negative binomial.
USE_RAW_RESPONSES <- FALSE
if (USE_RAW_RESPONSES) {
  OUT_DIR <- file.path(PROJECT_DIR, "Combined_Reefs_GAMMS_RAW")
}
# AICc comparisons are ONLY within exactly the same response & data rows.

# Spatial smooths in the baseline (keep consistent with existing study).
SPATIAL_TERMS <- "s(longitude, k=3) + s(latitude, k=3) + s(habitat, bs = 're') "
if (INCLUDE_SITE_RANDOM_EFFECT) {
  SPATIAL_TERMS <- paste0(SPATIAL_TERMS, " + s(site, bs='re')")
}

CONTINUOUS_PREDICTORS <- c(
  "coral_cover", "coral_shannon", "distance_to_estuary", "depth",
  "relief_max_500m", "slope_mean_500m",
  "prox_index_rocky_reef_500m_log"
)

RESPONSES <- c(
  "targeted_abundance_log", "targeted_biomass_log",
  "targeted_richness", "targeted_shannon",
  "nontargeted_abundance_log", "nontargeted_biomass_log",
  "nontargeted_richness", "nontargeted_shannon",
  "carangidae_biomass_log", "sparidae_biomass_log",
  "lutjanidae_biomass_log", "serranidae_biomass_log",
  "lethrinidae_biomass_log"
)

# ---- 1. PACKAGES AND DATA ----------------------------------------------------
NEEDED <- c("FSSgam", "mgcv", "MuMIn", "ggplot2", "cowplot")
missing <- NEEDED[!vapply(NEEDED, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}
if (RUN_DHARMA && !requireNamespace("DHARMa", quietly = TRUE)) {
  stop("Install DHARMa with install.packages('DHARMa') or set RUN_DHARMA <- FALSE")
}
set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(DATA_FILE) && RUN_MODELS) stop("Data file not found: ", DATA_FILE)
if (RUN_MODELS) {
  dat <- read.csv(DATA_FILE, check.names = TRUE, stringsAsFactors = FALSE)
  if (!"habitat" %in% names(dat)) stop("Column 'habitat' missing")
  dat$habitat <- factor(dat$habitat)
  if (INCLUDE_SITE_RANDOM_EFFECT) {
    if (!"site" %in% names(dat)) stop("site column needed for site random effect")
    dat$site <- factor(dat$site)
    if (nlevels(dat$site) < 3L) stop("Too few site levels for random effect")
  }
  if (nlevels(dat$habitat) < 2L) warning("Habitat has only one level")
  numeric_columns <- unique(c("latitude", "longitude", CONTINUOUS_PREDICTORS,
                              RESPONSES, sub("_log$", "", RESPONSES)))
  for (nm in intersect(numeric_columns, names(dat))) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  }
  message("Habitat breakdown: ")
  print(table(dat$habitat, useNA = "ifany"))
  missing_pred <- setdiff(CONTINUOUS_PREDICTORS, names(dat))
  if (length(missing_pred)) {
    warning("Missing predictors (excluded, not silently substituted): ",
            paste(missing_pred, collapse = ", "))
  }
}

# ---- 2. HELPERS --------------------------------------------------------------
safe_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
write_table <- function(x, path) write.csv(x, path, row.names = FALSE)

# Human-readable headings for variable-importance CSVs. Variable keys remain
# unchanged for joining back to the original data and model summaries.
PREDICTOR_HEADINGS <- c(
  coral_cover = "Coral cover",
  coral_shannon = "Coral Shannon diversity",
  distance_to_estuary = "Distance to estuary (m)",
  depth = "Depth (m)",
  relief_max_500m = "Maximum relief (500 m)",
  slope_mean_500m = "Mean slope (500 m)",
  prox_index_rocky_reef_500m_log = "Rocky reef proximity index",
  New_Log_RS = "Fishing Pressure Index"
)
pretty_predictor <- function(v) {
  ifelse(v %in% names(PREDICTOR_HEADINGS),
         unname(PREDICTOR_HEADINGS[v]), gsub("_", " ", v))
}

# Importances are the sums of normalized AICc model weights for ALL fitted,
# converged candidate models that contain a given ecological predictor.
# Spatial terms and the mandatory s(habitat, bs='re') are deliberately excluded.
# These describe importance within the candidate model set, not causal effects.
write_aicc_importance <- function(tab, fits, predictors, response, family, dest) {
  if (!length(predictors) || !nrow(tab)) return(invisible(NULL))
  model_weights <- as.numeric(tab$weight_AICc)
  present <- vapply(seq_len(nrow(tab)), function(i) {
    model <- fits[[as.character(tab$model_id[i])]]
    if (is.null(model)) stop("Missing model for importance calculation: ", tab$model_id[i])
    predictors %in% all.vars(stats::formula(model))
  }, logical(length(predictors)))
  # vapply makes a predictor-by-model matrix, even for a single model.
  present <- matrix(present, nrow = length(predictors), ncol = nrow(tab))
  rownames(present) <- predictors
  scores <- as.numeric(present %*% model_weights)
  selected_fit <- fits[[as.character(tab$model_id[1])]]
  selected_vars <- intersect(predictors, all.vars(stats::formula(selected_fit)))
  imp <- data.frame(
    response = response, family = family,
    predictor = predictors,
    predictor_heading = as.character(pretty_predictor(predictors)),
    AICc_importance = scores,
    candidate_models_including_predictor = rowSums(present),
    candidate_models_compared = nrow(tab),
    retained_in_best_model = predictors %in% selected_vars,
    stringsAsFactors = FALSE
  )
  imp <- imp[order(-imp$AICc_importance, imp$predictor), , drop = FALSE]
  imp$importance_rank <- seq_len(nrow(imp))
  write_table(imp, file.path(dest, "VARIABLE_IMPORTANCE_AICc.csv"))
  wide <- as.data.frame(as.list(setNames(imp$AICc_importance, imp$predictor)),
                        check.names = FALSE)
  wide <- cbind(data.frame(response = response, family = family), wide)
  write_table(wide, file.path(dest, "VARIABLE_IMPORTANCE_AICc_WIDE.csv"))
  wide_labelled <- wide
  names(wide_labelled)[-(1:2)] <- as.character(pretty_predictor(names(wide)[-(1:2)]))
  write_table(wide_labelled, file.path(dest, "VARIABLE_IMPORTANCE_AICc_WIDE_LABELLED.csv"))
  invisible(imp)
}

# FSSgam stores variable importance as a nested object, not as a plain table.
# Keep the .rds AND export its native AIC (not AICc) weights separately.
write_fss_native_importance <- function(vi, response, family, dest) {
  w <- tryCatch(vi$aic$variable.weights.raw, error = function(e) NULL)
  if (is.null(w)) return(invisible(FALSE))
  if (is.matrix(w) || is.data.frame(w)) {
    mm <- as.matrix(w)
    if (nrow(mm) == 1L && !is.null(colnames(mm))) {
      vals <- as.numeric(mm[1L, ])
      keys <- colnames(mm)
    } else if (ncol(mm) == 1L && !is.null(rownames(mm))) {
      vals <- as.numeric(mm[, 1L])
      keys <- rownames(mm)
    } else return(invisible(FALSE))
  } else if (is.numeric(w) && !is.null(names(w))) {
    vals <- as.numeric(w)
    keys <- names(w)
  } else return(invisible(FALSE))
  result <- data.frame(
    response = response, family = family, predictor = keys,
    predictor_heading = as.character(pretty_predictor(keys)),
    FSSgam_AIC_importance = vals, stringsAsFactors = FALSE
  )
  write_table(result, file.path(dest, "FSSgam_VARIABLE_IMPORTANCE_AIC.csv"))
  invisible(TRUE)
}
response_families <- function(resp, y) {
  z <- y[is.finite(y)]
  if (!length(z)) return(list())
  nonnegative <- all(z >= 0)
  positive <- all(z > 0)
  integerish <- nonnegative && all(abs(z - round(z)) < 1e-8)
  if (grepl("richness$", resp)) {
    if (!integerish) stop("Noninteger or negative richness: ", resp)
    return(list(poisson = poisson(link = "log"), nb = mgcv::nb()))
  }
  if (grepl("shannon$", resp)) return(list(gaussian = gaussian()))
  if (grepl("_log$", resp)) {
    ans <- list(gaussian = gaussian())
    # A Gamma model is INVALID if even ONE exact zero occurs.
    if (TEST_GAMMA_ON_LOGGED && positive) ans$gamma <- Gamma(link = "log")
    return(ans)
  }
  if (grepl("abundance$", resp)) {
    if (!integerish) stop("Noninteger/negative raw abundance: ", resp)
    return(list(poisson = poisson(link = "log"), nb = mgcv::nb()))
  }
  if (grepl("biomass$", resp)) {
    if (!nonnegative) stop("Negative raw biomass: ", resp)
    ans <- list(tweedie = mgcv::tw(link = "log"))
    if (positive) ans$gamma <- Gamma(link = "log")
    return(ans)
  }
  list(gaussian = gaussian())
}

# Do not select on different datasets. Generate one complete-case frame per
# response, then reuse exactly that frame across its candidate families.
get_analysis_frame <- function(data, resp) {
  preds <- intersect(CONTINUOUS_PREDICTORS, names(data))
  if (!all(c("longitude", "latitude", "habitat", resp) %in% names(data))) {
    stop("Required response / spatial / habitat column missing for ", resp)
  }
  # Exclude predictors with >25% missingness BEFORE complete-case selection.
  preds <- preds[vapply(preds, function(v) mean(!is.finite(data[[v]])) <= .25,
                        logical(1))]
  preds <- preds[vapply(preds, function(v) {
    a <- data[[v]][is.finite(data[[v]])]
    length(unique(a)) >= 4L
  }, logical(1))]
  cols <- unique(c(resp, "latitude", "longitude", "habitat",
                   if (INCLUDE_SITE_RANDOM_EFFECT) "site", preds))
  use.dat <- data[, cols, drop = FALSE]
  for (v in setdiff(cols, c(resp, "habitat", "site"))) {
    use.dat[[v]][!is.finite(use.dat[[v]])] <- NA_real_
  }
  use.dat[[resp]][!is.finite(use.dat[[resp]])] <- NA_real_
  use.dat <- droplevels(use.dat[stats::complete.cases(use.dat), , drop = FALSE])
  if (INCLUDE_SITE_RANDOM_EFFECT && nlevels(use.dat$site) < 3L) {
    stop("Too few site levels after complete cases")
  }
  if (nrow(use.dat) < MIN_ROWS) {
    stop(sprintf("%s: n=%d below minimum %d", resp, nrow(use.dat), MIN_ROWS))
  }
  if (length(unique(use.dat[[resp]])) < 3L) stop("Response lacks variation: ", resp)
  if (length(unique(use.dat$latitude)) < 4L ||
      length(unique(use.dat$longitude)) < 4L) {
    stop("Not enough distinct coordinates for spatial baseline: ", resp)
  }
  # Reprune after missingness removal, including factor with one level.
  preds <- preds[vapply(preds, function(v) length(unique(use.dat[[v]])) >= 4L,
                        logical(1))]
  use.dat <- use.dat[, unique(c(resp, "latitude", "longitude", "habitat",
                                if (INCLUDE_SITE_RANDOM_EFFECT) "site", preds)),
                     drop = FALSE]
  use.dat$response <- use.dat[[resp]]
  rownames(use.dat) <- NULL
  # Habitat remains in use.dat for s(habitat, bs='re') in the baseline,
  # but must NOT be passed to FSSgam as a candidate fixed factor.
  list(df = use.dat, predictors = preds, factors = character(0))
}

# FSSgam >=1.2.0 exports underscore functions; aliases remain for older code.
fss_generate <- function(...) {
  fn <- if ("generate_model_set" %in% getNamespaceExports("FSSgam")) {
    getExportedValue("FSSgam", "generate_model_set")
  } else getExportedValue("FSSgam", "generate.model.set")
  fn(...)
}
fss_fit <- function(...) {
  fn <- if ("fit_model_set" %in% getNamespaceExports("FSSgam")) {
    getExportedValue("FSSgam", "fit_model_set")
  } else getExportedValue("FSSgam", "fit.model.set")
  fn(...)
}

# ---- 3. FIT A SINGLE RESPONSE / FAMILY --------------------------------------
fit_one_family <- function(frame, response, family_name, family_object) {
  dd <- frame$df
  pred <- frame$predictors
  fac <- frame$factors
  dest <- file.path(OUT_DIR, response, family_name)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest, "diagnostics"), showWarnings = FALSE)
  dir.create(file.path(dest, "models"), showWarnings = FALSE)
  
  null_formula <- as.formula(paste("response ~", SPATIAL_TERMS))
  null_fit <- mgcv::gam(null_formula, family = family_object, data = dd,
                        method = "ML", na.action = na.fail)
  if (!isTRUE(null_fit$converged)) stop("Null failed to converge")
  # FSSgam updates the original gam call. Put a literal family object in it.
  null_fit$call$family <- family_object
  
  result <- NULL
  if (length(pred) + length(fac)) {
    model_set <- fss_generate(
      use.dat = dd, test.fit = null_fit,
      pred.vars.cont = if (length(pred)) pred else NA,
      pred.vars.fact = NULL,
      null.terms = SPATIAL_TERMS,
      max.predictors = min(MAX_PREDICTORS, length(pred) + length(fac)),
      k = SMOOTH_K, cov.cutoff = CORRELATION_CUTOFF,
      factor.smooth.interactions = FALSE,
      factor.factor.interactions = FALSE,
      smooth.smooth.interactions = FALSE
    )
    n_candidates <- length(model_set$formulas)
    message("  ", response, " [", family_name, "]: ", n_candidates,
            " FSS candidate formulas")
    result <- fss_fit(model_set, parallel = FALSE,
                      max.models = max(200L, n_candidates + 10L),
                      save.model.fits = TRUE)
    if (!is.null(result$mod.data.out)) {
      write_table(result$mod.data.out,
                  file.path(dest, "FSSgam_original_model_table.csv"))
    }
    if (!is.null(result$variable.importance)) {
      saveRDS(result$variable.importance,
              file.path(dest, "FSSgam_variable_importance.rds"))
    }
    if (!is.null(result$variable.importance)) {
      tryCatch(
        write_fss_native_importance(result$variable.importance, response, family_name, dest),
        error = function(e) warning("Could not export native FSSgam importance: ", conditionMessage(e))
      )
    }
    if (length(result$failed.models)) {
      writeLines(capture.output(str(result$failed.models)),
                 file.path(dest, "FSSgam_failed_models.txt"))
    }
  }
  
  # Explicitly include spatial null even if FSSgam did not save it.
  fits <- list(SPATIAL_NULL = null_fit)
  if (!is.null(result) && length(result$success.models)) {
    fits <- c(fits, result$success.models)
  }
  fits <- fits[vapply(fits, inherits, logical(1), what = "gam")]
  if (!length(fits)) stop("No usable GAM models")
  
  # Recompute ALL candidate AICc from their actual saved ML model objects.
  # Ensures the best-fit RDS truly matches the ranked model.
  tab <- do.call(rbind, lapply(seq_along(fits), function(i) {
    m <- fits[[i]]
    aicc <- tryCatch(as.numeric(MuMIn::AICc(m)), error = function(e) NA_real_)
    data.frame(model_id = names(fits)[i], family = family_name,
               AICc = aicc, AIC = stats::AIC(m), n = stats::nobs(m),
               n_parameters = attr(stats::logLik(m), "df"),
               converged = isTRUE(m$converged),
               formula = paste(deparse(stats::formula(m)), collapse = " "),
               deviance_explained_ML = summary(m)$dev.expl,
               stringsAsFactors = FALSE)
  }))
  tab <- tab[is.finite(tab$AICc) & tab$converged, , drop = FALSE]
  if (!nrow(tab)) stop("No converged, finite-AICc models")
  if (length(unique(tab$n)) != 1L || unique(tab$n) != nrow(dd)) {
    stop("Model n differs across candidates (AICc comparison invalid)")
  }
  # Remove null duplicates, if they occurred, by exact formula string.
  tab <- tab[!duplicated(tab$formula), , drop = FALSE]
  tab <- tab[order(tab$AICc, tab$n_parameters), , drop = FALSE]
  tab$delta_AICc <- tab$AICc - min(tab$AICc)
  tab$weight_AICc <- exp(-0.5 * tab$delta_AICc)
  tab$weight_AICc <- tab$weight_AICc / sum(tab$weight_AICc)
  tab$within_delta2 <- tab$delta_AICc <= 2
  tab$selection <- ifelse(seq_len(nrow(tab)) == 1L, "minimum_AICc", "")
  # Alternative simplest candidate among those with Delta AICc <= 2.
  compatible <- which(tab$within_delta2)
  simpler <- compatible[which.min(tab$n_parameters[compatible])]
  tab$parsimonious_within_2 <- seq_len(nrow(tab)) == simpler
  write_table(tab, file.path(dest, "ALL_CANDIDATES_RANKED_AICc.csv"))
  write_table(subset(tab, within_delta2),
              file.path(dest, "TOP_MODELS_DELTA_AICc_LE_2.csv"))
  
  # Readable variable-importance exports for every response and family.
  # Uses the SAME deduplicated, converged ML candidate set as model selection.
  write_aicc_importance(tab, fits, pred, response, family_name, dest)
  
  winner <- tab[1, , drop = FALSE]
  best_ml <- fits[[winner$model_id]]
  if (is.null(best_ml) || !inherits(best_ml, "gam")) {
    stop("Selected model was not found: ", winner$model_id)
  }
  saveRDS(best_ml, file.path(dest, "models", "BEST_ML.rds"))
  best_reml <- tryCatch(update(best_ml, method = "REML"),
                        error = function(e) {
                          warning("REML refit failed for ", response, " / ", family_name,
                                  ": ", conditionMessage(e))
                          NULL
                        })
  if (!is.null(best_reml) && isTRUE(best_reml$converged)) {
    saveRDS(best_reml, file.path(dest, "models", "BEST_REML.rds"))
  } else {
    # Keep ML model as a usable fallback, with accurate metadata.
    best_reml <- best_ml
    saveRDS(best_reml, file.path(dest, "models", "BEST_REML.rds"))
  }
  
  model_summary <- capture.output(summary(best_reml))
  writeLines(c("Selected ML candidate (minimum AICc):", winner$formula,
               "", model_summary), file.path(dest, "BEST_MODEL_SUMMARY.txt"))
  
  # Concurvity; k-diagnostics and DHARMa on selected presentation model.
  cv <- tryCatch(mgcv::concurvity(best_reml, full = TRUE),
                 error = function(e) NULL)
  if (!is.null(cv)) {
    write.csv(cv, file.path(dest, "diagnostics", "CONCURVITY.csv"))
  }
  ktab <- tryCatch(mgcv::k.check(best_reml), error = function(e) NULL)
  if (!is.null(ktab)) {
    write.csv(ktab, file.path(dest, "diagnostics", "K_CHECK.csv"))
  }
  grDevices::png(file.path(dest, "diagnostics", "GAM_CHECK.png"),
                 width = 1350, height = 1000, res = 150)
  try(mgcv::gam.check(best_reml), silent = TRUE)
  grDevices::dev.off()
  
  if (RUN_DHARMA) {
    set.seed(SEED)
    sim <- tryCatch(DHARMa::simulateResiduals(best_reml, n = 1000),
                    error = function(e) {
                      warning("DHARMa failure: ", conditionMessage(e))
                      NULL
                    })
    if (!is.null(sim)) {
      grDevices::png(file.path(dest, "diagnostics", "DHARMA_OVERVIEW.png"),
                     width = 1300, height = 950, res = 150)
      plot(sim)
      grDevices::dev.off()
      fmt_test <- function(fun) tryCatch(
        capture.output(print(fun(sim))),
        error = function(e) paste("Error:", conditionMessage(e)))
      diagnostics <- c("Uniformity:", fmt_test(DHARMa::testUniformity),
                       "", "Dispersion:", fmt_test(DHARMa::testDispersion),
                       "", "Outliers:", fmt_test(DHARMa::testOutliers),
                       "", "Zero inflation:", fmt_test(DHARMa::testZeroInflation))
      if (RUN_SPATIAL_TEST && requireNamespace("sf", quietly = TRUE)) {
        coords <- best_reml$model[, c("longitude", "latitude"), drop = FALSE]
        if (all(stats::complete.cases(coords)) &&
            !anyDuplicated(paste(coords$longitude, coords$latitude))) {
          locs <- sf::st_as_sf(coords, coords = c("longitude", "latitude"), crs = 4326)
          locs <- sf::st_transform(locs, 32756)
          xy <- sf::st_coordinates(locs)
          diagnostics <- c(diagnostics, "", "Spatial autocorrelation (exploratory):",
                           tryCatch(capture.output(print(DHARMa::testSpatialAutocorrelation(
                             sim, x = xy[, 1], y = xy[, 2]))),
                             error = function(e) conditionMessage(e)))
        } else {
          diagnostics <- c(diagnostics, "Spatial test skipped: duplicate/missing coordinates")
        }
      }
      writeLines(diagnostics, file.path(dest, "diagnostics", "DHARMA_TESTS.txt"))
    }
  }
  
  data.frame(response = response, family = family_name,
             model_id = winner$model_id, formula = winner$formula,
             n = winner$n, AICc = winner$AICc,
             delta_within_family = 0,
             within_family_AICc_weight = winner$weight_AICc,
             deviance_explained_REML = summary(best_reml)$dev.expl,
             null_is_best = identical(as.character(winner$model_id), "SPATIAL_NULL"),
             fishing_predictor_selected = grepl(FISHING_PREDICTOR, winner$formula, fixed = TRUE),
             n_candidate_models = nrow(tab),
             parsimonious_within_2 = tab$model_id[simpler],
             file = file.path(dest, "models", "BEST_REML.rds"),
             stringsAsFactors = FALSE)
}

# ---- 4. FIT ALL RESPONSES; SELECT FAMILY BY COMPARABLE ML AICc ----------------
if (RUN_MODELS) {
  if (USE_RAW_RESPONSES) {
    RESPONSES <- unique(ifelse(grepl("_log$", RESPONSES),
                               sub("_log$", "", RESPONSES), RESPONSES))
  }
  overall <- list()
  failed <- character()
  for (resp in RESPONSES) {
    message("\n============= RESPONSE: ", resp, " =============")
    if (!resp %in% names(dat)) {
      failed <- c(failed, paste(resp, "missing response column"))
      next
    }
    frame <- tryCatch(get_analysis_frame(dat, resp), error = function(e) e)
    if (inherits(frame, "error")) {
      failed <- c(failed, paste(resp, conditionMessage(frame)))
      next
    }
    write_table(data.frame(response = resp, n = nrow(frame$df),
                           missing_response = sum(is.na(dat[[resp]])),
                           selected_predictors = paste(frame$predictors, collapse = "; "),
                           stringsAsFactors = FALSE),
                file.path(OUT_DIR, paste0("INPUT_AUDIT_", safe_name(resp), ".csv")))
    families <- tryCatch(response_families(resp, frame$df$response),
                         error = function(e) e)
    if (inherits(families, "error")) {
      failed <- c(failed, paste(resp, conditionMessage(families)))
      next
    }
    results <- list()
    for (fn in names(families)) {
      message("Testing ", resp, " as ", fn)
      answer <- tryCatch(fit_one_family(frame, resp, fn, families[[fn]]),
                         error = function(e) e)
      if (inherits(answer, "error")) {
        failed <- c(failed, paste(resp, fn, conditionMessage(answer)))
        warning("FAILED: ", resp, " / ", fn, ": ", conditionMessage(answer))
      } else {
        results[[fn]] <- answer
      }
    }
    if (!length(results)) next
    ff <- do.call(rbind, results)
    # Valid because SAME raw response values, transformation and observation
    # rows across the candidate families. Never compare logged-vs-raw AICc.
    ff <- ff[order(ff$AICc), , drop = FALSE]
    ff$delta_family_AICc <- ff$AICc - min(ff$AICc)
    ff$family_AICc_weight <- exp(-.5 * ff$delta_family_AICc)
    ff$family_AICc_weight <- ff$family_AICc_weight / sum(ff$family_AICc_weight)
    ff$best_family <- seq_len(nrow(ff)) == 1L
    write_table(ff, file.path(OUT_DIR, resp, "FAMILY_COMPARISON.csv"))
    overall[[resp]] <- ff
  }
  if (length(overall)) {
    all_results <- do.call(rbind, overall)
    write_table(all_results, file.path(OUT_DIR, "ALL_FAMILY_COMPARISONS.csv"))
    best <- all_results[all_results$best_family, , drop = FALSE]
    rownames(best) <- NULL
    write_table(best, file.path(OUT_DIR, "BEST_FAMILY_AND_MODEL_BY_RESPONSE.csv"))
    
    # MASTER importance outputs: one row per predictor/response and one
    # response-by-predictor matrix. Only the selected distribution for each
    # response is used; alternate family comparisons remain in their folders.
    all_imp <- list()
    for (i in seq_len(nrow(best))) {
      imp_file <- file.path(OUT_DIR, best$response[i], best$family[i],
                            "VARIABLE_IMPORTANCE_AICc.csv")
      if (file.exists(imp_file)) {
        all_imp[[length(all_imp) + 1L]] <- read.csv(imp_file, check.names = FALSE)
      }
    }
    if (length(all_imp)) {
      master <- do.call(rbind, all_imp)
      write_table(master, file.path(OUT_DIR, "VARIABLE_IMPORTANCE_BEST_FAMILY_LONG.csv"))
      
      # NA means predictor was not in that response's candidate set,
      # whereas 0 means it was eligible but had no AICc weight support.
      keys <- unique(CONTINUOUS_PREDICTORS)
      wide_rows <- lapply(seq_len(nrow(best)), function(i) {
        sub <- master[master$response == best$response[i] &
                        master$family == best$family[i], , drop = FALSE]
        values <- stats::setNames(rep(NA_real_, length(keys)), keys)
        values[sub$predictor] <- sub$AICc_importance
        data.frame(response = best$response[i], family = best$family[i],
                   as.list(values), check.names = FALSE)
      })
      master_wide <- do.call(rbind, wide_rows)
      write_table(master_wide,
                  file.path(OUT_DIR, "VARIABLE_IMPORTANCE_BEST_FAMILY_WIDE.csv"))
      master_labelled <- master_wide
      names(master_labelled)[-(1:2)] <- as.character(pretty_predictor(names(master_wide)[-(1:2)]))
      write_table(master_labelled,
                  file.path(OUT_DIR, "VARIABLE_IMPORTANCE_BEST_FAMILY_WIDE_LABELLED.csv"))
      write_table(data.frame(
        predictor = keys, predictor_heading = as.character(pretty_predictor(keys)),
        stringsAsFactors = FALSE),
        file.path(OUT_DIR, "VARIABLE_IMPORTANCE_PREDICTOR_HEADINGS.csv"))
    }
  }
  writeLines(if (length(failed)) failed else "No failed/skipped combinations",
             file.path(OUT_DIR, "MODEL_FAILURES_AND_SKIPS.txt"))
}


# STRUCTURED FIGURES: SHOW ALL SELECTED ECOLOGICAL EFFECTS

# Purpose:
#   (1) Community metrics: Targets vs non-targets, 4 metric rows
#       Each response gets TWO side-by-side predictor slots.
#   (2) Fish families: five family biomass rows, TWO predictor slots per row.
#
# Uses EXISTING saved BEST_REML.rds model objects. No model refitting.
# Shows all selected non-spatial effects, including habitat if selected.
# Excludes the baseline spatial smooths and optional site random effect.
# Does NOT mistake habitat being present in model$model for it being selected.
# Stops rather than silently discarding an effect if a model has >2 effects.
#
# Files written to: <MODEL_RESULTS>/BEST_FIT_STRUCTURED_ALL_EFFECTS/


# ---- 1. USER SETTINGS --------------------------------------------------------
PROJECT_DIR <- "input" # Where users project directory is
MODEL_RESULTS <- "Combined_Reefs_GAMMS" # or "combined_reef_gamms" / "combined_reef_gamms_REVISED_RAW"
ROOT <- file.path(PROJECT_DIR, MODEL_RESULTS)
FIG_DIR <- file.path(ROOT, "BEST_FIT_Summary_Figures")
SLOTS <- 2L               # original FSS script allows <= 2 ecological predictors
COLOUR <- "#ad1735"      # similar to supplied manuscript example
CI_LEVEL <- 0.95
SHOW_MISSING_SLOT_TEXT <- FALSE  # leave genuinely unselected second predictor blank

for (pkg in c("ggplot2", "cowplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install: ", pkg)
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 2. READ MODEL SELECTION -------------------------------------------------
selected_file <- file.path(ROOT, "BEST_FAMILY_AND_MODEL_BY_RESPONSE.csv")
if (file.exists(selected_file)) {
  selected <- read.csv(selected_file, stringsAsFactors = FALSE)
  if (!all(c("response", "family", "file") %in% names(selected))) {
    stop("Selection file needs columns: response, family, file")
  }
} else {
  # Compatibility with the original (pre-revision) FSS script.
  selected_file <- file.path(ROOT, "FAMILY_SELECTION.csv")
  if (!file.exists(selected_file)) stop("No model selection CSV in: ", ROOT)
  selected <- read.csv(selected_file, stringsAsFactors = FALSE)
  selected$file <- file.path(ROOT, selected$response, selected$family, "models",
                             paste0("best_REML_", selected$family, "_", selected$response, ".rds"))
  warning("Older FAMILY_SELECTION.csv used; previously selected models can reflect old ranking.")
}

if (anyDuplicated(selected$response)) stop("Multiple selected model rows for the same response")

resolve_model_file <- function(f) {
  if (is.na(f) || !nzchar(f)) return(NA_character_)
  if (file.exists(f)) return(f)
  rel <- file.path(ROOT, f)
  if (file.exists(rel)) return(rel)
  f
}
selected$file <- vapply(selected$file, resolve_model_file, character(1))

# ---- 3. PRESENTATION HELPERS -------------------------------------------------
labelling <- c(
  "depth" = "Depth (m)",
  "relief_max_500m" = "Maximum relief (500 m)",
  "slope_mean_500m" = "Mean slope (500 m)",
  "distance_to_estuary" = "Distance to estuary (m)",
  "coral_cover" = "Coral cover",
  "coral_shannon" = "Coral Shannon diversity",
  "prox_index_rocky_reef_500m_log" = "Rocky reef proximity index",
  "habitat" = "Habitat"
)
pretty_term <- function(x) {
  if (x %in% names(labelling)) return(unname(labelling[[x]]))
  gsub("_", " ", x)
}

blank <- function(note = NULL) {
  if (is.null(note) || !nzchar(note)) return(ggplot2::ggplot() + ggplot2::theme_void())
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
    ggplot2::geom_text(label = note, colour = "grey55", size = 3.4) +
    ggplot2::theme_void()
}

strip_text <- function(label, angle = 0, size = 5.2) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
    ggplot2::geom_text(label = label, fontface = "bold", angle = angle, size = size) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(clip = "off")
}

style_panel <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 9),
      axis.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(5, 8, 5, 5)
    )
}

# Selected terms are taken from predict(type="terms") for the FITTED MODEL,
# NOT inferred from columns stored in the input model frame.
get_ecological_terms <- function(model) {
  d <- model$model
  if (is.null(d) || !nrow(d)) stop("Fitted model has no model frame")
  term_fit <- stats::predict(model, newdata = d[1, , drop = FALSE], type = "terms")
  term_names <- colnames(as.matrix(term_fit))
  if (!length(term_names)) return(character(0))
  keep <- vapply(term_names, function(term) {
    # Check smooth metadata first; this correctly identifies spatial and RE terms.
    matched <- which(vapply(model$smooth, function(s) identical(s$label, term), logical(1)))
    if (length(matched)) {
      s <- model$smooth[[matched[1]]]
      if (inherits(s, "random.effect")) return(FALSE)
      if (any(s$term %in% c("latitude", "longitude", "site"))) return(FALSE)
      return(TRUE)
    }
    # Parametric baseline spatial terms are also excluded, if present.
    !(term %in% c("latitude", "longitude", "site"))
  }, logical(1))
  unique(term_names[keep])
}

make_reference_frame <- function(model, n = 150L) {
  d <- model$model
  nd <- d[rep(1L, n), , drop = FALSE]
  for (nm in names(d)) {
    if (is.factor(d[[nm]])) {
      tab <- sort(table(d[[nm]]), decreasing = TRUE)
      level <- names(tab)[1]
      nd[[nm]] <- factor(rep(level, n), levels = levels(d[[nm]]))
    } else if (is.numeric(d[[nm]])) {
      nd[[nm]] <- rep(stats::median(d[[nm]], na.rm = TRUE), n)
    } else if (is.character(d[[nm]])) {
      nd[[nm]] <- rep(d[[nm]][1], n)
    }
  }
  nd
}

# All panels plot the PARTIAL term effect on the LINEAR PREDICTOR scale.
# Bands/intervals are approximate 95% pointwise confidence intervals.
plot_effect <- function(model, term) {
  d <- model$model
  sm_idx <- which(vapply(model$smooth, function(s) identical(s$label, term), logical(1)))
  if (length(sm_idx)) {
    sm <- model$smooth[[sm_idx[1]]]
    if (length(sm$term) != 1L) return(blank(paste("Complex effect:", term)))
    var <- sm$term[[1]]
  } else {
    var <- term
  }
  if (!(var %in% names(d))) return(blank(paste("Unavailable term:", term)))
  
  ngrid <- 150L
  nd <- make_reference_frame(model, n = ngrid)
  var_data <- d[[var]]
  if (is.factor(var_data)) {
    levels_seen <- levels(droplevels(var_data))
    nd <- nd[seq_along(levels_seen), , drop = FALSE]
    nd[[var]] <- factor(levels_seen, levels = levels(var_data))
    pr <- stats::predict(model, newdata = nd, type = "terms", terms = term, se.fit = TRUE)
    if (!(term %in% colnames(pr$fit))) stop("Term not present in predictions: ", term)
    z <- data.frame(x = nd[[var]], fit = as.numeric(pr$fit[, term]),
                    se = as.numeric(pr$se.fit[, term]))
    ggplot2::ggplot(z, ggplot2::aes(x = x, y = fit)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      ggplot2::geom_pointrange(ggplot2::aes(
        ymin = fit - stats::qnorm(1 - (1 - CI_LEVEL) / 2) * se,
        ymax = fit + stats::qnorm(1 - (1 - CI_LEVEL) / 2) * se
      ), colour = COLOUR, linewidth = 0.5) +
      ggplot2::labs(x = pretty_term(var), y = "Estimated effect") + style_panel()
  } else if (is.numeric(var_data)) {
    x_range <- range(var_data[is.finite(var_data)])
    if (!all(is.finite(x_range)) || diff(x_range) == 0) {
      return(blank(paste("No variation:", pretty_term(var))))
    }
    nd[[var]] <- seq(x_range[1], x_range[2], length.out = ngrid)
    # Equivalent to mgcv's seWithMean = TRUE for the continuous smooth.
    # Do not change habitat/factor contrast behaviour.
    pr <- stats::predict(model, newdata = nd, type = "iterms",
                         terms = term, se.fit = TRUE, iterms.type = 1)
    if (!(term %in% colnames(pr$fit))) stop("Term not present in predictions: ", term)
    z <- data.frame(x = nd[[var]], fit = as.numeric(pr$fit[, term]),
                    se = as.numeric(pr$se.fit[, term]))
    half_width <- stats::qnorm(1 - (1 - CI_LEVEL) / 2) * z$se
    z$lower <- z$fit - half_width
    z$upper <- z$fit + half_width
    ggplot2::ggplot(z, ggplot2::aes(x = x, y = fit)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           fill = COLOUR, alpha = .24) +
      ggplot2::geom_line(colour = COLOUR, linewidth = .9) +
      ggplot2::geom_rug(data = data.frame(x = var_data),
                        ggplot2::aes(x = x), inherit.aes = FALSE,
                        sides = "b", alpha = .17, colour = COLOUR) +
      ggplot2::labs(x = pretty_term(var), y = "Estimated effect") + style_panel()
  } else {
    blank(paste("Unsupported term:", term))
  }
}

# ---- 4. LOAD SELECTED MODELS ONCE --------------------------------------------
responses_community <- list(
  Abundance = c("targeted_abundance_log", "nontargeted_abundance_log"),
  Biomass = c("targeted_biomass_log", "nontargeted_biomass_log"),
  Richness = c("targeted_richness", "nontargeted_richness"),
  Shannon = c("targeted_shannon", "nontargeted_shannon")
)
responses_family <- c(
  "carangidae_biomass_log", "sparidae_biomass_log", "lutjanidae_biomass_log",
  "serranidae_biomass_log", "lethrinidae_biomass_log"
)
needed_responses <- unique(c(unlist(responses_community), responses_family))
models <- list()
terms_by_response <- list()
audit <- list()
issues <- character(0)

for (resp in needed_responses) {
  sel <- selected[selected$response == resp, , drop = FALSE]
  if (!nrow(sel)) {
    issues <- c(issues, paste(resp, "not in selected-model summary"))
    next
  }
  file <- sel$file[1]
  if (!file.exists(file)) {
    issues <- c(issues, paste(resp, "RDS missing:", file))
    next
  }
  m <- tryCatch(readRDS(file), error = function(e) e)
  if (inherits(m, "error") || !inherits(m, "gam")) {
    issues <- c(issues, paste(resp, "could not read fitted gam model"))
    next
  }
  trms <- tryCatch(get_ecological_terms(m), error = function(e) e)
  if (inherits(trms, "error")) {
    issues <- c(issues, paste(resp, "term identification failed:", conditionMessage(trms)))
    next
  }
  if (length(trms) > SLOTS) {
    stop(resp, " has ", length(trms), " selected ecological effects: ",
         paste(trms, collapse = ", "), ". Increase SLOTS/layout; no term omitted.")
  }
  models[[resp]] <- m
  terms_by_response[[resp]] <- trms
  audit[[resp]] <- data.frame(response = resp, family = sel$family[1],
                              n_selected_terms = length(trms),
                              effect_1 = if (length(trms) >= 1) trms[1] else NA_character_,
                              effect_2 = if (length(trms) >= 2) trms[2] else NA_character_,
                              model_file = file, stringsAsFactors = FALSE)
}

if (length(audit)) {
  utils::write.csv(do.call(rbind, audit),
                   file.path(FIG_DIR, "SELECTED_ECOLOGICAL_EFFECTS_AUDIT.csv"),
                   row.names = FALSE)
}
writeLines(if (length(issues)) issues else "All requested model files located.",
           file.path(FIG_DIR, "MISSING_MODELS_AND_PLOT_WARNINGS.txt"))

# ---- 5. PANEL MAKER (ALL SELECTED EFFECTS, NOT JUST FIRST) ------------------
render_panel <- function(resp, slot) {
  if (!(resp %in% names(models))) return(blank("Selected model unavailable"))
  trms <- terms_by_response[[resp]]
  if (length(trms) == 0L) {
    return(if (slot == 1L) blank("Spatial-only best model") else blank())
  }
  if (slot > length(trms)) {
    return(if (SHOW_MISSING_SLOT_TEXT) blank("No second selected effect") else blank())
  }
  tryCatch(plot_effect(models[[resp]], trms[[slot]]), error = function(e) {
    warning(resp, " / ", trms[[slot]], ": ", conditionMessage(e))
    blank(paste("Plot error:", trms[[slot]]))
  })
}

# ---- 6. FIGURE A: COMMUNITY METRICS -----------------------------------------
# Each response has 2 side-by-side panels: one for EACH selected predictor.
# Layout: [metric] [targets term 1][targets term 2] [non-targets term 1][non-targets term 2]
community_header <- cowplot::plot_grid(
  blank(),
  cowplot::plot_grid(strip_text("Fisheries targets", size = 5.5),
                     strip_text("Non-targets", size = 5.5), ncol = 2),
  ncol = 2, rel_widths = c(.23, 4)
)

community_rows <- lapply(names(responses_community), function(metric) {
  resp <- responses_community[[metric]]
  cowplot::plot_grid(
    strip_text(metric, angle = 90, size = 5.0),
    render_panel(resp[1], 1), render_panel(resp[1], 2),
    render_panel(resp[2], 1), render_panel(resp[2], 2),
    ncol = 5, rel_widths = c(.23, 1, 1, 1, 1)
  )
})
community_figure <- cowplot::plot_grid(
  plotlist = c(list(community_header), community_rows),
  ncol = 1, rel_heights = c(.15, 1, 1, 1, 1)
)

# ---- 7. FIGURE B: FISH FAMILIES ---------------------------------------------
# One row per family, with the TWO selected effects side by side.
family_name <- function(resp) {
  z <- sub("_biomass(_log)?$", "", resp)
  tools::toTitleCase(z)
}
family_header <- cowplot::plot_grid(
  strip_text("Fish family", size = 5),
  strip_text("Selected effect 1", size = 4.3),
  strip_text("Second selected term", size = 4.3),
  ncol = 3, rel_widths = c(.26, 1, 1)
)
family_rows <- lapply(responses_family, function(resp) {
  cowplot::plot_grid(strip_text(family_name(resp), size = 4.4),
                     render_panel(resp, 1), render_panel(resp, 2),
                     ncol = 3, rel_widths = c(.26, 1, 1))
})
family_figure <- cowplot::plot_grid(
  plotlist = c(list(family_header), family_rows),
  ncol = 1, rel_heights = c(.17, 1, 1, 1, 1, 1)
)

# ---- 8. EXPORT PDF & HIGH-RES PNG -------------------------------------------
output <- list(
  COMMUNITY_METRICS_ALL_SELECTED_EFFECTS = list(fig = community_figure, w = 15.5, h = 13.5),
  FISH_FAMILIES_ALL_SELECTED_EFFECTS = list(fig = family_figure, w = 11.5, h = 16)
)
for (nm in names(output)) {
  x <- output[[nm]]
  ggplot2::ggsave(file.path(FIG_DIR, paste0(nm, ".pdf")), x$fig,
                  width = x$w, height = x$h, units = "in", limitsize = FALSE)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(nm, ".png")), x$fig,
                  width = x$w, height = x$h, units = "in", dpi = 300,
                  limitsize = FALSE)
}
message("\nDONE: BOTH FIGURES (all selected non-spatial effects) saved to: ", FIG_DIR)
message("Audit CSV confirms exactly which terms were plotted for each response.")
message("Separate y-axis scales are used across responses; all effects are on the link scale.")

