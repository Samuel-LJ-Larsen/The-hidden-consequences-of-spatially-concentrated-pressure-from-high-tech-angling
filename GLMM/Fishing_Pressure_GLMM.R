
# Table S1: Fishing pressure by reef habitat
# Ordered beta GLMM + pairwise habitat contrast + correctly scaled figure

# IMPORTANT DESIGN:
#   - Inference: glmmTMB ordered-beta (logit) GLMM, site random intercept.
#   - Figure: observed habitat means of SITE MEANS +/- SE among sites.
#   - Letters: the GLMM Rocky-versus-Indurated pairwise comparison.
# These are deliberately different quantities; caption MUST state this.
# Never plot emmeans' link-scale 'emmean' as a 0-1 response mean.
# In glmmTMB 1.1.14, even emmeans type='response' for ordbeta can differ from the exact ordered-beta expected response (glmmTMB issue #1332).



# Settings ----------------------------------------------------------------


DATA_FILE <- paste0(
"input" 
) # Input user's master csv
RESPONSE_COL <- "New_Log_RS"
SITE_COL <- "site"
HABITAT_COL <- "habitat"
ALPHA <- 0.05
N_SIM <- 500
SEED <- 123
OUT_DIR <- file.path(dirname(DATA_FILE), "Fishing_Pressure_Habitat_GLMM")

required <- c("glmmTMB", "emmeans", "ggplot2", "DHARMa")
missing_pkgs <- setdiff(required, rownames(installed.packages()))
if (length(missing_pkgs)) {
  stop("Install these packages first: ", paste(missing_pkgs, collapse = ", "))
}
if (!"ordbeta" %in% getNamespaceExports("glmmTMB")) {
  stop("Your glmmTMB version does not export ordbeta(). Update glmmTMB.")
}
if (!file.exists(DATA_FILE)) stop("Missing CSV: ", DATA_FILE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# Read and check data -----------------------------------------------------

raw <- read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c(RESPONSE_COL, SITE_COL, HABITAT_COL)
if (length(setdiff(needed, names(raw)))) {
  stop("Missing column(s): ", paste(setdiff(needed, names(raw)), collapse = ", "))
}

dat <- data.frame(
  fpi = suppressWarnings(as.numeric(as.character(raw[[RESPONSE_COL]]))),
  site = trimws(as.character(raw[[SITE_COL]])),
  habitat_original = trimws(as.character(raw[[HABITAT_COL]])),
  stringsAsFactors = FALSE
)
x <- tolower(dat$habitat_original)
dat$habitat <- ifelse(grepl("indurat", x), "Indurated reef",
               ifelse(grepl("rock", x), "Rocky reef", NA_character_))
unknown <- unique(dat$habitat_original[is.na(dat$habitat) &
                                       !is.na(dat$habitat_original) &
                                       nzchar(dat$habitat_original)])
if (length(unknown)) stop("Unexpected habitat labels: ", paste(unknown, collapse = ", "))

valid <- is.finite(dat$fpi) & !is.na(dat$habitat) &
         !is.na(dat$site) & nzchar(dat$site)
if (any(!valid)) message("Excluding ", sum(!valid), " rows missing model fields.")
dat <- dat[valid, , drop = FALSE]
dat$habitat <- factor(dat$habitat, levels = c("Rocky reef", "Indurated reef"))
dat$site <- factor(dat$site)
if (nlevels(droplevels(dat$habitat)) != 2) stop("Need both habitats.")
if (any(dat$fpi < 0 | dat$fpi > 1)) stop("Response outside [0,1].")

site_by_habitat <- table(dat$site, dat$habitat)
if (any(rowSums(site_by_habitat > 0) > 1)) {
  stop("At least one site occurs in both habitats; review site/habitat structure.")
}
nsites <- colSums(site_by_habitat > 0)
if (any(nsites < 2)) warning("Fewer than 2 sites in a habitat: between-site SE is unavailable.")
cat("\nResponse:", RESPONSE_COL,
    "\nRange:", paste(range(dat$fpi), collapse = " to "),
    "\nZeros:", sum(dat$fpi == 0),
    "\nOnes:", sum(dat$fpi == 1),
    "\nIndependent sites by habitat:\n")
print(nsites)
cat("\nReplicates by site/habitat:\n")
print(site_by_habitat)

# Ordered-beta glmm -------------------------------------------------------

set.seed(SEED)
fit <- glmmTMB::glmmTMB(
  fpi ~ habitat + (1 | site),
  data = dat,
  family = glmmTMB::ordbeta(link = "logit")
)
if (!isTRUE(fit$sdr$pdHess)) warning("Non-positive definite Hessian; check model fit.")
if (fit$fit$convergence != 0) warning("Model convergence code: ", fit$fit$convergence)
print(summary(fit))
saveRDS(fit, file.path(OUT_DIR, "Ordered_Beta_GLMM.rds"))



# Pairwise ----------------------------------------------------------------

# Test Rocky vs Indurated on the MODEL LINK SCALE, NOT response-scale means.
# There is exactly one pairwise comparison: Tukey and unadjusted p are equal.
emm <- emmeans::emmeans(fit, ~ habitat)
link_emm <- as.data.frame(summary(emm, infer = c(TRUE, FALSE)))
write.csv(link_emm, file.path(OUT_DIR, "Estimated_Marginal_Means_LINK_SCALE.csv"),
          row.names = FALSE)

pairwise <- pairs(emm, adjust = "tukey")
pairwise_df <- as.data.frame(summary(pairwise, infer = c(TRUE, TRUE),
                                      adjust = "tukey"))
if (nrow(pairwise_df) != 1L) stop("Expected one habitat contrast.")
stat_name <- intersect(c("z.ratio", "t.ratio"), names(pairwise_df))
if (length(stat_name) != 1L) stop("Expected z.ratio or t.ratio.")
p_pair <- pairwise_df$p.value[1]
cat("\n===== ROCKY - INDURATED: ordered beta model comparison =====\n")
print(pairwise_df)

comparison_table <- data.frame(
  Response = RESPONSE_COL,
  Family = "Ordered beta (logit)",
  Comparison = as.character(pairwise_df$contrast),
  Link_scale_estimate = pairwise_df$estimate,
  SE_link_scale = pairwise_df$SE,
  Test_statistic = pairwise_df[[stat_name]],
  Statistic_type = stat_name,
  p_value = p_pair,
  Adjustment = "Tukey (single comparison)",
  row.names = NULL
)
write.csv(comparison_table,
          file.path(OUT_DIR, "Table_S1_Ordered_Beta_Habitat_Contrast.csv"),
          row.names = FALSE)


# Figure data -------------------------------------------------------------

# Compute an observed mean for each independent REEF SITE first.
# Summarise across sites, thereby avoiding treating video replicates as
# independent reef-level replicates in the descriptive SE.
site_means <- stats::aggregate(fpi ~ habitat + site, dat, mean)
names(site_means)[names(site_means) == "fpi"] <- "site_mean"
write.csv(site_means, file.path(OUT_DIR, "Observed_Per_Site_Means.csv"),
          row.names = FALSE)

by_habitat <- split(site_means$site_mean, site_means$habitat)
plot_df <- data.frame(
  habitat = factor(names(by_habitat), levels = levels(dat$habitat)),
  mean = vapply(by_habitat, mean, numeric(1)),
  SE = vapply(by_habitat, function(v) {
    if (length(v) < 2) return(NA_real_)
    stats::sd(v) / sqrt(length(v))
  }, numeric(1)),
  n_sites = lengths(by_habitat),
  row.names = NULL
)
plot_df <- plot_df[order(plot_df$habitat), , drop = FALSE]
if (is.finite(p_pair) && p_pair < ALPHA) {
  plot_df$letter <- ifelse(plot_df$mean == max(plot_df$mean), "a", "b")
} else {
  plot_df$letter <- "a"
}

upper <- plot_df$mean + ifelse(is.na(plot_df$SE), 0, plot_df$SE)
y_top <- max(0.02, upper, na.rm = TRUE) * 1.35
plot_df$letter_y <- upper + y_top * 0.035
write.csv(plot_df, file.path(OUT_DIR, "Figure_Plotted_Values_SITE_MEANS.csv"),
          row.names = FALSE)

fig <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = habitat, y = mean, fill = habitat)) +
  ggplot2::geom_col(width = 0.48) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, mean - SE), ymax = mean + SE),
    width = 0.12, linewidth = 0.55, color = "black", na.rm = TRUE
  ) +
  ggplot2::geom_text(ggplot2::aes(y = letter_y, label = letter),
                     size = 4.3, color = "black") +
  ggplot2::annotate("text", x = 0.58, y = y_top * 0.96,
                    label = "c)", hjust = 0, fontface = "bold", size = 5) +
  ggplot2::scale_fill_manual(values = c(
    "Rocky reef" = "#0072B2", "Indurated reef" = "#D55E00")) +
  ggplot2::scale_y_continuous(limits = c(0, y_top),
                              expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(x = NULL, y = "Fishing pressure index") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(colour = "black", size = 10),
    axis.text.y = ggplot2::element_text(colour = "black"),
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.line = ggplot2::element_line(linewidth = 0.5),
    plot.margin = ggplot2::margin(12, 12, 8, 8)
  )
print(fig)
ggplot2::ggsave(file.path(OUT_DIR, "Figure_Fishing_Pressure_Habitats.png"),
                fig, width = 4.5, height = 4.0, dpi = 320, bg = "white")
ggplot2::ggsave(file.path(OUT_DIR, "Figure_Fishing_Pressure_Habitats.pdf"),
                fig, width = 4.5, height = 4.0, bg = "white")


# DHARMa diagnostics ------------------------------------------------------

set.seed(SEED)
sim <- DHARMa::simulateResiduals(fittedModel = fit, n = N_SIM, plot = FALSE)
grDevices::png(file.path(OUT_DIR, "DHARMa_Diagnostics.png"),
               width = 1800, height = 1000, res = 175)
plot(sim)
grDevices::dev.off()

report <- capture.output({
  cat("== Response and site counts ==\n")
  print(summary(dat$fpi))
  print(nsites)
  cat("\n== Ordered beta model ==\n")
  print(summary(fit))
  cat("\n== Link-scale habitat contrast ==\n")
  print(comparison_table)
  cat("\n== Figure means (observed, site-level) ==\n")
  print(plot_df)
  cat("\n== DHARMa uniformity ==\n")
  print(DHARMa::testUniformity(sim))
  cat("\n== DHARMa dispersion ==\n")
  print(DHARMa::testDispersion(sim))
  cat("\n== DHARMa outliers ==\n")
  print(DHARMa::testOutliers(sim))
  cat("\n== Session info ==\n")
  print(sessionInfo())
})
writeLines(report, file.path(OUT_DIR, "Model_Results_and_Diagnostics.txt"))

cat("\nComplete. Outputs saved to: ", OUT_DIR, "\n", sep = "")
cat("NOTE: Figure bars are observed site-level means +/- between-site SE; ",
    "letters use the GLMM habitat test.\n", sep = "")
