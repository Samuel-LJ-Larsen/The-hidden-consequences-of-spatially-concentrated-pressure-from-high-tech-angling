# Predictor Collinearity Diagnostic (Pearson Correlations)

# Purpose:
#   Screens a set of continuous environmental/habitat predictor variables for
#   pairwise collinearity before they go into a downstream model. Produces a
#   full correlation matrix, a lower-triangle heatmap highlighting |r| > 0.7,
#   and a flat list of the offending pairs.

# Required input file (in the working directory):
#   - master csv file
#       One row per observation/replicate. Required columns:
#         habitat
#       Plus whichever of the `want_preds` columns (below) are present;
#       the script only errors if NONE of them are found.
#       User change the habitats relevant to their analysis

# What this script does:
#   1. Load master csv
#   2. Coerce predictors to numeric and drop any row with a missing value in any predictor (complete-case analysis)
#   3. Compute the Pearson correlation matrix on the complete cases.
#   4. Export the full matrix as CSV, and a lower-triangle heatmap (PDF) that visually mutes |r| <= 0.7 to white so only concerning pairs stand out.
#   5. Flag and export every pair with |r| > 0.7 (upper triangle only, so each pair is listed once rather than twice).
#   6. Print summary statistics (median/mean/max |r|, % of pairs over threshold) to the console.

# Output files (written to ./correlation_analysis_outputs/):
#   01_correlation_matrix_full.csv  - full predictor x predictor Pearson matrix
#   01_correlation_matrix_plot.pdf  - lower-triangle heatmap, |r| > 0.7 in red
#   01_high_correlations_list.csv   - pairs with |r| > threshold (only written if at least one such pair exists)
#

# Load packages -----------------------------------------------------------


library(corrplot) 
library(dplyr)

csv_path <- "input" # User's master csv file
if (!file.exists(csv_path)) {
  stop("Error: final.csv not found in current directory")
}

output_dir <- "correlation_analysis_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

dat_raw <- read.csv(csv_path, check.names = FALSE)

# Filter to relevant habitats
habitats_keep <- c("input") # Change for user's relevant habitats
dat_raw <- dat_raw %>% dplyr::filter(habitat %in% habitats_keep)

# Ensure habitat is a factor
dat_raw$habitat <- factor(dat_raw$habitat)

# Define predictor variables (continuous)
want_preds <- c(
"input") # User's predictor variables

# Get only the predictors that exist in the data
env_pred_cols <- intersect(want_preds, names(dat_raw))

if (length(env_pred_cols) == 0) {
  stop("No predictor variables found in the dataset!")
}

cat("DIAGNOSTIC ANALYSIS FOR", length(env_pred_cols), "CONTINUOUS PREDICTORS\n")
cat("Variables included:\n")
for (i in seq_along(env_pred_cols)) {
  cat(sprintf("  %2d. %s\n", i, env_pred_cols[i]))
}
cat("\n")

# Extract predictor data
pred_data <- dat_raw[, env_pred_cols, drop = FALSE]

# Convert to numeric if needed
pred_data[] <- lapply(pred_data, function(x) as.numeric(as.character(x)))

# Remove rows with NA values
pred_data_complete <- pred_data[complete.cases(pred_data), ]

cat("Sample size:", nrow(pred_data_complete), "observations\n")
cat("Complete cases:", nrow(pred_data_complete), "/", nrow(pred_data), "\n\n")

# 1. PEARSON CORRELATION ANALYSIS
cat(rep("=", 80), "\n", sep="")
cat("1. PEARSON CORRELATION ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep="")

# Calculate correlation matrix
cor_matrix <- cor(pred_data_complete, use = "complete.obs")
cor_matrix_rounded <- round(cor_matrix, 3)

# Save full correlation matrix
write.csv(cor_matrix_rounded, 
          file.path(output_dir, "01_correlation_matrix_full.csv"), 
          row.names = TRUE)

# Create custom color palette - heatmap style
# White for values between -0.7 and 0.7, red for values outside this range
n_colors <- 200
col_palette <- colorRampPalette(c("#8B0000", "#FF6666", "white", "white", "#FF6666", "#8B0000"))(n_colors)

# Find which colors correspond to values outside ±0.7
breaks_vec <- seq(-1, 1, length.out = n_colors)
threshold_idx_lower <- max(which(breaks_vec <= -0.7))
threshold_idx_upper <- min(which(breaks_vec >= 0.7))

# Set middle values to white
col_palette[(threshold_idx_lower+1):(threshold_idx_upper-1)] <- "white"

# Create correlation plot as PDF (lower triangle only)
pdf(file.path(output_dir, "01_correlation_matrix_plot.pdf"), 
    width = 14, height = 14)
corrplot(cor_matrix, 
         method = "color", 
         type = "lower",
         diag = FALSE,
         tl.col = "black", 
         tl.srt = 45, 
         tl.cex = 0.55,
         addCoef.col = "black",
         number.cex = 0.35,
         number.font = 1,  # 1 = plain (not bold)
         col = col_palette,  # USE col_palette instead of colors
         is.corr = TRUE,
         title = "Pearson Correlation Matrix (|r| > 0.7 highlighted in red)",
         mar = c(0,0,2,0),
         addgrid.col = "gray80",
         cl.pos = "r")
dev.off()

# Identify high correlations
threshold <- 0.7
upper_tri <- cor_matrix
upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA
high_cor_indices <- which(abs(upper_tri) > threshold, arr.ind = TRUE)

if (nrow(high_cor_indices) > 0) {
  high_cor_list <- data.frame(
    Variable1 = rownames(cor_matrix)[high_cor_indices[, 1]],
    Variable2 = colnames(cor_matrix)[high_cor_indices[, 2]],
    Correlation = cor_matrix[high_cor_indices]
  )
  high_cor_list <- high_cor_list[order(-abs(high_cor_list$Correlation)), ]
  
  write.csv(high_cor_list, 
            file.path(output_dir, "01_high_correlations_list.csv"), 
            row.names = FALSE)
  
  cat("High correlations (|r| >", threshold, "):\n")
  for (i in 1:min(30, nrow(high_cor_list))) {
    cat(sprintf("  %2d. %-45s <-> %-45s  r = %7.3f\n",
                i,
                high_cor_list$Variable1[i],
                high_cor_list$Variable2[i],
                high_cor_list$Correlation[i]))
  }
  if (nrow(high_cor_list) > 30) {
    cat(sprintf("  ... and %d more pairs (see CSV file)\n", nrow(high_cor_list) - 30))
  }
  cat("\n")
} else {
  cat("No pairs with |r| >", threshold, "\n\n")
}

# Summary statistics
abs_cors <- abs(upper_tri[!is.na(upper_tri)])
cat("Summary statistics:\n")
cat(sprintf("  Median |r|: %.3f\n", median(abs_cors)))
cat(sprintf("  Mean |r|:   %.3f\n", mean(abs_cors)))
cat(sprintf("  Max |r|:    %.3f\n", max(abs_cors)))
cat(sprintf("  Pairs with |r| > 0.7: %d (%.1f%%)\n\n", 
            sum(abs_cors > 0.7), 100*mean(abs_cors > 0.7)))

