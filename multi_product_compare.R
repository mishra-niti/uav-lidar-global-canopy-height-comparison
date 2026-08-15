# ============================================================
# Compare UAV LiDAR CHM with global canopy height products
# Version: 1.0.0
# Author: Niti B. Mishra
# Products:
#   1. GCH 2020
#   2. GFCH 2019
#
# UAV reference:
#   95th percentile of valid UAV CHM pixels inside each
#   global product pixel.
#
# Error convention:
#   error = Global product height - UAV P95
#
# Interpretation:
#   Negative error = global product underestimates canopy height
#                    relative to UAV LiDAR P95.
#
#   Positive error = global product overestimates canopy height
#                    relative to UAV LiDAR P95.
#
# Main outputs:
#   1. Pixel-level comparison CSV for each product
#   2. Summary statistics CSV for each product
#   3. Histogram of product - UAV P95 for each product
#   4. Combined summary table
#   5. Combined histogram figure
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

# Install these once if needed:
# install.packages(c("terra", "dplyr", "ggplot2", "readr"))

library(terra)
library(dplyr)
library(ggplot2)
library(readr)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------
# In R, Windows paths are easier to write with forward slashes.

uav_chm_path <- "data/uav_chm.tif"

gch_path <- "data/gch_2020.tif"

gfch_path <- "data/gfch_2019.tif"

# Stop early with a clear message if any input is missing.
# Change the three paths above if your files use different names or locations.
input_paths <- c(
  UAV_CHM = uav_chm_path,
  GCH_2020 = gch_path,
  GFCH_2019 = gfch_path
)

missing_inputs <- input_paths[!file.exists(input_paths)]

if (length(missing_inputs) > 0) {
  stop(
    paste0(
      "Input file(s) not found:\n  ",
      paste(missing_inputs, collapse = "\n  "),
      "\n\nPlace the rasters in data/ or edit the paths in section 2."
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 3. Output folder
# ------------------------------------------------------------

out_dir <- "outputs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 4. Read UAV LiDAR CHM
# ------------------------------------------------------------

uav_chm <- rast(uav_chm_path)

names(uav_chm) <- "uav_chm"

# Basic CHM cleanup.
# Negative CHM values usually represent artifacts.
# The upper limit can be modified if needed.
uav_chm[uav_chm < 0] <- NA
uav_chm[uav_chm > 80] <- NA

cat("\nUAV CHM information:\n")
print(uav_chm)


# ------------------------------------------------------------
# 5. Function to compare one global product with UAV CHM
# ------------------------------------------------------------

compare_product_to_uav <- function(
    product_path,
    product_name,
    uav_chm,
    out_dir,
    min_uav_pixels = 10
) {
  
  cat("\n============================================================\n")
  cat("Processing:", product_name, "\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # 5.1 Read global canopy height product
  # ----------------------------------------------------------
  
  product <- rast(product_path)
  
  # Make sure product has a simple name
  names(product) <- "product_height"
  
  # Remove invalid or unreasonable values.
  # Many products use 0 or NA for no-data / no-canopy areas.
  product[product <= 0] <- NA
  product[product > 80] <- NA
  
  cat("\nProduct raster information:\n")
  print(product)
  
  
  # ----------------------------------------------------------
  # 5.2 Check CRS and reproject if needed
  # ----------------------------------------------------------
  #
  # The UAV CHM and global product must be in the same CRS
  # for the analysis.
  #
  # If the CRS differs, this reprojects the product to the UAV CRS.
  # Nearest-neighbor is used to avoid smoothing the product heights.
  
  if (!same.crs(product, uav_chm)) {
    
    cat("\nCRS differs. Projecting product to UAV CHM CRS...\n")
    
    product <- project(
      product,
      crs(uav_chm),
      method = "near"
    )
    
    names(product) <- "product_height"
  }
  
  
  # ----------------------------------------------------------
  # 5.3 Crop UAV CHM to the product extent
  # ----------------------------------------------------------
  #
  # We crop the UAV CHM to the global product extent.
  #
  # IMPORTANT:
  # Do not directly run:
  #
  #   mask(uav_crop, product)
  #
  # because the two rasters may overlap spatially but still have
  # different grid geometry. That causes:
  #
  #   Error: [mask] extents do not match
  #
  # Instead, we resample the product to the UAV grid first.
  
  uav_crop <- crop(uav_chm, product, snap = "out")
  
  
  # ----------------------------------------------------------
  # 5.4 Resample product to UAV grid for masking only
  # ----------------------------------------------------------
  #
  # This raster is used only to identify UAV pixels that overlap
  # valid global product pixels.
  
  product_on_uav_grid <- resample(
    product,
    uav_crop,
    method = "near"
  )
  
  names(product_on_uav_grid) <- "product_height_on_uav_grid"
  
  # Keep UAV CHM pixels only where the global product has valid data.
  # Now mask works because product_on_uav_grid and uav_crop have
  # identical grid geometry.
  uav_crop <- mask(uav_crop, product_on_uav_grid)
  
  n_uav_pixels <- global(!is.na(uav_crop), "sum", na.rm = TRUE)[1, 1]
  
  cat("\nValid UAV pixels after crop/mask:", n_uav_pixels, "\n")
  
  
  # ----------------------------------------------------------
  # 5.5 Create a unique cell ID for every product pixel
  # ----------------------------------------------------------
  #
  # Each product pixel gets a unique ID.
  # We will use this ID to group UAV pixels by product pixel.
  
  product_cell_id <- product
  
  values(product_cell_id) <- seq_len(ncell(product_cell_id))
  
  # Keep IDs only where product height is valid.
  product_cell_id <- mask(product_cell_id, product)
  
  names(product_cell_id) <- "cell_id"
  
  
  # ----------------------------------------------------------
  # 5.6 Resample product cell IDs to the UAV grid
  # ----------------------------------------------------------
  #
  # Cell IDs are categorical, so use nearest neighbor.
  # Do not use bilinear interpolation for ID rasters.
  
  cell_id_on_uav_grid <- resample(
    product_cell_id,
    uav_crop,
    method = "near"
  )
  
  names(cell_id_on_uav_grid) <- "cell_id"
  
  
  # ----------------------------------------------------------
  # 5.7 Calculate UAV P95 inside each product pixel
  # ----------------------------------------------------------
  #
  # zonal() groups UAV CHM pixels by product-pixel ID.
  #
  # IMPORTANT:
  # The custom function must include "..." because terra::zonal()
  # may pass additional arguments such as na.rm.
  #
  # Without "...", you may get:
  #
  #   Error in FUN(X[[i]], ...) : unused argument (na.rm = TRUE)
  
  cat("\nCalculating UAV P95 inside each", product_name, "pixel...\n")
  
  uav_p95_by_cell <- zonal(
    uav_crop,
    cell_id_on_uav_grid,
    fun = function(x, ...) {
      quantile(x, probs = 0.95, na.rm = TRUE)
    },
    na.rm = TRUE
  )
  
  names(uav_p95_by_cell) <- c("cell_id", "uav_p95")
  
  
  # ----------------------------------------------------------
  # 5.8 Count valid UAV pixels inside each product pixel
  # ----------------------------------------------------------
  #
  # This is a quality-control variable.
  # Product pixels with too few UAV CHM pixels should be excluded,
  # especially near the UAV CHM edge.
  
  uav_n_by_cell <- zonal(
    uav_crop,
    cell_id_on_uav_grid,
    fun = function(x, ...) {
      sum(!is.na(x))
    },
    na.rm = FALSE
  )
  
  names(uav_n_by_cell) <- c("cell_id", "uav_n")
  
  
  # ----------------------------------------------------------
  # 5.9 Extract global product height for each product pixel
  # ----------------------------------------------------------
  
  product_df <- as.data.frame(
    product,
    cells = TRUE,
    na.rm = TRUE
  )
  
  product_df <- product_df %>%
    rename(
      cell_id = cell,
      product_height_m = product_height
    )
  
  
  # ----------------------------------------------------------
  # 5.10 Combine product height with UAV P95
  # ----------------------------------------------------------
  
  comparison_df <- product_df %>%
    left_join(uav_p95_by_cell, by = "cell_id") %>%
    left_join(uav_n_by_cell, by = "cell_id") %>%
    filter(
      !is.na(product_height_m),
      !is.na(uav_p95),
      uav_n >= min_uav_pixels
    ) %>%
    mutate(
      product_name = product_name,
      
      # Main error metric
      error_product_minus_uav_p95 = product_height_m - uav_p95,
      
      # Absolute error
      abs_error_m = abs(error_product_minus_uav_p95)
    )
  
  cat("\nNumber of product pixels retained:", nrow(comparison_df), "\n")
  cat("Minimum UAV pixels required per product pixel:", min_uav_pixels, "\n")
  
  
  # ----------------------------------------------------------
  # 5.11 Summary statistics
  # ----------------------------------------------------------
  
  summary_df <- comparison_df %>%
    summarize(
      product_name = first(product_name),
      n_pixels = n(),
      
      product_median_m = median(product_height_m, na.rm = TRUE),
      uav_p95_median_m = median(uav_p95, na.rm = TRUE),
      
      mean_bias_m = mean(error_product_minus_uav_p95, na.rm = TRUE),
      median_bias_m = median(error_product_minus_uav_p95, na.rm = TRUE),
      
      mae_m = mean(abs_error_m, na.rm = TRUE),
      
      rmse_m = sqrt(
        mean(error_product_minus_uav_p95^2, na.rm = TRUE)
      ),
      
      error_p05_m = as.numeric(
        quantile(error_product_minus_uav_p95, 0.05, na.rm = TRUE)
      ),
      
      error_p25_m = as.numeric(
        quantile(error_product_minus_uav_p95, 0.25, na.rm = TRUE)
      ),
      
      error_p75_m = as.numeric(
        quantile(error_product_minus_uav_p95, 0.75, na.rm = TRUE)
      ),
      
      error_p95_m = as.numeric(
        quantile(error_product_minus_uav_p95, 0.95, na.rm = TRUE)
      ),
      
      pearson_r = cor(
        product_height_m,
        uav_p95,
        use = "complete.obs"
      )
    )
  
  cat("\nSummary statistics:\n")
  print(summary_df)
  
  
  # ----------------------------------------------------------
  # 5.12 Save comparison and summary tables
  # ----------------------------------------------------------
  
  comparison_csv <- file.path(
    out_dir,
    paste0(product_name, "_vs_UAV_P95_pixel_comparison.csv")
  )
  
  summary_csv <- file.path(
    out_dir,
    paste0(product_name, "_vs_UAV_P95_summary.csv")
  )
  
  write_csv(comparison_df, comparison_csv)
  write_csv(summary_df, summary_csv)
  
  cat("\nSaved pixel comparison table to:\n", comparison_csv, "\n")
  cat("Saved summary table to:\n", summary_csv, "\n")
  
  
  # ----------------------------------------------------------
  # 5.13 Histogram of product minus UAV P95
  # ----------------------------------------------------------
  #
  # This follows the style of the example figure.
  #
  # x-axis:
  #   product height - UAV P95
  #
  # y-axis:
  #   density
  #
  # dashed vertical line:
  #   median error
  
  median_error <- median(
    comparison_df$error_product_minus_uav_p95,
    na.rm = TRUE
  )
  
  hist_plot <- ggplot(
    comparison_df,
    aes(x = error_product_minus_uav_p95)
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 50,
      fill = "grey70",
      color = "grey70"
    ) +
    geom_vline(
      xintercept = median_error,
      linetype = "dashed",
      linewidth = 0.9
    ) +
    labs(
      title = paste0(product_name, " - UAV P95"),
      subtitle = paste0(
        "Median error = ", round(median_error, 2), " m; ",
        "n = ", nrow(comparison_df), " product pixels"
      ),
      x = paste0(product_name, " - UAV P95 (m)"),
      y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  print(hist_plot)
  
  hist_png <- file.path(
    out_dir,
    paste0(product_name, "_minus_UAV_P95_histogram.png")
  )
  
  hist_pdf <- file.path(
    out_dir,
    paste0(product_name, "_minus_UAV_P95_histogram.pdf")
  )
  
  ggsave(
    filename = hist_png,
    plot = hist_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = hist_pdf,
    plot = hist_plot,
    width = 7,
    height = 5
  )
  
  cat("\nSaved histogram PNG to:\n", hist_png, "\n")
  cat("Saved histogram PDF to:\n", hist_pdf, "\n")
  
  
  # ----------------------------------------------------------
  # 5.14 Return results to R environment
  # ----------------------------------------------------------
  
  return(
    list(
      comparison = comparison_df,
      summary = summary_df,
      histogram = hist_plot
    )
  )
}


# ------------------------------------------------------------
# 6. Run comparison for GCH 2020
# ------------------------------------------------------------

gch_results <- compare_product_to_uav(
  product_path = gch_path,
  product_name = "GCH_2020",
  uav_chm = uav_chm,
  out_dir = out_dir,
  min_uav_pixels = 10
)


# ------------------------------------------------------------
# 7. Run comparison for GFCH 2019
# ------------------------------------------------------------

gfch_results <- compare_product_to_uav(
  product_path = gfch_path,
  product_name = "GFCH_2019",
  uav_chm = uav_chm,
  out_dir = out_dir,
  min_uav_pixels = 10
)


# ------------------------------------------------------------
# 8. Combine summaries into one table
# ------------------------------------------------------------

combined_summary <- bind_rows(
  gch_results$summary,
  gfch_results$summary
)

cat("\nCombined summary:\n")
print(combined_summary)

combined_summary_csv <- file.path(
  out_dir,
  "GCH_GFCH_vs_UAV_P95_combined_summary.csv"
)

write_csv(
  combined_summary,
  combined_summary_csv
)

cat("\nSaved combined summary to:\n", combined_summary_csv, "\n")


# ------------------------------------------------------------
# 9. Combine pixel-level comparison tables
# ------------------------------------------------------------

combined_comparison <- bind_rows(
  gch_results$comparison,
  gfch_results$comparison
)

combined_comparison_csv <- file.path(
  out_dir,
  "GCH_GFCH_vs_UAV_P95_combined_pixel_comparison.csv"
)

write_csv(
  combined_comparison,
  combined_comparison_csv
)

cat("\nSaved combined pixel comparison table to:\n", combined_comparison_csv, "\n")


# ------------------------------------------------------------
# 10. Combined histogram for both products
# ------------------------------------------------------------

median_errors <- combined_comparison %>%
  group_by(product_name) %>%
  summarize(
    median_error = median(error_product_minus_uav_p95, na.rm = TRUE),
    n_pixels = n(),
    .groups = "drop"
  )

combined_hist <- ggplot(
  combined_comparison,
  aes(x = error_product_minus_uav_p95)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "grey70",
    color = "grey70"
  ) +
  geom_vline(
    data = median_errors,
    aes(xintercept = median_error),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ product_name,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Global Canopy Height Products Compared with UAV LiDAR P95",
    subtitle = "Error calculated as global product canopy height minus UAV CHM P95",
    x = "Product height - UAV P95 (m)",
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(combined_hist)

combined_hist_png <- file.path(
  out_dir,
  "GCH_GFCH_minus_UAV_P95_combined_histogram.png"
)

combined_hist_pdf <- file.path(
  out_dir,
  "GCH_GFCH_minus_UAV_P95_combined_histogram.pdf"
)

ggsave(
  filename = combined_hist_png,
  plot = combined_hist,
  width = 7,
  height = 8,
  dpi = 300
)

ggsave(
  filename = combined_hist_pdf,
  plot = combined_hist,
  width = 7,
  height = 8
)

cat("\nSaved combined histogram PNG to:\n", combined_hist_png, "\n")
cat("Saved combined histogram PDF to:\n", combined_hist_pdf, "\n")


# ------------------------------------------------------------
# 11. Save R and package version information
# ------------------------------------------------------------

session_info_path <- file.path(
  out_dir,
  "session_info.txt"
)

writeLines(
  capture.output(sessionInfo()),
  con = session_info_path
)

cat("\nSaved R session information to:\n", session_info_path, "\n")


# ------------------------------------------------------------
# 12. Interpretation guide
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("INTERPRETATION GUIDE\n")
cat("============================================================\n")
cat("Main error column:\n")
cat("  error_product_minus_uav_p95 = product_height_m - uav_p95\n\n")
cat("Negative values:\n")
cat("  Global product is lower than UAV LiDAR P95.\n\n")
cat("Positive values:\n")
cat("  Global product is higher than UAV LiDAR P95.\n\n")
cat("Median bias near zero:\n")
cat("  Little systematic bias.\n\n")
cat("Large RMSE or wide histogram:\n")
cat("  Large pixel-level disagreement even if median bias is small.\n")
cat("============================================================\n")


# ============================================================
# End of script
# ============================================================
