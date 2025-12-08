library(tidyverse)
library(pracma)

source("R/reconstruction_model.R")
source("R/error_analysis.R")

if (!file.exists("data/prosail_data.RData")) {
    stop("data/prosail_data.RData not found.")
}
load("data/prosail_data.RData", envir = globalenv())

# --- Run Reconstruction on All Samples ---
wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_grid <- wavelengths[wl_mask]
spectra_mat <- spectra_mat[, wl_mask]

n_samples <- nrow(spectra_mat)
rec_mat <- matrix(NA, nrow = n_samples, ncol = length(wl_grid))

message("Running reconstruction on ", n_samples, " samples...")
target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)

for (i in 1:n_samples) {
    obs <- spectra_mat[i, ]
    idx <- vapply(target_bands, function(t) which.min(abs(wl_grid - t)), integer(1))
    obs_anchors <- obs[idx]
    names(obs_anchors) <- paste0("R", target_bands)

    res <- reconstruct_spectrum(obs_anchors, wl_grid)
    rec_mat[i, ] <- res$spectrum
}

# --- Analyze Errors by Range ---
error_df <- calculate_band_error(spectra_mat, rec_mat, wl_grid)

# Define Ranges
error_df <- error_df %>%
    mutate(Region = case_when(
        Wavelength <= 680 ~ "Visible",
        Wavelength > 680 & Wavelength <= 780 ~ "Red Edge",
        Wavelength > 780 ~ "NIR"
    ))

# Summarize RMSE by Region
region_summary <- error_df %>%
    group_by(Region) %>%
    summarise(
        Mean_RMSE = mean(RMSE),
        Max_RMSE = max(RMSE)
    )

print(region_summary)

# Plot RMSE
png("output/error_analysis_ranges.png", width = 800, height = 600)
ggplot(error_df, aes(x = Wavelength, y = RMSE, color = Region)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = c(680, 780), linetype = "dashed", color = "gray") +
    labs(title = "Reconstruction Error by Wavelength", y = "RMSE", x = "Wavelength (nm)") +
    theme_minimal()
dev.off()

message("Error analysis plot saved to output/error_analysis_ranges.png")
