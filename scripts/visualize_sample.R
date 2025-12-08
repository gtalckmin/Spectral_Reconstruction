library(tidyverse)
library(pracma)

# --- Load Logic ---
source("R/reconstruction_model.R")
source("R/error_analysis.R")

# --- Load Data ---
if (!file.exists("data/prosail_data.RData")) {
    stop("data/prosail_data.RData not found.")
}
load("data/prosail_data.RData", envir = globalenv())

# --- Configuration ---
# Select a random sample or specify one
sample_idx <- sample(1:nrow(spectra_mat), 1)
# sample_idx <- 100 # Uncomment to fix sample

message(paste("Visualizing Sample:", sample_idx))

# --- Reconstruction ---
obs <- as.numeric(spectra_mat[sample_idx, ])
wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_grid <- wavelengths[wl_mask]
obs <- obs[wl_mask]

# Prepare input
target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
idx <- vapply(target_bands, function(t) which.min(abs(wl_grid - t)), integer(1))
obs_values <- obs[idx]
names(obs_values) <- paste0("R", target_bands)

# Run Model
res <- reconstruct_spectrum(obs_values, wl_grid)

# --- Plotting ---
png("output/sample_visualization.png", width = 800, height = 600)

# Layout: Top = Spectrum, Bottom = Error
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# 1. Spectrum
plot(wl_grid, obs,
    type = "l", col = "darkgreen", lwd = 3,
    main = paste("Reconstruction of Sample", sample_idx),
    xlab = "Wavelength (nm)", ylab = "Reflectance", ylim = c(0, max(obs, res$spectrum) * 1.1)
)
lines(wl_grid, res$spectrum, col = "red", lwd = 2, lty = 2)

# Add anchor points
points(wl_grid[idx], obs[idx], col = "blue", pch = 19, cex = 1.2)
legend("topleft",
    legend = c("Original (PROSAIL)", "Reconstructed", "Anchor Bands"),
    col = c("darkgreen", "red", "blue"), lty = c(1, 2, NA), pch = c(NA, NA, 19), lwd = c(3, 2, NA)
)

# Add regions
abline(v = c(680, 780), lty = 3, col = "gray")
text(c(540, 730, 940), c(0.02, 0.02, 0.02), c("Visible", "Red Edge", "NIR"), col = "gray40")

# 2. Residuals (Error)
residuals <- res$spectrum - obs
plot(wl_grid, residuals,
    type = "h", col = "purple", lwd = 2,
    main = "Residuals (Reconstructed - Original)",
    xlab = "Wavelength (nm)", ylab = "Error"
)
abline(h = 0, col = "black")

# Highlight max error
max_err_idx <- which.max(abs(residuals))
points(wl_grid[max_err_idx], residuals[max_err_idx], col = "red", pch = 8, cex = 1.5)
text(wl_grid[max_err_idx], residuals[max_err_idx],
    labels = paste0("Max Err: ", round(residuals[max_err_idx], 4), " @ ", wl_grid[max_err_idx], "nm"),
    pos = 3, col = "red"
)

dev.off()

message("Visualization saved to output/sample_visualization.png")
a
