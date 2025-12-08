library(tidyverse)
library(pracma)

# --- 1. Load Logic ---
source("R/reconstruction_model.R")
source("R/error_analysis.R")

# --- 2. Load Data ---
if (!file.exists("data/prosail_data.RData")) {
  stop("data/prosail_data.RData not found. Please run R/generate_prosail_data.R first.")
}
load("data/prosail_data.RData", envir = globalenv())
message("Loaded PROSAIL reference data.")

wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_grid <- wavelengths[wl_mask]
spectra_mat <- spectra_mat[, wl_mask, drop = FALSE]

# --- 3. Run Reconstruction on a Sample ---
sample_idx <- 100
obs <- spectra_mat[sample_idx, ]

# Prepare input for the model
target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
idx <- vapply(target_bands, function(t) which.min(abs(wl_grid - t)), integer(1))
obs_values <- obs[idx]
names(obs_values) <- paste0("R", target_bands)

# Reconstruct
res <- reconstruct_spectrum(obs_values, wl_grid)

# --- 4. Visualize ---
plot(wl_grid, obs, type="l", col="darkgreen", lwd=2, 
     main=paste("Sample", sample_idx), xlab="Wavelength (nm)", ylab="Reflectance")
lines(wl_grid, res, col="red", lwd=2, lty=2)
legend("topleft", legend=c("Original", "Reconstructed"), col=c("darkgreen", "red"), lty=c(1, 2), lwd=2)

# --- 5. Error Analysis (Example) ---
# Run on first 10 samples
n_test <- 10
rec_mat <- matrix(NA, nrow=n_test, ncol=length(wl_grid))

for(i in 1:n_test) {
  obs_i <- spectra_mat[i, ]
  idx_i <- vapply(target_bands, function(t) which.min(abs(wl_grid - t)), integer(1))
  vals_i <- obs_i[idx_i]
  names(vals_i) <- paste0("R", target_bands)
  
  rec_mat[i, ] <- reconstruct_spectrum(vals_i, wl_grid)
}

# Calculate error
error_df <- calculate_band_error(spectra_mat[1:n_test, ], rec_mat, wl_grid)
print(head(error_df))

# Find worst regions
worst <- get_max_error_regions(error_df)
print("Regions with highest error:")
print(worst)
