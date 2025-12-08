library(tidyverse)
library(pracma)

# --- 1. Load Data ---
if (!all(c("spectra_mat", "wavelengths", "mean_spectra", "params") %in% ls())) {
  if (file.exists("app/data/prosail_data.RData")) {
    load("app/data/prosail_data.RData", envir = globalenv())
    message("Loaded PROSAIL reference data.")
  } else {
    stop("prosail_data.RData not found. Please run generate_prosail_data.R first.")
  }
}

wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_400_1100 <- wavelengths[wl_mask]
spectra_mat <- spectra_mat[, wl_mask, drop = FALSE]

# --- 2. Helper Functions ---

# Guyot & Baret Formula for Inflection Point
calc_lambda_i <- function(r670, r700, r740, r780) {
  rho_i <- (r670 + r780) / 2
  denom <- r740 - r700
  if (abs(denom) < 1e-6) denom <- ifelse(denom < 0, -1, 1) * 1e-6
  700 + 40 * ((rho_i - r700) / denom)
}

# Logistic Function
predict_logistic <- function(lambda, C, r_min, r_max, lambda_i) {
  r_min + (r_max - r_min) / (1 + exp(C * (lambda_i - lambda)))
}

# --- 3. Main Reconstruction Function ---

reconstruct_spec_parametric <- function(obs_r, wl) {
  # Define required anchor bands
  # We need: 400 (base), 550 (green), 670 (red min), 700, 740, 750 (base end), 780 (nir start), 980 (water), 1100 (nir end)
  target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)

  # Find indices (nearest neighbor if exact match missing)
  idx <- vapply(target_bands, function(t) which.min(abs(wl - t)), integer(1))
  r_anchors <- obs_r[idx]
  names(r_anchors) <- paste0("R", target_bands)

  # --- Segment 1: Visible (400 - 680 nm) ---
  # Model: Inverted Gaussian on Linear Baseline
  # Baseline connects R400 and R750 (as per instructions)

  slope_base <- (r_anchors["R750"] - r_anchors["R400"]) / (750 - 400)
  intercept_base <- r_anchors["R400"] - slope_base * 400
  calc_baseline <- function(l) slope_base * l + intercept_base

  mu <- 670 # Fixed center of absorption

  # Calculate Amplitude (A) at 670nm
  # R_obs(670) = Baseline(670) - A  =>  A = Baseline(670) - R_obs(670)
  base_670 <- calc_baseline(670)
  A <- base_670 - r_anchors["R670"]

  # Calculate Sigma at 550nm
  # R_obs(550) = Baseline(550) - A * exp(...)
  base_550 <- calc_baseline(550)
  delta_550 <- base_550 - r_anchors["R550"]

  # Solve for sigma: ratio = delta / A = exp(...)
  ratio <- delta_550 / A

  # Safety checks for log
  if (is.na(ratio) || ratio <= 0 || ratio >= 1) {
    sigma <- 40 # Fallback width
  } else {
    # ratio = exp( - (550-670)^2 / 2sigma^2 )
    # ln(ratio) = - (120)^2 / 2sigma^2
    # sigma^2 = - (120)^2 / (2 * ln(ratio))
    sigma <- sqrt(-(550 - 670)^2 / (2 * log(ratio)))
  }

  vis_fun <- function(l) {
    calc_baseline(l) - A * exp(-(l - mu)^2 / (2 * sigma^2))
  }

  # --- Segment 2: Red Edge (680 - 780 nm) ---
  # Model: 4-Parameter Logistic

  lambda_i <- calc_lambda_i(r_anchors["R670"], r_anchors["R700"], r_anchors["R740"], r_anchors["R780"])

  # Optimize C to fit R700 and R740
  rss_logistic <- function(C) {
    preds <- predict_logistic(c(700, 740), C, r_anchors["R670"], r_anchors["R780"], lambda_i)
    sum((preds - c(r_anchors["R700"], r_anchors["R740"]))^2)
  }

  opt_res <- optim(0.05, rss_logistic, method = "L-BFGS-B", lower = 0.001, upper = 2.0)
  C_opt <- opt_res$par

  re_fun <- function(l) {
    predict_logistic(l, C_opt, r_anchors["R670"], r_anchors["R780"], lambda_i)
  }

  # --- Segment 3: NIR Plateau (780 - 1100 nm) ---
  # Model: Quadratic Polynomial (ax^2 + bx + c)
  # Constrained by R780, R980, R1100

  nir_x <- c(780, 980, 1100)
  nir_y <- c(r_anchors["R780"], r_anchors["R980"], r_anchors["R1100"])

  # Solve system of linear equations for a, b, c
  X_mat <- cbind(nir_x^2, nir_x, 1)
  coeffs <- tryCatch(solve(X_mat, nir_y), error = function(e) c(0, 0, mean(nir_y)))
  a_nir <- coeffs[1]
  b_nir <- coeffs[2]
  c_nir <- coeffs[3]

  nir_fun <- function(l) {
    a_nir * l^2 + b_nir * l + c_nir
  }

  # --- Feature: AUC (780 - 1000 nm) ---
  # Integral of quadratic: (a/3)x^3 + (b/2)x^2 + cx
  auc_integ <- function(l) (a_nir / 3) * l^3 + (b_nir / 2) * l^2 + c_nir * l
  auc_val <- auc_integ(1000) - auc_integ(780)

  # --- Construct Full Spectrum ---
  rec_spec <- numeric(length(wl))

  mask_vis <- wl <= 680
  rec_spec[mask_vis] <- vis_fun(wl[mask_vis])

  mask_re <- wl > 680 & wl <= 780
  rec_spec[mask_re] <- re_fun(wl[mask_re])

  mask_nir <- wl > 780
  rec_spec[mask_nir] <- nir_fun(wl[mask_nir])

  list(
    spectrum = rec_spec,
    features = list(
      A = A,
      sigma = sigma,
      lambda_i = lambda_i,
      C = C_opt,
      a_nir = a_nir,
      AUC = auc_val
    )
  )
}

# --- 4. Execution & Feature Extraction ---

message("Processing samples for feature extraction...")
n_proc <- nrow(spectra_mat)
features_list <- vector("list", n_proc)
rmse_vals <- numeric(n_proc)

# Loop through all samples
for (i in 1:n_proc) {
  res <- reconstruct_spec_parametric(spectra_mat[i, ], wl_400_1100)
  features_list[[i]] <- res$features
  rmse_vals[i] <- sqrt(mean((res$spectrum - spectra_mat[i, ])^2))
}

features_df <- bind_rows(features_list)
features_df$RMSE <- rmse_vals

message(sprintf("Average Reconstruction RMSE: %.4f", mean(rmse_vals)))

# Combine with Ground Truth Parameters
analysis_df <- bind_cols(features_df, params[1:n_proc, ])

# --- 5. Validation Models ---

message("\n--- Feature Validation Models ---")

# Model 1: Predict Chlorophyll (Cab)
# Hypothesis: A (Amplitude) and C (Slope) should be strong predictors
model_cab <- lm(Cab ~ A + C + lambda_i + AUC, data = analysis_df)
message("\nPredicting Chlorophyll (Cab):")
print(summary(model_cab)$coefficients)
message(sprintf("R-squared: %.4f", summary(model_cab)$r.squared))

# Model 2: Predict LAI
# Hypothesis: AUC and C should be strong predictors
model_lai <- lm(LAI ~ A + C + lambda_i + AUC, data = analysis_df)
message("\nPredicting LAI:")
print(summary(model_lai)$coefficients)
message(sprintf("R-squared: %.4f", summary(model_lai)$r.squared))

# --- 6. Visualization (Sample) ---
sample_idx <- 100
res_sample <- reconstruct_spec_parametric(spectra_mat[sample_idx, ], wl_400_1100)

plot_data <- tibble(
  WL = wl_400_1100,
  Original = spectra_mat[sample_idx, ],
  Reconstructed = res_sample$spectrum
) %>% pivot_longer(-WL, names_to = "Type", values_to = "Reflectance")

p <- ggplot(plot_data, aes(WL, Reflectance, color = Type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Original" = "darkgreen", "Reconstructed" = "red")) +
  geom_vline(xintercept = c(680, 780), linetype = "dashed", color = "gray") +
  labs(
    title = "Parametric Reconstruction (Feature-Based)",
    subtitle = sprintf(
      "Sample %d | Cab: %.1f, LAI: %.1f | Rec RMSE: %.4f",
      sample_idx, params$Cab[sample_idx], params$LAI[sample_idx], rmse_vals[sample_idx]
    ),
    x = "Wavelength (nm)", y = "Reflectance"
  ) +
  theme_minimal()

print(p)
