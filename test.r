library(tidyverse)
library(pracma)

# Load data
if (!all(c("spectra_mat", "wavelengths", "mean_spectra", "params") %in% ls())) {
    if (file.exists("prosail_data.RData")) {
        load("prosail_data.RData", envir = globalenv())
        message("Loaded PROSAIL reference data.")
    } else {
        stop("prosail_data.RData not found.")
    }
}

wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_400_1100 <- wavelengths[wl_mask]
spectra_mat <- spectra_mat[, wl_mask, drop = FALSE]

# --- Helper Functions ---

calc_lambda_i <- function(r670, r700, r740, r780) {
    rho_i <- (r670 + r780) / 2
    denom <- r740 - r700
    if (abs(denom) < 1e-6) denom <- ifelse(denom < 0, -1, 1) * 1e-6
    700 + 40 * ((rho_i - r700) / denom)
}

predict_logistic <- function(lambda, C, r_min, r_max, lambda_i) {
    r_min + (r_max - r_min) / (1 + exp(C * (lambda_i - lambda)))
}

# --- Reconstruction & Feature Extraction ---

reconstruct_spec_parametric <- function(obs_r, wl) {
    # Required Anchor Bands
    target_bands <- c(400, 550, 670, 700, 740, 780, 980, 1100)
    idx <- match(target_bands, wl)

    if (any(is.na(idx))) {
        # Fallback: Find nearest indices
        idx <- vapply(target_bands, function(t) which.min(abs(wl - t)), integer(1))
    }

    r_anchors <- obs_r[idx]
    names(r_anchors) <- paste0("R", target_bands)

    # 1. Visible Segment (400-680 nm)
    # Model: R(l) = Baseline(l) - A * exp(-(l-mu)^2 / 2sigma^2)
    # Baseline: Linear between R400 and R780

    slope_base <- (r_anchors["R780"] - r_anchors["R400"]) / (780 - 400)
    intercept_base <- r_anchors["R400"] - slope_base * 400
    calc_baseline <- function(l) slope_base * l + intercept_base

    mu <- 670 # Fixed

    # Calculate A from R670
    base_670 <- calc_baseline(670)
    A <- base_670 - r_anchors["R670"]

    # Calculate sigma from R550
    base_550 <- calc_baseline(550)
    delta_550 <- base_550 - r_anchors["R550"]

    # Avoid log of negative or zero
    ratio <- delta_550 / A
    if (is.na(ratio) || ratio <= 0 || ratio >= 1) {
        sigma <- 50 # Default width if fit fails
    } else {
        sigma <- sqrt(-(550 - 670)^2 / (2 * log(ratio)))
    }

    vis_fun <- function(l) {
        val <- calc_baseline(l) - A * exp(-(l - mu)^2 / (2 * sigma^2))
        return(val)
    }

    # 2. Red Edge Segment (680-780 nm)
    # Model: Logistic

    lambda_i <- calc_lambda_i(r_anchors["R670"], r_anchors["R700"], r_anchors["R740"], r_anchors["R780"])

    rss_logistic <- function(C) {
        preds <- predict_logistic(c(700, 740), C, r_anchors["R670"], r_anchors["R780"], lambda_i)
        sum((preds - c(r_anchors["R700"], r_anchors["R740"]))^2)
    }

    opt_res <- optim(0.05, rss_logistic, method = "L-BFGS-B", lower = 0.001, upper = 1.0)
    C_opt <- opt_res$par

    re_fun <- function(l) {
        predict_logistic(l, C_opt, r_anchors["R670"], r_anchors["R780"], lambda_i)
    }

    # 3. NIR Segment (780-1100 nm)
    # Model: Quadratic a*l^2 + b*l + c
    # Anchors: 780, 980, 1100

    nir_x <- c(780, 980, 1100)
    nir_y <- c(r_anchors["R780"], r_anchors["R980"], r_anchors["R1100"])

    X_mat <- cbind(nir_x^2, nir_x, 1)
    # Solve for [a, b, c]
    coeffs <- tryCatch(solve(X_mat, nir_y), error = function(e) c(0, 0, mean(nir_y)))
    a_nir <- coeffs[1]
    b_nir <- coeffs[2]
    c_nir <- coeffs[3]

    nir_fun <- function(l) {
        a_nir * l^2 + b_nir * l + c_nir
    }

    # AUC (780 - 1000)
    auc_integ <- function(l) (a_nir / 3) * l^3 + (b_nir / 2) * l^2 + c_nir * l
    auc_val <- auc_integ(1000) - auc_integ(780)

    # Construct Spectrum
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

# --- Execution & Validation ---

message("Processing samples for feature extraction...")
n_proc <- nrow(spectra_mat)
features_list <- vector("list", n_proc)
rmse_vals <- numeric(n_proc)

for (i in 1:n_proc) {
    res <- reconstruct_spec_parametric(spectra_mat[i, ], wl_400_1100)
    features_list[[i]] <- res$features
    rmse_vals[i] <- sqrt(mean((res$spectrum - spectra_mat[i, ])^2))
}

features_df <- bind_rows(features_list)
features_df$RMSE <- rmse_vals

message(sprintf("Average Reconstruction RMSE: %.4f", mean(rmse_vals)))

# Combine with Ground Truth
# Ensure params match the spectra (assuming order is preserved from generation)
analysis_df <- bind_cols(features_df, params[1:n_proc, ])

message("\n--- Feature Validation Models ---")

# Model 1: Predict Chlorophyll (Cab)
model_cab <- lm(Cab ~ A + C + lambda_i + AUC, data = analysis_df)
message("\nPredicting Chlorophyll (Cab):")
print(summary(model_cab)$coefficients)
message(sprintf("R-squared: %.4f", summary(model_cab)$r.squared))

# Model 2: Predict LAI
model_lai <- lm(LAI ~ A + C + lambda_i + AUC, data = analysis_df)
message("\nPredicting LAI:")
print(summary(model_lai)$coefficients)
message(sprintf("R-squared: %.4f", summary(model_lai)$r.squared))

# --- Visualization of one sample ---
sample_idx <- 100
res_sample <- reconstruct_spec_parametric(spectra_mat[sample_idx, ], wl_400_1100)

plot_df <- tibble(
    WL = wl_400_1100,
    Original = spectra_mat[sample_idx, ],
    Reconstructed = res_sample$spectrum
) %>% pivot_longer(-WL, names_to = "Type", values_to = "Reflectance")

p <- ggplot(plot_df, aes(WL, Reflectance, color = Type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Original" = "darkgreen", "Reconstructed" = "red")) +
    labs(
        title = "Parametric Reconstruction (Feature-Based)",
        subtitle = sprintf(
            "Sample %d | Cab: %.1f, LAI: %.1f | Rec RMSE: %.4f",
            sample_idx, params$Cab[sample_idx], params$LAI[sample_idx], rmse_vals[sample_idx]
        )
    ) +
    theme_minimal()

print(p)
