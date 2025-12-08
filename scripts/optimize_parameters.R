library(tidyverse)
library(pracma)

source("R/reconstruction_model.R")
source("R/cost_functions.R")

if (!file.exists("data/prosail_data.RData")) {
    stop("data/prosail_data.RData not found.")
}
load("data/prosail_data.RData", envir = globalenv())

set.seed(123)
train_idx <- sample(seq_len(nrow(spectra_mat)), 50)
train_obs <- spectra_mat[train_idx, ]
wl_mask <- wavelengths >= 400 & wavelengths <= 1100
wl_grid <- wavelengths[wl_mask]
train_obs <- train_obs[, wl_mask]

optimize_sample <- function(obs_full, wl) {
    target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
    idx <- vapply(target_bands, function(t) which.min(abs(wl - t)), integer(1))
    obs_anchors <- obs_full[idx]
    names(obs_anchors) <- paste0("R", target_bands)

    cost_fn <- function(p) {
        params <- list(
            sigma = p[1],
            C = p[2],
            A_blue = p[3],
            sigma_water = p[4]
        )

        res <- reconstruct_spectrum(obs_anchors, wl, params)
        rec_spec <- res$spectrum

        # DEBUG
        if (!is.numeric(rec_spec) || !is.numeric(obs_full)) {
            message("DEBUG: Non-numeric detected")
            message("rec_spec class: ", class(rec_spec))
            message("obs_full class: ", class(obs_full))
            if (is.list(rec_spec)) message("rec_spec is a list")
        }

        sse <- sum((rec_spec - obs_full)^2)
        return(sse)
    }

    init_par <- c(40, 0.5, 0.01, 40)

    opt <- optim(init_par, cost_fn,
        method = "L-BFGS-B",
        lower = c(10, 0.01, 0, 10), upper = c(100, 5.0, 0.2, 100)
    )

    return(opt$par)
}

test_idx <- 100
obs_test <- spectra_mat[test_idx, wl_mask]

# Ensure obs_test is numeric vector
obs_test <- as.numeric(obs_test)

message("Optimizing Sample ", test_idx, "...")
start_time <- Sys.time()
best_pars <- optimize_sample(obs_test, wl_grid)
end_time <- Sys.time()

message("Optimization finished in ", round(end_time - start_time, 2), "s")
print(best_pars)

message("Best Sigma (Red): ", round(best_pars[1], 2))
message("Best C (Red Edge): ", round(best_pars[2], 4))
message("Best A_blue: ", round(best_pars[3], 4))
message("Best Sigma (Water): ", round(best_pars[4], 2))

target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
idx <- vapply(target_bands, function(t) which.min(abs(wl_grid - t)), integer(1))
obs_anchors <- obs_test[idx]
names(obs_anchors) <- paste0("R", target_bands)

res_default <- reconstruct_spectrum(obs_anchors, wl_grid)
res_opt <- reconstruct_spectrum(obs_anchors, wl_grid, params = list(
    sigma = best_pars[1],
    C = best_pars[2],
    A_blue = best_pars[3],
    sigma_water = best_pars[4]
))

sse_default <- sum((res_default$spectrum - obs_test)^2)
sse_opt <- sum((res_opt$spectrum - obs_test)^2)

message("SSE Default: ", round(sse_default, 4))
message("SSE Optimized: ", round(sse_opt, 4))
message("Improvement: ", round((sse_default - sse_opt) / sse_default * 100, 2), "%")

png("output/optimization_comparison.png", width = 800, height = 600)
plot(wl_grid, obs_test, type = "l", col = "black", lwd = 3, main = "Optimization Effect", ylab = "Reflectance")
lines(wl_grid, res_default$spectrum, col = "red", lty = 2, lwd = 2)
lines(wl_grid, res_opt$spectrum, col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = c("Original", "Default", "Optimized"), col = c("black", "red", "blue"), lty = c(1, 2, 2), lwd = c(3, 2, 2))
dev.off()
message("Comparison plot saved to output/optimization_comparison.png")
