# PROSAIL reference data generator (400-1100 nm)
# Requires the hsdar package (includes PROSAIL implementation)

suppressPackageStartupMessages({
    library(hsdar)
    library(tidyverse)
})

n_samples <- 500
set.seed(123)

params <- tibble(
    N = runif(n_samples, 1.0, 2.5),
    Cab = runif(n_samples, 10, 80),
    Car = runif(n_samples, 0, 20),
    Cbrown = runif(n_samples, 0, 1),
    Cw = runif(n_samples, 0.005, 0.02),
    Cm = runif(n_samples, 0.005, 0.02),
    LAI = runif(n_samples, 0.5, 8),
    psoil = runif(n_samples, 0, 1),
    typelidf = 1
)

spectra_list <- vector("list", n_samples)
message("Starting PROSAIL simulation (", n_samples, " samples)...")

for (i in seq_len(n_samples)) {
    p <- params[i, ]
    spectra_list[[i]] <- tryCatch(
        {
            spec_obj <- PROSAIL(
                N = p$N,
                Cab = p$Cab,
                Car = p$Car,
                Cbrown = p$Cbrown,
                Cw = p$Cw,
                Cm = p$Cm,
                LAI = p$LAI,
                psoil = p$psoil,
                TypeLidf = p$typelidf
            )

            vals <- if (inherits(spec_obj, "Speclib")) {
                m <- spectra(spec_obj)
                if (is.matrix(m)) m[1, ] else as.numeric(m)
            } else {
                as.numeric(spec_obj)
            }

            vals
        },
        error = function(e) {
            message(sprintf("Sample %d failed: %s", i, e$message))
            NULL
        }
    )
}

valid <- !sapply(spectra_list, is.null)
params <- params[valid, ]
spectra_list <- spectra_list[valid]

if (!length(spectra_list)) {
    stop("No successful PROSAIL simulations. Check hsdar installation.")
}

spectra_full <- do.call(rbind, spectra_list)
full_wavelengths <- 400:2500
mask <- full_wavelengths >= 400 & full_wavelengths <= 1100

spectra_mat <- spectra_full[, mask, drop = FALSE]
wavelengths <- full_wavelengths[mask]
mean_spectra <- colMeans(spectra_mat, na.rm = TRUE)

save(spectra_mat, wavelengths, mean_spectra, params, file = "prosail_data.RData")
message("Saved PROSAIL reference set to prosail_data.RData (", nrow(spectra_mat), " samples).")
