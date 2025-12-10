suppressPackageStartupMessages({
    library(hsdar)
    library(e1071)
    library(tibble)
})

# Method 1 anchor bands (nm)
method1_band_centers <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)

# Gaussian weights for a center wavelength
.gaussian_weights <- function(center, wl, fwhm = 10) {
    sigma <- fwhm / (2 * sqrt(2 * log(2)))
    w <- exp(-0.5 * ((wl - center) / sigma)^2)
    w / sum(w)
}

# Convolve full spectrum to the Method 1 band set
convolve_to_bands <- function(spec, wl, centers = method1_band_centers, fwhm = 10) {
    if (is.null(dim(spec))) {
        spec <- matrix(spec, nrow = 1)
    }

    weights <- lapply(centers, .gaussian_weights, wl = wl, fwhm = fwhm)
    band_mat <- matrix(NA_real_, nrow = nrow(spec), ncol = length(centers))

    for (j in seq_along(centers)) {
        band_mat[, j] <- as.numeric(spec %*% weights[[j]])
    }

    colnames(band_mat) <- paste0("B", centers)
    if (nrow(spec) == 1) {
        setNames(as.numeric(band_mat[1, ]), colnames(band_mat))
    } else {
        band_mat
    }
}

# Sample parameter ranges (Option B: fix less sensitive terms)
.sample_params <- function(n) {
    tibble(
        N = runif(n, 1.0, 2.5),
        Cab = runif(n, 10, 80),
        Car = runif(n, 0, 20),
        Cbrown = 0, # fixed
        Cw = runif(n, 0.005, 0.02),
        Cm = 0.01, # fixed
        LAI = runif(n, 0.5, 8),
        psoil = runif(n, 0, 1),
        typelidf = 1
    )
}

# Generate PROSAIL LUT and downsample to method1 bands
build_prosail_lut <- function(n_samples = 2000, wl_full = 400:1100, fwhm = 20, seed = NULL) {
    params <- .sample_params(n_samples)
    spectra_list <- vector("list", n_samples)

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
            error = function(e) NULL
        )
    }

    valid <- !sapply(spectra_list, is.null)
    params <- params[valid, ]
    spectra_list <- spectra_list[valid]

    if (!length(spectra_list)) stop("PROSAIL LUT generation failed.")

    spectra_full <- do.call(rbind, spectra_list)
    all_wl <- 400:2500
    mask <- all_wl %in% wl_full
    spectra_trim <- spectra_full[, mask, drop = FALSE]

    bands <- convolve_to_bands(spectra_trim, wl_full, centers = method1_band_centers, fwhm = fwhm)

    list(
        bands = bands,
        params = params,
        wl_full = wl_full,
        spectra_trim = spectra_trim
    )
}

# Train one SVR per parameter
train_svr_models <- function(band_mat, params_df) {
    targets <- c("N", "Cab", "Car", "Cw", "LAI", "psoil")
    models <- list()
    band_df <- as.data.frame(band_mat)

    for (tgt in targets) {
        df <- cbind(band_df, y = params_df[[tgt]])
        models[[tgt]] <- svm(y ~ ., data = df, type = "eps-regression")
    }

    models
}

predict_params_svr <- function(models, band_vec) {
    band_df <- as.data.frame(t(band_vec))
    preds <- vapply(models, predict, numeric(1), newdata = band_df)
    as.list(preds)
}

forward_prosail_spectrum <- function(pars, wl_full = 400:1100) {
    spec_obj <- PROSAIL(
        N = pars$N,
        Cab = pars$Cab,
        Car = pars$Car,
        Cbrown = 0,
        Cw = pars$Cw,
        Cm = 0.01,
        LAI = pars$LAI,
        psoil = pars$psoil,
        TypeLidf = 1
    )

    vals <- if (inherits(spec_obj, "Speclib")) {
        m <- spectra(spec_obj)
        if (is.matrix(m)) m[1, ] else as.numeric(m)
    } else {
        as.numeric(spec_obj)
    }

    all_wl <- 400:2500
    mask <- all_wl %in% wl_full
    vals[mask]
}

# Convenience wrapper: build models from scratch
build_prosail_inversion <- function(n_lut = 2000, wl_full = 400:1100, fwhm = 20, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    lut <- build_prosail_lut(n_samples = n_lut, wl_full = wl_full, fwhm = fwhm, seed = seed)
    models <- train_svr_models(lut$bands, lut$params)
    list(models = models, wl_full = wl_full, fwhm = fwhm, n_lut = n_lut, seed = seed)
}

# Invert a single observation (band-level) then forward simulate full spectrum
invert_and_forward <- function(obs_band_vec, inversion, wl_full = 400:1100) {
    preds <- predict_params_svr(inversion$models, obs_band_vec)
    spec <- forward_prosail_spectrum(preds, wl_full = wl_full)
    list(params = preds, spectrum = spec)
}

# Save/Load helpers
save_prosail_inversion <- function(inversion, file_path) {
    dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(inversion, file = file_path)
    invisible(file_path)
}

load_prosail_inversion <- function(file_path) {
    if (!file.exists(file_path)) stop("Inversion model file not found: ", file_path)
    readRDS(file_path)
}

# Construct default filename for model saving
default_prosail_model_filename <- function(dir = "models", fwhm = 20, seed = 123, n_lut = 2000) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    file.path(dir, sprintf("prosail_inversion_fwhm-%d_seed-%d_nlut-%d.rds", fwhm, seed, n_lut))
}
