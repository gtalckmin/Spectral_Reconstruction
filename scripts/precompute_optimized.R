#!/usr/bin/env Rscript
# Precompute optimized reconstructions (no data leakage) for all samples
library(jsonlite)

# Load PROSAIL data
data_paths <- c('app/data/prosail_data.RData', 'data/prosail_data.RData', '../app/data/prosail_data.RData')
found <- FALSE
for (p in data_paths) {
  if (file.exists(p)) {
    load(p)
    found <- TRUE
    break
  }
}
if (!found) stop('prosail_data.RData not found in expected locations')

# Expect objects: spectra_mat (n x m) and wavelengths (m)
if (!exists('spectra_mat') || !exists('wavelengths')) stop('spectra_mat or wavelengths missing in RData')

target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
idx_anchors_fun <- function(wl) vapply(target_bands, function(t) which.min(abs(wl - t)), integer(1))

reconstruct_optimized_local <- function(obs_anchors, wl, pars) {
  sigma_red <- pars$sigma
  C <- pars$C
  A_blue <- pars$A_blue
  sigma_water <- pars$sigma_water

  r400 <- obs_anchors['R400']
  r670 <- obs_anchors['R670']
  r700 <- obs_anchors['R700']
  r740 <- obs_anchors['R740']
  r750 <- obs_anchors['R750']
  r780 <- obs_anchors['R780']
  r980 <- obs_anchors['R980']
  r1100 <- obs_anchors['R1100']

  slope_vis <- (r750 - r400) / (750 - 400)
  int_vis <- r400 - slope_vis * 400
  base_vis <- function(l) slope_vis * l + int_vis
  A_red <- base_vis(670) - r670
  sigma_blue <- 20
  vis_fun <- function(l) {
    base_vis(l) - A_red * exp(-(l - 670)^2 / (2 * sigma_red^2)) - A_blue * exp(-(l - 450)^2 / (2 * sigma_blue^2))
  }

  rho_i <- (r670 + r780) / 2
  denom <- r740 - r700
  if (abs(denom) < 1e-6) denom <- 1e-6
  lambda_i <- 700 + 40 * ((rho_i - r700) / denom)
  re_fun <- function(l) r670 + (r780 - r670) / (1 + exp(C * (lambda_i - l)))

  slope_nir <- (r1100 - r780) / (1100 - 780)
  int_nir <- r780 - slope_nir * 780
  base_nir <- function(l) slope_nir * l + int_nir
  A_water <- base_nir(980) - r980
  nir_fun <- function(l) base_nir(l) - A_water * exp(-(l - 980)^2 / (2 * sigma_water^2))

  rec <- numeric(length(wl))
  m_vis <- wl <= 680
  m_re <- wl > 680 & wl <= 780
  m_nir <- wl > 780
  rec[m_vis] <- vis_fun(wl[m_vis])
  rec[m_re] <- re_fun(wl[m_re])
  rec[m_nir] <- nir_fun(wl[m_nir])

  list(spectrum = rec, params = list(A_red = A_red, sigma_red = sigma_red, A_blue = A_blue, sigma_blue = sigma_blue, C = C, lambda_i = lambda_i, A_water = A_water, sigma_water = sigma_water))
}

optimize_sample <- function(obs, wl) {
  idx_anchors <- idx_anchors_fun(wl)
  obs_anchors <- obs[idx_anchors]
  names(obs_anchors) <- paste0('R', target_bands)

  cost_fn <- function(p) {
    C_val <- p[1]
    A_blue_val <- p[2]

    slope_vis <- (obs_anchors['R750'] - obs_anchors['R400']) / (750 - 400)
    int_vis <- obs_anchors['R400'] - slope_vis * 400
    base_vis <- function(l) slope_vis * l + int_vis
    A_red_est <- base_vis(670) - obs_anchors['R670']
    A700_est <- base_vis(700) - obs_anchors['R700']
    ratio <- A700_est / (A_red_est + 1e-12)
    if (!is.finite(ratio) || ratio <= 0 || ratio >= 1) sigma_red_est <- 40 else sigma_red_est <- sqrt((700 - 670)^2 / (-2 * log(ratio)))

    slope_nir <- (obs_anchors['R1100'] - obs_anchors['R780']) / (1100 - 780)
    int_nir <- obs_anchors['R780'] - slope_nir * 780
    base_nir <- function(l) slope_nir * l + int_nir
    A_water_est <- base_nir(980) - obs_anchors['R980']
    A1100_est <- base_nir(1100) - obs_anchors['R1100']
    ratio_w <- A1100_est / (A_water_est + 1e-12)
    if (!is.finite(ratio_w) || ratio_w <= 0 || ratio_w >= 1) sigma_water_est <- 40 else sigma_water_est <- sqrt((1100 - 980)^2 / (-2 * log(ratio_w)))

    pars <- list(sigma = sigma_red_est, C = C_val, A_blue = A_blue_val, sigma_water = sigma_water_est)
    rec <- reconstruct_optimized_local(obs_anchors, wl, pars)
    rec_anchors <- rec$spectrum[idx_anchors]
    sum((rec_anchors - obs_anchors)^2)
  }

  init <- c(0.5, 0.01)
  opt <- optim(init, cost_fn, method = 'L-BFGS-B', lower = c(0.01, 0), upper = c(5, 0.2))
  best_C <- opt$par[1]
  best_A_blue <- opt$par[2]

  # recompute sigmas deterministically
  slope_vis <- (obs_anchors['R750'] - obs_anchors['R400']) / (750 - 400)
  int_vis <- obs_anchors['R400'] - slope_vis * 400
  base_vis <- function(l) slope_vis * l + int_vis
  A_red_est <- base_vis(670) - obs_anchors['R670']
  A700_est <- base_vis(700) - obs_anchors['R700']
  ratio <- A700_est / (A_red_est + 1e-12)
  if (!is.finite(ratio) || ratio <= 0 || ratio >= 1) sigma_red_est <- 40 else sigma_red_est <- sqrt((700 - 670)^2 / (-2 * log(ratio)))
  slope_nir <- (obs_anchors['R1100'] - obs_anchors['R780']) / (1100 - 780)
  int_nir <- obs_anchors['R780'] - slope_nir * 780
  base_nir <- function(l) slope_nir * l + int_nir
  A_water_est <- base_nir(980) - obs_anchors['R980']
  A1100_est <- base_nir(1100) - obs_anchors['R1100']
  ratio_w <- A1100_est / (A_water_est + 1e-12)
  if (!is.finite(ratio_w) || ratio_w <= 0 || ratio_w >= 1) sigma_water_est <- 40 else sigma_water_est <- sqrt((1100 - 980)^2 / (-2 * log(ratio_w)))

  final_pars <- list(sigma = sigma_red_est, C = best_C, A_blue = best_A_blue, sigma_water = sigma_water_est)
  final_rec <- reconstruct_optimized_local(obs_anchors, wl, final_pars)
  list(spectrum = final_rec$spectrum, params = final_rec$params)
}

# Iterate samples
n_samples <- nrow(spectra_mat)
results_spectra <- vector('list', n_samples)
results_params <- vector('list', n_samples)
for (i in seq_len(n_samples)) {
  obs <- spectra_mat[i, ]
  res <- optimize_sample(obs, wavelengths)
  results_spectra[[i]] <- res$spectrum
  results_params[[i]] <- res$params
  if (i %% 50 == 0) message('Processed sample ', i, ' / ', n_samples)
}

out <- list(wavelengths = wavelengths, spectra = results_spectra, params = results_params)

dir.create('docs/data', recursive = TRUE, showWarnings = FALSE)
dir.create('app/data', recursive = TRUE, showWarnings = FALSE)
json_path_docs <- file.path('docs', 'data', 'optimized_reconstructions.json')
json_path_app <- file.path('app', 'data', 'optimized_reconstructions.json')
write(toJSON(out, digits = 8, auto_unbox = TRUE), json_path_docs)
write(toJSON(out, digits = 8, auto_unbox = TRUE), json_path_app)
message('Wrote optimized reconstructions to: ', json_path_docs)
