#' Cost Functions for Spectral Reconstruction
#'
#' This module defines cost functions to evaluate the fit of the reconstruction.
#' It supports both local (segment-specific) and global (full spectrum) costs.

#' Calculate Total Cost (Global + Local)
#'
#' @param params Named vector of parameters to optimize (e.g., c(C=0.5, sigma=40))
#' @param obs_anchors Named vector of anchor band reflectances
#' @param obs_validation Optional: Named vector of additional validation bands for global fitting
#' @param wl_grid Vector of wavelengths
#' @param weights List of weights for different cost components (default: equal)
#' @return Scalar cost value
calculate_cost <- function(params, obs_anchors, obs_validation = NULL, wl_grid, weights = list(global = 1, local = 1)) {
    # 1. Reconstruct using current parameters
    # Note: We need a version of reconstruct_spectrum that accepts explicit parameters
    # For now, we assume we can pass them or the function adapts.
    # Since reconstruct_spectrum currently derives most parameters,
    # we might need to refactor it to accept overrides.

    # Placeholder for parameter injection logic:
    # If we are optimizing 'C', we override the internal optimization.

    # For this implementation, let's assume we are evaluating the *result* of a reconstruction
    # against a set of validation bands.

    # ... (This requires refactoring reconstruction_model.R to accept params)

    return(0) # Placeholder
}

#' Local Cost: Red Edge Logistic Fit
#'
#' @param C Slope parameter of the logistic function
#' @param anchors Named vector containing R670, R700, R740, R780
#' @param lambda_i Inflection point (calculated externally or passed)
#' @return Sum of Squared Errors for the Red Edge anchors
cost_red_edge <- function(C, anchors, lambda_i) {
    # Predict at 700 and 740 nm
    # Formula: Rmin + (Rmax - Rmin) / (1 + exp(C * (lambda_i - lambda)))

    pred_700 <- anchors["R670"] + (anchors["R780"] - anchors["R670"]) / (1 + exp(C * (lambda_i - 700)))
    pred_740 <- anchors["R670"] + (anchors["R780"] - anchors["R670"]) / (1 + exp(C * (lambda_i - 740)))

    # RSS
    err_700 <- pred_700 - anchors["R700"]
    err_740 <- pred_740 - anchors["R740"]

    return(err_700^2 + err_740^2)
}

#' Global Cost: Full Spectrum Deviation
#'
#' @param reconstructed_spectrum Vector of reconstructed values
#' @param true_spectrum Vector of true values (must match length)
#' @return RMSE or SSE
cost_global <- function(reconstructed_spectrum, true_spectrum) {
    sum((reconstructed_spectrum - true_spectrum)^2, na.rm = TRUE)
}

#' Combined Cost Function for Optimization
#'
#' This function can be passed to optim().
#' It calculates the reconstruction with given parameters and compares to ALL available bands.
#'
#' @param x Vector of parameters being optimized (e.g. x[1] = C)
#' @param param_names Names of the parameters (e.g. "C")
#' @param reconstruction_fn Function that takes parameters and returns a spectrum
#' @param target_values Vector of observed values to fit against
#' @param target_indices Indices of the target values in the output spectrum
#' @return Cost value
cost_combined <- function(x, param_names, reconstruction_fn, target_values, target_indices) {
    # Assign names to x
    names(x) <- param_names

    # Generate spectrum
    # We assume reconstruction_fn takes a list of parameters
    spec <- reconstruction_fn(x)

    # Extract predicted values at target indices
    preds <- spec[target_indices]

    # Calculate Residuals
    residuals <- preds - target_values

    # Global SSE
    sse <- sum(residuals^2)

    return(sse)
}
