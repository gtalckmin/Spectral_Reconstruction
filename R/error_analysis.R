#' Calculate Reconstruction Error per Band
#' @param true_spectra Matrix of true spectra (samples x bands)
#' @param rec_spectra Matrix of reconstructed spectra (samples x bands)
#' @param wavelengths Vector of wavelengths corresponding to columns
#' @return Data frame with RMSE per wavelength
calculate_band_error <- function(true_spectra, rec_spectra, wavelengths) {
    # Calculate squared errors
    sq_diff <- (true_spectra - rec_spectra)^2

    # Mean Squared Error per band
    mse <- colMeans(sq_diff, na.rm = TRUE)
    rmse <- sqrt(mse)

    return(data.frame(
        Wavelength = wavelengths,
        RMSE = rmse
    ))
}

#' Find Max Error Regions
#' @param error_df Output from calculate_band_error
#' @param top_n Number of peak error regions to return
get_max_error_regions <- function(error_df, top_n = 3) {
    # Simple peak finding or sorting
    error_df[order(-error_df$RMSE), ][1:top_n, ]
}
