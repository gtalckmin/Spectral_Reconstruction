# Build a demo PROSAIL inversion model and save as an .rds
# Adjust n_lut for runtime/quality tradeoff

options(width = 200)

required <- c("hsdar", "e1071")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing)) {
    stop(sprintf("Missing required packages: %s. Install them before running this script.", paste(missing, collapse = ", ")))
}

# Source inversion helpers
script_paths <- c("R/prosail_inversion.R", "app/R/prosail_inversion.R")
src <- script_paths[file.exists(script_paths)][1]
if (is.na(src) || !file.exists(src)) stop("Could not find R/prosail_inversion.R in repo. Ensure file exists.")
source(src)

# Build LUT and train SVR models
n_lut <- 1000L
wl_full <- 400:1100
fwhm <- 10
seed <- 123

cat(sprintf("Building PROSAIL inversion model with n_lut=%d, fwhm=%d, seed=%d\n", n_lut, fwhm, seed))
inv <- build_prosail_inversion(n_lut = n_lut, wl_full = wl_full, fwhm = fwhm, seed = seed)

# Ensure models directory exists
dir.create("models", showWarnings = FALSE)
out_file <- file.path("models", sprintf("prosail_inversion_fwhm-%d_seed-%d_nlut-%d.rds", fwhm, seed, n_lut))
cat(sprintf("Saving model to %s\n", out_file))
saveRDS(inv, out_file)
cat("Done.\n")
