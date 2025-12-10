# Spectral Reconstruction

## Overview
This repository hosts a Shiny application for reconstructing hyperspectral signatures from multispectral samples using two complementary approaches:

1. **Method 1: Parametric Reconstruction**
   - Fits a piecewise spectral model that emulates the visible, red-edge, and NIR regions with gaussian, logistic, and quadratic baselines.
   - Designed for fast, analytic reconstruction and for visualizing parameter estimates (chlorophyll amplitude, red-edge inflection, NIR slope, etc.).

2. **Method 2: PROSAIL Inversion (hybrid LUT + SVM)**
   - Builds a PROSAIL LUT sampled over key biophysical parameters (N, Cab, Car, Cw, LAI, psoil) and converts those spectra to the multispectral anchor bands via 10 nm Gaussian SRFs.
   - Trains one Support Vector Regression (SVR) model per PROSAIL parameter using the convolved band responses.
   - At runtime the Shiny app loads a pretrained `.rds` inversion model, convolves the selected sample to the same multispectral bands, predicts the PROSAIL parameters, runs PROSAIL forward, and displays the 400–1100 nm spectral curve.

## Data and Pretrained Model
- `data/prosail_data.RData` contains the reference dataset used by Method 1 (observed spectra, wavelength axis, and parameter table).
- `models/prosail_inversion_fwhm-10_seed-123_nlut-1000.rds` is the pretrained PROSAIL inversion model (n_lut=1000, fwhm=10). The Shiny app auto-loads this file if it exists.

## Running the App
1. Ensure dependencies are installed (R packages `shiny`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `pracma`, `e1071`, `hsdar`).
2. From the repository root run:
   ```bash
   R -e "shiny::runApp('app')"
   ```
3. Navigate to `http://127.0.0.1:7395/` in your browser. The app presents two tabs:
   - **Method 1**: Choose a sample index to see the reference and reconstructed spectra along with the parameter table.
   - **Method 2**: The pretrained PROSAIL inversion model is auto-loaded. Select a sample index and click *Run PROSAIL inversion* to compare the observed hyperspectral curve against the PROSAIL (SVR) prediction, inspect estimated PROSAIL parameters, and view band-wise residuals.

## Building or Retraining the PROSAIL Inversion Model
If you want to retrain the inversion model with more samples or different seeds:
1. Install the same R dependencies plus `hsdar` for PROSAIL.
2. Run:
   ```bash
   Rscript scripts/build_prosail_inversion_model.R
   ```
   This will generate `models/prosail_inversion_fwhm-10_seed-123_nlut-1000.rds`.
3. Restart the Shiny app so it auto-loads the fresh model.

## Repository Structure
```
app/                       # Shiny UI/server scripts
R/                         # Supporting R functions (parametric and inversion helpers)
scripts/                   # Data/model generation utilities
models/                    # Pretrained PROSAIL inversion models
data/                      # Shared datasets (prosail_data.RData)
README.md                  # This file
```

## Testing
- The app can be launched locally with `shiny::runApp('app')`. Running Method 2 requires the pretrained `.rds`, which is provided in `models/`.
- For automated workflows, call `build_prosail_inversion()` from `R/prosail_inversion.R` to generate a model for offline training before deploying it to the app.
