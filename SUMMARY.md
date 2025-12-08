# Project Summary

## Achievements
1.  **Restructured Project**: Organized into `R/`, `data/`, `output/`.
2.  **Modularized Logic**: `R/reconstruction_model.R` now contains the core logic, separated from the UI.
3.  **Global Cost Component**: Implemented `R/cost_functions.R` and updated `reconstruct_spectrum` to accept external parameters (`params`). This allows optimizing the model against the full spectrum.
4.  **Model Enhancements**:
    *   **Visible**: Added a "Blue Absorption" term (Gaussian at 450nm) to better model the 400-600nm range.
    *   **NIR**: Switched from a simple Quadratic to a **Linear Baseline + Water Absorption** model (Gaussian at 980nm) to better reflect physical properties.
5.  **Error Analysis**: Created `R/analyze_errors_ranges.R` to quantify performance across spectral regions.

## Next Steps
1.  **Debug Optimization**: The `R/optimize_parameters.R` script is set up to find the best global parameters but requires debugging of data types.
2.  **Refine Constraints**: Ensure the new Blue and Water absorption terms are constrained to physically realistic values during optimization.
3.  **Deploy**: Update `app.R` to expose these new parameters (e.g., sliders for `A_blue` or `sigma_water`) for interactive exploration.
