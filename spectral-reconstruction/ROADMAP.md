# Spectral Reconstruction Roadmap

## Goal
Reconstruct hyperspectral data (400-1100nm) from multispectral observations using piecewise parametric models.

## Phases

### Phase 1: Architecture & Cleanup (Current)
- [x] Restructure project into `R/`, `data/`, and `output/`.
- [x] Modularize reconstruction logic into `R/reconstruction_model.R`.
- [x] Implement error analysis tools.

### Phase 2: Advanced Fitting & Optimization
- [ ] Implement `R/cost_functions.R` for global and range-specific error calculation.
- [ ] Replace analytical solutions with numerical optimization (`optim`) where beneficial.
- [ ] Investigate "Piecewise Polynomial" (Spline) alternatives for specific regions.

### Phase 3: Band Selection & Sensitivity
- [ ] Use error analysis to identify spectral regions with highest reconstruction error.
- [ ] Test adding new bands to the input set to reduce error.

### Phase 4: Visualization & Deployment
- [ ] Update Shiny App to use the modular `R/` functions.
- [ ] Deploy or package the solution.
