library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pracma)

# Load data and model
if (file.exists("data/prosail_data.RData")) {
    load("data/prosail_data.RData")
} else {
    # Fallback for local testing
    if (file.exists("../data/prosail_data.RData")) load("../data/prosail_data.RData")
}

# Source existing model for "Default" reconstruction
if (file.exists("R/reconstruction_model.R")) source("R/reconstruction_model.R")

# --- Extended Model for Optimization (Internal) ---
# This implements the "New" model with Blue and Water absorption terms
reconstruct_optimized <- function(obs_anchors, wl, params) {
    # Unpack params
    sigma_red <- params$sigma
    C <- params$C
    A_blue <- params$A_blue
    sigma_water <- params$sigma_water

    # Anchors
    r400 <- obs_anchors["R400"]
    r670 <- obs_anchors["R670"]
    r700 <- obs_anchors["R700"]
    r740 <- obs_anchors["R740"]
    r750 <- obs_anchors["R750"]
    r780 <- obs_anchors["R780"]
    r980 <- obs_anchors["R980"]
    r1100 <- obs_anchors["R1100"]

    # --- 1. Visible (400-680) ---
    # Baseline: Linear 400 -> 750
    slope_vis <- (r750 - r400) / (750 - 400)
    int_vis <- r400 - slope_vis * 400
    base_vis <- function(l) slope_vis * l + int_vis

    # Red Absorption (A derived)
    A_red <- base_vis(670) - r670

    # Blue Absorption (A_blue optimized, sigma_blue fixed)
    sigma_blue <- 20 # Fixed assumption

    vis_fun <- function(l) {
        base_vis(l) -
            A_red * exp(-(l - 670)^2 / (2 * sigma_red^2)) -
            A_blue * exp(-(l - 450)^2 / (2 * sigma_blue^2))
    }

    # --- 2. Red Edge (680-780) ---
    # Logistic
    rho_i <- (r670 + r780) / 2
    denom <- r740 - r700
    if (abs(denom) < 1e-6) denom <- 1e-6
    lambda_i <- 700 + 40 * ((rho_i - r700) / denom)

    re_fun <- function(l) {
        r670 + (r780 - r670) / (1 + exp(C * (lambda_i - l)))
    }

    # --- 3. NIR (780-1100) ---
    # Baseline: Linear 780 -> 1100
    slope_nir <- (r1100 - r780) / (1100 - 780)
    int_nir <- r780 - slope_nir * 780
    base_nir <- function(l) slope_nir * l + int_nir

    # Water Absorption at 980
    A_water <- base_nir(980) - r980

    nir_fun <- function(l) {
        base_nir(l) - A_water * exp(-(l - 980)^2 / (2 * sigma_water^2))
    }

    # Stitch
    rec <- numeric(length(wl))
    m_vis <- wl <= 680
    m_re <- wl > 680 & wl <= 780
    m_nir <- wl > 780

    rec[m_vis] <- vis_fun(wl[m_vis])
    rec[m_re] <- re_fun(wl[m_re])
    rec[m_nir] <- nir_fun(wl[m_nir])

    list(
        spectrum = rec,
        params = list(
            A_red = A_red, sigma_red = sigma_red,
            A_blue = A_blue, sigma_blue = sigma_blue,
            C = C, lambda_i = lambda_i,
            A_water = A_water, sigma_water = sigma_water
        )
    )
}

# Define UI
ui <- fluidPage(
    withMathJax(),
    titlePanel("Spectral Reconstruction Analysis"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("sample_idx", "Select Sample Index:",
                min = 1, max = nrow(spectra_mat), value = 1, step = 1
            ),
            hr(),
            h4("Optimized Parameters"),
            tableOutput("params_table"),
            hr(),
            helpText("Blue line shows the optimized reconstruction with additional Blue/Water absorption terms.")
        ),
        mainPanel(
            # Stacked Plots
            plotOutput("main_plot", height = "400px"),
            plotOutput("error_plot", height = "250px"),
            hr(),

            # Equations Section
            h3("Mathematical Model & Parameters"),
            fluidRow(
                column(
                    12,
                    h4("1. Visible Region (400 - 680 nm)"),
                    p("Modeled as a linear baseline minus two Gaussian absorption features (Chlorophyll at 670nm and Blue absorption at 450nm)."),
                    uiOutput("eq_vis"),
                    verbatimTextOutput("vals_vis")
                )
            ),
            fluidRow(
                column(
                    12,
                    h4("2. Red Edge (680 - 780 nm)"),
                    p("Modeled using a 4-parameter logistic function."),
                    uiOutput("eq_re"),
                    verbatimTextOutput("vals_re")
                )
            ),
            fluidRow(
                column(
                    12,
                    h4("3. NIR Plateau (780 - 1100 nm)"),
                    p("Modeled as a linear baseline minus a Gaussian water absorption feature at 980nm."),
                    uiOutput("eq_nir"),
                    verbatimTextOutput("vals_nir")
                )
            )
        )
    )
)

# Define Server
server <- function(input, output) {
    # Reactive: Perform Reconstruction & Optimization
    data_res <- reactive({
        idx <- input$sample_idx
        obs <- spectra_mat[idx, ]

        # 1. Default Reconstruction (Old Model)
        res_def <- reconstruct_spec_parametric(obs, wavelengths)

        # 2. Optimization (New Model)
        # Prepare anchors
        target_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
        idx_anchors <- vapply(target_bands, function(t) which.min(abs(wavelengths - t)), integer(1))
        obs_anchors <- obs[idx_anchors]
        names(obs_anchors) <- paste0("R", target_bands)

        # Cost Function (NO DATA-LEAKAGE)
        # Only compare reconstructed values at anchor wavelengths.
        # Sigma parameters are estimated from anchor points (not fitted using full spectrum).
        cost_fn <- function(p) {
            # p: [C, A_blue]
            C_val <- p[1]
            A_blue_val <- p[2]

            # --- derive sigma_red from anchor-based amplitudes ---
            slope_vis <- (obs_anchors["R750"] - obs_anchors["R400"]) / (750 - 400)
            int_vis <- obs_anchors["R400"] - slope_vis * 400
            base_vis <- function(l) slope_vis * l + int_vis
            A_red_est <- base_vis(670) - obs_anchors["R670"]
            A700_est <- base_vis(700) - obs_anchors["R700"]
            ratio <- A700_est / (A_red_est + 1e-12)
            if (!is.finite(ratio) || ratio <= 0 || ratio >= 1) {
                sigma_red_est <- 40
            } else {
                sigma_red_est <- sqrt((700 - 670)^2 / (-2 * log(ratio)))
            }

            # sigma_blue fixed (cannot be robustly estimated from sparse anchors)
            sigma_blue <- 20

            # --- derive sigma_water from anchor-based amplitudes ---
            slope_nir <- (obs_anchors["R1100"] - obs_anchors["R780"]) / (1100 - 780)
            int_nir <- obs_anchors["R780"] - slope_nir * 780
            base_nir <- function(l) slope_nir * l + int_nir
            A_water_est <- base_nir(980) - obs_anchors["R980"]
            A1100_est <- base_nir(1100) - obs_anchors["R1100"]
            ratio_w <- A1100_est / (A_water_est + 1e-12)
            if (!is.finite(ratio_w) || ratio_w <= 0 || ratio_w >= 1) {
                sigma_water_est <- 40
            } else {
                sigma_water_est <- sqrt((1100 - 980)^2 / (-2 * log(ratio_w)))
            }

            pars <- list(sigma = sigma_red_est, C = C_val, A_blue = A_blue_val, sigma_water = sigma_water_est)
            rec <- reconstruct_optimized(obs_anchors, wavelengths, pars)

            # Compare only at anchor indices -> prevents full-spectrum leakage
            rec_anchors <- rec$spectrum[idx_anchors]
            sum((rec_anchors - obs_anchors)^2)
        }

        # Initial Guesses (only for parameters that are optimized: C and A_blue)
        init <- c(0.5, 0.01)

        # Run Optim (optimize C and A_blue only)
        opt <- optim(init, cost_fn,
            method = "L-BFGS-B",
            lower = c(0.01, 0.0), upper = c(5.0, 0.2)
        )

        # Recompute sigma estimates with the best-fit parameters and produce final reconstruction
        best_C <- opt$par[1]
        best_A_blue <- opt$par[2]

        # recompute sigmas deterministically (same logic as in cost_fn)
        slope_vis <- (obs_anchors["R750"] - obs_anchors["R400"]) / (750 - 400)
        int_vis <- obs_anchors["R400"] - slope_vis * 400
        base_vis <- function(l) slope_vis * l + int_vis
        A_red_est <- base_vis(670) - obs_anchors["R670"]
        A700_est <- base_vis(700) - obs_anchors["R700"]
        ratio <- A700_est / (A_red_est + 1e-12)
        if (!is.finite(ratio) || ratio <= 0 || ratio >= 1) {
            sigma_red_est <- 40
        } else {
            sigma_red_est <- sqrt((700 - 670)^2 / (-2 * log(ratio)))
        }
        sigma_blue <- 20
        slope_nir <- (obs_anchors["R1100"] - obs_anchors["R780"]) / (1100 - 780)
        int_nir <- obs_anchors["R780"] - slope_nir * 780
        base_nir <- function(l) slope_nir * l + int_nir
        A_water_est <- base_nir(980) - obs_anchors["R980"]
        A1100_est <- base_nir(1100) - obs_anchors["R1100"]
        ratio_w <- A1100_est / (A_water_est + 1e-12)
        if (!is.finite(ratio_w) || ratio_w <= 0 || ratio_w >= 1) {
            sigma_water_est <- 40
        } else {
            sigma_water_est <- sqrt((1100 - 980)^2 / (-2 * log(ratio_w)))
        }

        best_pars <- list(sigma = sigma_red_est, C = best_C, A_blue = best_A_blue, sigma_water = sigma_water_est)
        res_opt <- reconstruct_optimized(obs_anchors, wavelengths, best_pars)

        list(
            obs = obs,
            def = res_def$spectrum,
            opt = res_opt$spectrum,
            params = res_opt$params
        )
    })

    # Main Plot
    output$main_plot <- renderPlot({
        d <- data_res()

        df <- data.frame(
            WL = wavelengths,
            Reflectance = as.numeric(d$obs),
            Type = "Original"
        )
        df_def <- data.frame(
            WL = wavelengths,
            Reflectance = d$def,
            Type = "Default"
        )
        df_opt <- data.frame(
            WL = wavelengths,
            Reflectance = d$opt,
            Type = "Optimized"
        )

        plot_data <- rbind(df, df_def, df_opt)

        ggplot(plot_data, aes(x = WL, y = Reflectance, color = Type, linetype = Type)) +
            geom_line(size = 1) +
            scale_color_manual(values = c("Original" = "black", "Default" = "red", "Optimized" = "blue")) +
            scale_linetype_manual(values = c("Original" = "solid", "Default" = "dashed", "Optimized" = "dashed")) +
            theme_minimal() +
            theme(
                legend.position = "top",
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                plot.margin = margin(b = 0),
                text = element_text(size = 14)
            ) +
            labs(title = paste("Sample", input$sample_idx, "Reconstruction"), y = "Reflectance")
    })

    # Error Plot
    output$error_plot <- renderPlot({
        d <- data_res()

        err_def <- d$def - as.numeric(d$obs)
        err_opt <- d$opt - as.numeric(d$obs)

        df_err <- rbind(
            data.frame(WL = wavelengths, Error = err_def, Type = "Default"),
            data.frame(WL = wavelengths, Error = err_opt, Type = "Optimized")
        )

        ggplot(df_err, aes(x = WL, y = Error, color = Type)) +
            geom_line(size = 1) +
            geom_hline(yintercept = 0, color = "gray") +
            scale_color_manual(values = c("Default" = "red", "Optimized" = "blue")) +
            theme_minimal() +
            theme(
                legend.position = "none",
                plot.margin = margin(t = 0),
                text = element_text(size = 14)
            ) +
            labs(y = "Error (Rec - Obs)", x = "Wavelength (nm)")
    })

    # Equations & Values
    output$eq_vis <- renderUI({
        withMathJax("$$ R_{vis}(\\lambda) = R_{base}(\\lambda) - A \\cdot \\exp\\left(-\\frac{(\\lambda - 670)^2}{2\\sigma^2}\\right) - A_{blue} \\cdot \\exp\\left(-\\frac{(\\lambda - 450)^2}{2\\sigma_{blue}^2}\\right) $$")
    })
    output$vals_vis <- renderText({
        p <- data_res()$params
        sprintf(
            "Parameters: A = %.4f, sigma = %.2f, A_blue = %.4f, sigma_blue = %.2f",
            p$A_red, p$sigma_red, p$A_blue, p$sigma_blue
        )
    })

    output$eq_re <- renderUI({
        withMathJax("$$ R_{re}(\\lambda) = R_{min} + \\frac{R_{max} - R_{min}}{1 + \\exp(C(\\lambda_i - \\lambda))} $$")
    })
    output$vals_re <- renderText({
        p <- data_res()$params
        sprintf("Parameters: C = %.4f, lambda_i = %.2f", p$C, p$lambda_i)
    })

    output$eq_nir <- renderUI({
        withMathJax("$$ R_{nir}(\\lambda) = R_{nir\\_base}(\\lambda) - A_{water} \\cdot \\exp\\left(-\\frac{(\\lambda - 980)^2}{2\\sigma_{water}^2}\\right) $$")
    })
    output$vals_nir <- renderText({
        p <- data_res()$params
        sprintf("Parameters: A_water = %.4f, sigma_water = %.2f", p$A_water, p$sigma_water)
    })

    # Sidebar Table
    output$params_table <- renderTable(
        {
            p <- data_res()$params
            data.frame(
                Param = c("C", "Sigma (Red)", "A (Blue)", "Sigma (Water)"),
                Value = sprintf("%.4f", c(p$C, p$sigma_red, p$A_blue, p$sigma_water))
            )
        },
        colnames = FALSE
    )
}

# Run the application
shinyApp(ui = ui, server = server)
