library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pracma)
library(e1071)
library(hsdar)

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("wavelengths", "spectra_mat", "Wavelength", "Reflectance", "Type", "Error"))
}

# Load data and model
if (file.exists("data/prosail_data.RData")) {
    load("data/prosail_data.RData")
} else {
    stop("Required data file 'data/prosail_data.RData' not found. Please ensure it is included in the deployment.")
}
source("R/reconstruction_model.R")
if (file.exists("R/prosail_inversion.R")) {
    source("R/prosail_inversion.R")
} else if (file.exists("../R/prosail_inversion.R")) {
    source("../R/prosail_inversion.R")
} else {
    # If the file isn't present locally, we error early with a helpful message
    stop("Required 'prosail_inversion.R' not found in app/R or ../R. Please add the file or run the project from the repository root.")
}

# Define UI
ui <- fluidPage(
    withMathJax(),
    titlePanel("Spectral Reconstruction Analysis"),
    tabsetPanel(
        tabPanel(
            "Method 1: Parametric",
            fluidRow(
                column(
                    width = 3,
                    sliderInput("sample_idx_m1", "Select Sample Index:",
                        min = 1, max = nrow(spectra_mat), value = 1, step = 1
                    ),
                    hr(),
                    h4("Reconstructed Parameters"),
                    tableOutput("params_table"),
                    hr(),
                    tags$details(
                        tags$summary("Show equations"),
                        h4("1. Visible Region (400 - 680 nm)"),
                        p("Linear baseline minus chlorophyll/blue Gaussian absorption."),
                        uiOutput("eq_vis"),
                        h4("2. Red Edge (680 - 780 nm)"),
                        p("4-parameter logistic."),
                        uiOutput("eq_re"),
                        h4("3. NIR Plateau (780 - 1100 nm)"),
                        p("Quadratic baseline with water absorption."),
                        uiOutput("eq_nir")
                    ),
                    helpText("All outputs shown on a single page for Method 1.")
                ),
                column(
                    width = 9,
                    plotOutput("main_plot", height = "320px"),
                    verbatimTextOutput("metrics"),
                    plotOutput("error_plot", height = "220px")
                )
            )
        ),
        tabPanel(
            "Method 2: PROSAIL Inversion",
            fluidRow(
                column(
                    width = 3,
                    sliderInput("sample_idx_m2", "Select Sample Index:",
                        min = 1, max = nrow(spectra_mat), value = 1, step = 1
                    ),
                    # Gaussian FWHM is fixed to 10 nm for the app (multispectral SRF)
                    # numeric input is intentionally removed; kept for compatibility and read-only display
                    tags$div(tags$strong("FWHM (Gaussian)", ": 10 nm (fixed)")),
                    verbatimTextOutput("m2_model_seed"),
                    actionButton("run_m2", "Run PROSAIL inversion", class = "btn-primary"),
                    fileInput("upload_m2_model", "Upload trained model (.rds):", accept = c(".rds")),
                    br(),
                    hr(),
                    helpText("SVM(Multispectral data) = PROSAIL parameters -> PROSAIL Forward = spectral signature."),
                    hr(),
                    h4("Estimated PROSAIL Parameters (SVM outputs)"),
                    tableOutput("m2_params_table"),
                    hr(),
                    h4("Observed vs PROSAIL Predicted Band Reflectance"),
                    tableOutput("m2_band_comp_table"),
                    hr(),
                    textOutput("m2_status")
                ),
                column(
                    width = 9,
                    plotOutput("m2_plot", height = "320px"),
                    verbatimTextOutput("m2_metrics"),
                    plotOutput("m2_error_plot", height = "220px")
                )
            )
        )
    )
)

# Define Server
server <- function(input, output, session) {
    # check for required packages at runtime and show a friendly message
    required_pkgs <- c("hsdar", "e1071", "tibble")
    missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs)) {
        stop(paste(
            "Missing required packages:", paste(missing_pkgs, collapse = ", "),
            "
Please install them with: install.packages(c('hsdar','e1071','tibble')) or via devtools if needed."
        ))
    }
    # --- Method 1: Parametric reconstruction ---
    reconstruction <- reactive({
        idx <- input$sample_idx_m1
        obs_spectrum <- spectra_mat[idx, ]

        names(obs_spectrum) <- paste0("R", wavelengths)
        res <- reconstruct_spec_parametric(obs_spectrum, wavelengths)

        list(
            obs = obs_spectrum,
            rec = res$spectrum,
            features = res$features
        )
    })

    # Main Plot
    output$main_plot <- renderPlot({
        data <- reconstruction()

        df <- data.frame(
            Wavelength = wavelengths,
            Reflectance = as.numeric(data$obs),
            Type = "Reference"
        )

        df_rec <- data.frame(
            Wavelength = wavelengths,
            Reflectance = data$rec,
            Type = "Reconstructed"
        )

        plot_df <- rbind(df, df_rec)

        ggplot(plot_df, aes(x = Wavelength, y = Reflectance, color = Type)) +
            geom_line(size = 1) +
            scale_color_manual(values = c("Reference" = "black", "Reconstructed" = "red")) +
            theme_minimal() +
            labs(
                title = paste("Sample", input$sample_idx_m1, "- Reference vs Reconstruction"),
                y = "Reflectance", x = "Wavelength (nm)"
            ) +
            theme(legend.position = "top", text = element_text(size = 14))
    })

    # Error Plot
    output$error_plot <- renderPlot({
        data <- reconstruction()

        error <- data$rec - as.numeric(data$obs)

        df_err <- data.frame(
            Wavelength = wavelengths,
            Error = error
        )

        ggplot(df_err, aes(x = Wavelength, y = Error)) +
            geom_line(color = "blue", size = 1) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
            theme_minimal() +
            labs(
                title = "Reconstruction Error (Reconstructed - Reference)",
                y = "Error", x = "Wavelength (nm)"
            ) +
            theme(text = element_text(size = 14))
    })

    # Parameters Table
    output$params_table <- renderTable(
        {
            feats <- reconstruction()$features

            data.frame(
                Parameter = c(
                    "A (Chl)", "Sigma (Chl)",
                    "Lambda_i", "C (Slope)", "AUC (NIR)"
                ),
                Value = sprintf("%.4f", c(
                    feats$A, feats$sigma,
                    feats$lambda_i, feats$C, feats$AUC
                ))
            )
        },
        colnames = FALSE
    )

    # Metrics
    output$metrics <- renderText({
        data <- reconstruction()
        rmse <- sqrt(mean((data$rec - as.numeric(data$obs))^2))
        paste("RMSE:", round(rmse, 5))
    })

    # Equations
    output$eq_vis <- renderUI({
        withMathJax(
            "$$ R_{vis}(\\lambda) = R_{base}(\\lambda) - A \\cdot \\exp\\left(-\\frac{(\\lambda - 670)^2}{2\\sigma^2}\\right) - A_{blue} \\cdot \\exp\\left(-\\frac{(\\lambda - 450)^2}{2\\sigma_{blue}^2}\\right) $$"
        )
    })

    output$eq_re <- renderUI({
        withMathJax(
            "$$ R_{re}(\\lambda) = R_{min} + \\frac{R_{max} - R_{min}}{1 + \\exp(C(\\lambda_i - \\lambda))} $$"
        )
    })

    output$eq_nir <- renderUI({
        withMathJax(
            "$$ R_{nir}(\\lambda) = R_{nir\\_base}(\\lambda) - A_{water} \\cdot \\exp\\left(-\\frac{(\\lambda - 980)^2}{2\\sigma_{water}^2}\\right) $$"
        )
    })

    # --- Method 2: PROSAIL inversion (SVR + forward) ---
    inversion_models <- reactiveVal(NULL)
    inversion_status <- reactiveVal("No PROSAIL model loaded. Upload a trained .rds model.")

    # Load uploaded trained inversion model (.rds) — training is offline
    observeEvent(input$upload_m2_model, {
        req(input$upload_m2_model)
        infile <- input$upload_m2_model$datapath
        tryCatch(
            {
                inv <- load_prosail_inversion(infile)
                inversion_models(inv)
                inversion_status(paste("Uploaded model loaded:", input$upload_m2_model$name))
                output$m2_model_seed <- renderText({
                    if (!is.null(inv$seed)) as.character(inv$seed) else "NA"
                })
            },
            error = function(e) {
                inversion_status(paste("Failed to load uploaded model:", e$message))
            }
        )
    })

    output$m2_status <- renderText({
        inversion_status()
    })

    # Try to auto-load demo model from `models/` if present
    default_model_path <- "models/prosail_inversion_fwhm-10_seed-123_nlut-1000.rds"
    if (file.exists(default_model_path)) {
        tryCatch(
            {
                inv0 <- load_prosail_inversion(default_model_path)
                inversion_models(inv0)
                inversion_status(paste("Auto-loaded model:", basename(default_model_path)))
                output$m2_model_seed <- renderText({
                    if (!is.null(inv0$seed)) as.character(inv0$seed) else "NA"
                })
            },
            error = function(e) {
                inversion_status(paste("Failed to auto-load default model:", e$message))
                output$m2_model_seed <- renderText({
                    "No model loaded"
                })
            }
        )
    } else {
        output$m2_model_seed <- renderText({
            "No model loaded"
        })
    }

    # (removed) manual load from saved models; upload is the supported flow

    method2_results <- eventReactive(input$run_m2, {
        inv <- inversion_models()
        req(!is.null(inv))

        idx <- input$sample_idx_m2
        obs_spectrum <- spectra_mat[idx, ]
        band_obs <- convolve_to_bands(obs_spectrum, wavelengths, fwhm = 10)

        res <- invert_and_forward(band_obs, inv, wl_full = wavelengths)

        list(
            obs = obs_spectrum,
            rec = res$spectrum,
            params = res$params,
            bands = band_obs
        )
    })

    # Observed vs Predicted band reflectance comparison table
    output$m2_band_comp_table <- renderTable(
        {
            data <- method2_results()
            req(!is.null(data))
            band_obs <- data$bands
            inv <- inversion_models()
            pred_rec <- data$rec
            pred_bands <- convolve_to_bands(pred_rec, wavelengths, fwhm = 10)

            if (is.null(names(band_obs))) {
                centers <- method1_band_centers
            } else {
                centers <- as.numeric(gsub("^B", "", names(band_obs)))
            }

            obs_vals <- as.numeric(band_obs)
            pred_vals <- as.numeric(pred_bands)
            df <- data.frame(
                Band = centers,
                Observed = sprintf("%.5f", obs_vals),
                Predicted = sprintf("%.5f", pred_vals),
                Residual = sprintf("%.5f", obs_vals - pred_vals)
            )
            df
        },
        colnames = TRUE
    )


    output$m2_plot <- renderPlot({
        data <- method2_results()
        req(!is.null(data))

        df <- data.frame(
            Wavelength = wavelengths,
            Reflectance = as.numeric(data$obs),
            Type = "Reference"
        )

        df_rec <- data.frame(
            Wavelength = wavelengths,
            Reflectance = data$rec,
            Type = "PROSAIL (SVR)"
        )

        plot_df <- rbind(df, df_rec)

        p <- ggplot(plot_df, aes(x = Wavelength, y = Reflectance, color = Type)) +
            geom_line(size = 1) +
            scale_color_manual(values = c("Reference" = "black", "PROSAIL (SVR)" = "darkgreen")) +
            theme_minimal() +
            labs(
                title = paste("Sample", input$sample_idx_m2, "- Reference vs PROSAIL"),
                y = "Reflectance", x = "Wavelength (nm)"
            ) +
            theme(legend.position = "top", text = element_text(size = 14))

        # Add observed/predicted band points and labels
        inv <- inversion_models()
        if (!is.null(inv)) {
            band_obs <- data$bands
            pred_rec <- data$rec
            pred_bands <- convolve_to_bands(pred_rec, wavelengths, fwhm = 10)
            centers <- if (is.null(names(band_obs))) method1_band_centers else as.numeric(gsub("^B", "", names(band_obs)))
            obs_vals <- as.numeric(band_obs)
            pred_vals <- as.numeric(pred_bands)
            bands_df_obs <- data.frame(Wavelength = centers, Reflectance = obs_vals)
            bands_df_pred <- data.frame(Wavelength = centers, Reflectance = pred_vals)
            p <- p +
                geom_point(data = bands_df_obs, aes(x = Wavelength, y = Reflectance), color = "blue", size = 3, shape = 16) +
                geom_point(data = bands_df_pred, aes(x = Wavelength, y = Reflectance), color = "red", size = 3, shape = 17) +
                geom_text(data = bands_df_obs, aes(x = Wavelength, y = Reflectance, label = sprintf("%.3f", Reflectance)), vjust = -1.2, color = "blue", size = 3) +
                geom_text(data = bands_df_pred, aes(x = Wavelength, y = Reflectance, label = sprintf("%.3f", Reflectance)), vjust = 1.6, color = "red", size = 3)
        }

        p
    })

    output$m2_error_plot <- renderPlot({
        data <- method2_results()
        req(!is.null(data))

        error <- data$rec - as.numeric(data$obs)

        df_err <- data.frame(
            Wavelength = wavelengths,
            Error = error
        )

        ggplot(df_err, aes(x = Wavelength, y = Error)) +
            geom_line(color = "blue", size = 1) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
            theme_minimal() +
            labs(
                title = "PROSAIL Reconstruction Error",
                y = "Error", x = "Wavelength (nm)"
            ) +
            theme(text = element_text(size = 14))
    })

    output$m2_params_table <- renderTable(
        {
            data <- method2_results()
            req(!is.null(data))

            params <- data$params
            # Include fixed parameters used for forward PROSAIL
            params$Cbrown <- 0
            params$Cm <- 0.01

            data.frame(
                Parameter = c("N", "Cab", "Car", "Cbrown", "Cw", "Cm", "LAI", "psoil"),
                Value = sprintf("%.4f", c(params$N, params$Cab, params$Car, params$Cbrown, params$Cw, params$Cm, params$LAI, params$psoil))
            )
        },
        colnames = FALSE
    )

    output$m2_metrics <- renderText({
        data <- method2_results()
        req(!is.null(data))
        rmse <- sqrt(mean((data$rec - as.numeric(data$obs))^2))
        paste("RMSE:", round(rmse, 5))
    })
}

# Run the application
shinyApp(ui = ui, server = server)
