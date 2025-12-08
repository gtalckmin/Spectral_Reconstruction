library(shiny)
library(ggplot2)
library(tidyverse)
library(pracma)

# Load data and model
load("data/prosail_data.RData")
# source("R/reconstruction_model.R") # Shiny auto-sources files in R/ directory

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
            h4("Reconstructed Parameters"),
            tableOutput("params_table"),
            hr(),
            helpText("Adjust the slider to view different samples from the PROSAIL dataset.")
        ),
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "Visualisation",
                    plotOutput("main_plot", height = "500px"),
                    h4("Model Performance"),
                    verbatimTextOutput("metrics")
                ),
                tabPanel(
                    "Error Analysis",
                    plotOutput("error_plot", height = "500px")
                ),
                tabPanel(
                    "Equations",
                    h3("Mathematical Model"),
                    p("The spectrum is reconstructed using a piecewise approach with three segments:"),
                    h4("1. Visible Region (400 - 680 nm)"),
                    p("Modeled as a linear baseline minus two Gaussian absorption features (Chlorophyll at 670nm and Blue absorption at 450nm)."),
                    uiOutput("eq_vis"),
                    h4("2. Red Edge (680 - 780 nm)"),
                    p("Modeled using a 4-parameter logistic function."),
                    uiOutput("eq_re"),
                    h4("3. NIR Plateau (780 - 1100 nm)"),
                    p("Modeled as a linear baseline minus a Gaussian water absorption feature at 980nm."),
                    uiOutput("eq_nir")
                )
            )
        )
    )
)

# Define Server
server <- function(input, output) {
    # Reactive reconstruction
    reconstruction <- reactive({
        idx <- input$sample_idx
        obs_spectrum <- spectra_mat[idx, ]

        # Create named vector for the function
        names(obs_spectrum) <- paste0("R", wavelengths)

        # Reconstruct
        res <- reconstruct_spectrum(obs_spectrum, wavelengths)

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
                title = paste("Sample", input$sample_idx, "- Reference vs Reconstruction"),
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
                    "A (Chl)", "Sigma (Chl)", "A (Blue)", "Sigma (Blue)",
                    "Lambda_i", "C (Slope)", "A (Water)", "Sigma (Water)", "AUC (NIR)"
                ),
                Value = sprintf("%.4f", c(
                    feats$A, feats$sigma, feats$A_blue, feats$sigma_blue,
                    feats$lambda_i, feats$C, feats$A_water, feats$sigma_water, feats$AUC
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
}

# Run the application
shinyApp(ui = ui, server = server)
