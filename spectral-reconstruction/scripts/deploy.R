# Deploy script for Shinylive
# This script prepares the app in 'app_staging' and exports it to 'docs'.

if (!require("shinylive", quietly = TRUE)) {
    install.packages("shinylive", repos = "https://cran.rstudio.com/")
}
if (!require("fs", quietly = TRUE)) {
    install.packages("fs", repos = "https://cran.rstudio.com/")
}

# Define paths
root_dir <- "."
staging_dir <- "app_staging"
dest_dir <- "docs"

# Create/Clean staging directory
if (fs::dir_exists(staging_dir)) fs::dir_delete(staging_dir)
fs::dir_create(staging_dir)

# Copy App Files
message("Copying files to staging area...")
fs::file_copy("app.R", file.path(staging_dir, "app.R"))
fs::dir_copy("R", file.path(staging_dir, "R"))
fs::dir_copy("data", file.path(staging_dir, "data"))

# Export the app
message("Exporting Shiny app to ", dest_dir, "...")
shinylive::export(staging_dir, dest_dir)

message("Export complete!")
message("To deploy to GitHub Pages:")
message("1. Commit the 'docs' directory.")
message("2. Push to GitHub.")
message("3. Go to GitHub Repository Settings -> Pages.")
message("4. Select 'main' branch and '/docs' folder.")
