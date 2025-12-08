# Script to prepare the app for GitHub Pages deployment using shinylive

# 1. Install shinylive if needed
if (!requireNamespace("shinylive", quietly = TRUE)) {
    message("Installing shinylive...")
    install.packages("shinylive")
}

# 2. Create a temporary staging directory
# This ensures we only deploy the necessary files, not the whole project (like .git, output, etc.)
staging_dir <- "app_staging"
if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE)
dir.create(staging_dir)

# 3. Copy essential files to staging
file.copy("app.R", staging_dir)
file.copy("DESCRIPTION", staging_dir)

dir.create(file.path(staging_dir, "R"))
file.copy(list.files("R", full.names = TRUE), file.path(staging_dir, "R"))

dir.create(file.path(staging_dir, "data"))
file.copy(list.files("data", full.names = TRUE), file.path(staging_dir, "data"))

message("Files copied to staging area.")

# 4. Export to docs/
if (dir.exists("docs")) unlink("docs", recursive = TRUE)
shinylive::export(appdir = staging_dir, destdir = "docs")

# 5. Cleanup
unlink(staging_dir, recursive = TRUE)

message("Success! App exported to 'docs/'.")
