# Script to prepare the app for GitHub Pages deployment using shinylive

# 1. Install required packages
# Note: On Linux, you may need to install 'libarchive-dev' first:
# sudo apt-get install libarchive-dev

if (!requireNamespace("shinylive", quietly = TRUE)) {
  install.packages("shinylive")
}
if (!requireNamespace("httpuv", quietly = TRUE)) {
  install.packages("httpuv")
}

# 2. Export the app to the 'docs' directory
# This compiles the app to WebAssembly-ready static files
shinylive::export(appdir = ".", destdir = "docs")

# 3. Test locally (optional)
# httpuv::runStaticServer("docs")

message("App exported to 'docs/' directory.")
message("Now push the changes to GitHub and enable GitHub Pages for the /docs folder.")
