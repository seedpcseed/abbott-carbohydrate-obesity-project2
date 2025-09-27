# Use native R with user libraries
user_lib <- "/root/R/x86_64-pc-linux-gnu-library"
if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE)
}

.libPaths(c(user_lib, .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Welcome message
cat("Native R 4.5.1 loaded with user libraries\n")
