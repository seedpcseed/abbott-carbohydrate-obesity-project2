# Force R to use only conda libraries to avoid symbol conflicts
# This prevents mixing system R packages with conda R packages

# Set library paths to only use conda R
.libPaths(c("/root/miniconda/lib/R/library"))

# Set R_LIBS environment variable
Sys.setenv(R_LIBS = "/root/miniconda/lib/R/library")

# Optional: Set R_LIBS_USER to avoid conflicts
Sys.setenv(R_LIBS_USER = "/root/miniconda/lib/R/library")

cat("R profile loaded: Using conda R libraries only\n")
cat("Library paths:", paste(.libPaths(), collapse = ", "), "\n")