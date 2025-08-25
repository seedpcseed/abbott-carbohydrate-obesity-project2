# Enhanced autocomplete for R in VS Code/Cursor
if (interactive()) {
  # Enable better autocomplete
  options(completion.install = TRUE)
  
  # Set up better tab completion
  utils::rc.settings(ipck = TRUE)
  
  # Enable fuzzy matching for autocomplete
  options(completion.install = TRUE)
  
  # Set up better history
  if (!exists(".First")) {
    .First <- function() {
      if (interactive()) {
        cat("\nWelcome to R!\n")
        cat("Enhanced autocomplete is enabled.\n")
      }
    }
  }
}
