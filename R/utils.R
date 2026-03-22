# Suppress R CMD check NOTEs for ggplot2 aes variables
utils::globalVariables(c("Group", "Score", "Dimension"))

#' Check if a suggested package is available
#' @param pkg Package name.
#' @param reason Why the package is needed.
#' @noRd
.check_suggests <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(reason)) msg <- paste0(msg, " ", reason)
    stop(paste0(msg, ". Install with: install.packages('", pkg, "')"),
         call. = FALSE)
  }
}
