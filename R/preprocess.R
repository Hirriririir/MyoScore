#' Preprocess Raw Counts to log2(CPM+1)
#'
#' @description
#' Normalize raw RNA-seq count data using CPM (Counts Per Million) followed
#' by log2 transformation. This is the standard preprocessing step before
#' MyoScore calculation.
#'
#' @param raw_counts A numeric matrix or data.frame of raw counts with genes
#'   as rows and samples as columns. Row names should be gene symbols.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A numeric matrix of log2(CPM+1) values with the same dimensions
#'   and names as the input.
#'
#' @details
#' The transformation pipeline is:
#' \enumerate{
#'   \item CPM: counts / total_counts * 1e6
#'   \item log2(CPM + 1)
#' }
#'
#' @export
#' @examples
#' # Create example count matrix
#' counts <- matrix(rpois(500, lambda = 100), nrow = 50, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:50)
#' colnames(counts) <- paste0("Sample", 1:10)
#'
#' log2cpm <- myoscore_preprocess(counts)
myoscore_preprocess <- function(raw_counts, verbose = TRUE) {
  # Coerce to numeric matrix
  mat <- as.matrix(raw_counts)
  storage.mode(mat) <- "double"

  if (any(mat < 0, na.rm = TRUE)) {
    stop("raw_counts contains negative values. Expected raw count data.")
  }

  if (verbose) {
    message("Preprocessing: Raw Count -> CPM -> log2(CPM+1)")
    message(sprintf("  Input: %d genes x %d samples", nrow(mat), ncol(mat)))
  }

  # CPM normalization
  total_counts <- colSums(mat, na.rm = TRUE)
  if (any(total_counts == 0)) {
    warning("Some samples have zero total counts and will produce NaN values.")
  }
  cpm <- sweep(mat, 2, total_counts, "/") * 1e6

  # log2 transform
  log2cpm <- log2(cpm + 1)

  if (verbose) {
    message(sprintf("  log2(CPM+1) range: %.2f - %.2f",
                    min(log2cpm, na.rm = TRUE), max(log2cpm, na.rm = TRUE)))
  }

  log2cpm
}
