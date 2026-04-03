#' Calculate MyoScore from Raw Count Data
#'
#' @description
#' Main entry point for computing MyoScore (Genetic Muscle Health Score).
#' Accepts either a file path or a count matrix, and returns per-sample
#' scores for all five dimensions plus the composite score.
#'
#' @param input Either a file path (character) to a raw count CSV/TSV, or a
#'   numeric matrix/data.frame with genes as rows and samples as columns.
#'   Gene symbols must be row names.
#' @param gene_weights Optional. A data.frame of gene weights with columns
#'   \code{ID}, \code{weight}, \code{direction_v3}, \code{dimension}.
#'   Default uses the built-in \code{myoscore_genes} dataset.
#' @param sep Separator for reading CSV files. Default \code{","}.
#'   Use \code{"\\t"} for tab-separated files.
#' @param min_coverage Minimum fraction (0-1) of genes required per dimension.
#'   Dimensions below this threshold return \code{NA}. Default \code{0.1}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A data.frame with samples as rows and columns:
#'   \code{Strength_score}, \code{Mass_score}, \code{LeanMuscle_score},
#'   \code{Youth_score}, \code{Resilience_score}, \code{MyoScore}.
#'
#' @details
#' ## Scoring Pipeline
#' \enumerate{
#'   \item Raw counts are normalized to log2(CPM+1).
#'   \item For each dimension, available genes are z-score standardized
#'         (gene-wise across all input samples).
#'   \item Z-scores are multiplied by gene direction and weight, then
#'         averaged (weighted mean).
#'   \item Raw dimension scores are min-max normalized to 0-100.
#'   \item Composite MyoScore is a weighted sum of the five dimensions.
#' }
#'
#' ## Interpretation
#' Higher scores indicate healthier muscle. The composite MyoScore
#' ranges from 0 (severe myopathy) to 100 (optimal muscle health).
#'
#' ## Important Notes
#' - Requires >= 20 samples for meaningful min-max normalization.
#' - Single-sample scoring is not recommended (use a reference cohort).
#' - Typical bulk RNA-seq datasets contain ~417 of the 1,116 scoring genes.
#'
#' @export
#' @examples
#' # Create a small example count matrix (50 genes x 10 samples)
#' set.seed(42)
#' genes <- head(MyoScore::myoscore_genes$ID, 50)
#' counts <- matrix(rpois(50 * 10, lambda = 100), nrow = 50,
#'                  dimnames = list(genes, paste0("S", 1:10)))
#' \donttest{
#' scores <- myoscore_score(counts, verbose = FALSE)
#' head(scores)
#' }
myoscore_score <- function(input,
                           gene_weights = NULL,
                           sep = ",",
                           min_coverage = 0.1,
                           verbose = TRUE) {

  if (verbose) {
    message("============================================================")
    message("MyoScore Calculator (v3.3)")
    message("============================================================")
  }

  # Load input data
  if (is.character(input)) {
    if (!file.exists(input)) {
      stop(sprintf("Input file not found: %s", input))
    }
    if (verbose) message(sprintf("\nLoading: %s", input))
    raw_counts <- utils::read.csv(input, row.names = 1,
                                  check.names = FALSE, sep = sep)
  } else if (is.matrix(input) || is.data.frame(input)) {
    raw_counts <- input
  } else {
    stop("input must be a file path, matrix, or data.frame.")
  }

  if (verbose) {
    message(sprintf("  Shape: %d genes x %d samples",
                    nrow(raw_counts), ncol(raw_counts)))
  }

  if (ncol(raw_counts) < 2) {
    stop("At least 2 samples are required for scoring.")
  }

  # Load gene weights
  if (is.null(gene_weights)) {
    gene_weights <- .load_builtin_genes(verbose)
  }
  .validate_gene_weights(gene_weights)

  # Preprocess
  log2cpm <- myoscore_preprocess(raw_counts, verbose = verbose)

  # Calculate dimension scores
  dims <- myoscore_dimensions()
  if (verbose) message("\nCalculating dimension scores...")

  scores <- data.frame(row.names = colnames(log2cpm))
  for (dim in dims) {
    scores[[paste0(dim, "_score")]] <- myoscore_score_dimension(
      log2cpm, gene_weights, dim,
      min_coverage = min_coverage, verbose = verbose
    )
  }

  # Composite score
  weights <- myoscore_weights()
  if (verbose) message("\nCalculating composite MyoScore...")

  composite <- rep(0, nrow(scores))
  weight_sum <- 0
  for (dim in dims) {
    col <- paste0(dim, "_score")
    if (!all(is.na(scores[[col]]))) {
      composite <- composite + weights[[dim]] * scores[[col]]
      weight_sum <- weight_sum + weights[[dim]]
    }
  }

  # Rescale if some dimensions are NA
  if (weight_sum > 0 && weight_sum < 1) {
    composite <- composite / weight_sum
  }
  scores$MyoScore <- composite

  if (verbose) {
    message(sprintf("\nMyoScore range: %.2f - %.2f",
                    min(scores$MyoScore, na.rm = TRUE),
                    max(scores$MyoScore, na.rm = TRUE)))
    message(sprintf("MyoScore mean: %.2f, sd: %.2f",
                    mean(scores$MyoScore, na.rm = TRUE),
                    stats::sd(scores$MyoScore, na.rm = TRUE)))
    message("\n============================================================")
    message(sprintf("Samples processed: %d", nrow(scores)))
    message("Done!")
  }

  scores
}


#' Calculate Score for a Single Dimension
#'
#' @param log2cpm Numeric matrix of log2(CPM+1) values (genes x samples).
#' @param gene_weights Data.frame with columns \code{ID}, \code{weight},
#'   \code{direction_v3}, \code{dimension}.
#' @param dimension Character. One of the five MyoScore dimensions.
#' @param min_coverage Minimum gene coverage fraction. Default \code{0.1}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return Numeric vector of dimension scores (0-100), one per sample.
#' @export
myoscore_score_dimension <- function(log2cpm,
                                     gene_weights = NULL,
                                     dimension,
                                     min_coverage = 0.1,
                                     verbose = TRUE) {
  dims <- myoscore_dimensions()
  if (!dimension %in% dims) {
    stop(sprintf("dimension must be one of: %s", paste(dims, collapse = ", ")))
  }

  if (is.null(gene_weights)) {
    gene_weights <- .load_builtin_genes(verbose = FALSE)
  }

  # Subset to this dimension

  dim_genes <- gene_weights[gene_weights$dimension == dimension, ]
  n_total <- nrow(dim_genes)

  # Find available genes
  available <- dim_genes[dim_genes$ID %in% rownames(log2cpm), ]
  n_available <- nrow(available)

  if (n_available == 0) {
    warning(sprintf("No genes found for %s dimension.", dimension))
    return(rep(NA_real_, ncol(log2cpm)))
  }

  coverage <- n_available / n_total
  if (verbose) {
    message(sprintf("  %s: %d/%d genes (%.1f%%)",
                    dimension, n_available, n_total, coverage * 100))
  }

  if (coverage < min_coverage) {
    warning(sprintf("%s coverage (%.1f%%) below threshold (%.1f%%). Returning NA.",
                    dimension, coverage * 100, min_coverage * 100))
    return(rep(NA_real_, ncol(log2cpm)))
  }

  # Extract expression for available genes
  expr <- log2cpm[available$ID, , drop = FALSE]

  # Gene-wise z-score standardization
  expr_mean <- rowMeans(expr)
  expr_sd <- apply(expr, 1, stats::sd)
  expr_sd[expr_sd == 0] <- 1  # avoid division by zero

  z_scores <- (expr - expr_mean) / expr_sd

  # Apply direction and weight
  directions <- available$direction_v3
  gene_w <- available$weight

  weighted_z <- sweep(z_scores, 1, directions * gene_w, "*")
  dim_score_raw <- colSums(weighted_z) / sum(gene_w)

  # Min-max normalize to 0-100
  min_val <- min(dim_score_raw)
  max_val <- max(dim_score_raw)

  if (max_val == min_val) {
    return(rep(50, ncol(log2cpm)))
  }

  (dim_score_raw - min_val) / (max_val - min_val) * 100
}


# ---- Internal helpers ----

#' Load built-in gene weights dataset
#' @noRd
.load_builtin_genes <- function(verbose = TRUE) {
  # Try package data first
  if (exists("myoscore_genes", where = asNamespace("MyoScore"),
             inherits = FALSE)) {
    genes <- get("myoscore_genes", envir = asNamespace("MyoScore"))
  } else {
    # Fallback to inst/extdata
    gene_file <- system.file("extdata", "myoscore_genes.csv",
                             package = "MyoScore")
    if (gene_file == "") {
      stop("Gene weights file not found. Provide gene_weights manually.")
    }
    genes <- utils::read.csv(gene_file, stringsAsFactors = FALSE)
  }

  if (verbose) {
    message(sprintf("Loaded %d gene weights (%d unique genes)",
                    nrow(genes), length(unique(genes$ID))))
  }
  genes
}

#' Validate gene weights data.frame
#' @noRd
.validate_gene_weights <- function(gene_weights) {
  required <- c("ID", "weight", "direction_v3", "dimension")
  missing <- setdiff(required, colnames(gene_weights))
  if (length(missing) > 0) {
    stop(sprintf("gene_weights missing columns: %s",
                 paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}
