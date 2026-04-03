#' Boxplot of MyoScore by Groups
#'
#' @description
#' Create grouped boxplots comparing MyoScore or individual dimension scores
#' across conditions. Uses base R graphics by default, or ggplot2 if available.
#'
#' @param scores A data.frame from [myoscore_score()].
#' @param groups A factor or character vector of group labels (one per sample).
#' @param which Which score to plot. One of \code{"MyoScore"},
#'   \code{"Strength"}, \code{"Mass"}, \code{"LeanMuscle"}, \code{"Youth"},
#'   \code{"Resilience"}, or \code{"all"} for a multi-panel figure.
#'   Default \code{"MyoScore"}.
#' @param colors Optional named or positional color vector.
#' @param use_ggplot Logical. Use ggplot2 if available. Default \code{TRUE}.
#' @param title Optional main title.
#' @param ... Additional arguments passed to \code{boxplot()} or
#'   \code{ggplot2::geom_boxplot()}.
#'
#' @return If \code{use_ggplot = TRUE} and ggplot2 is available, returns a
#'   ggplot object. Otherwise, invisible \code{NULL}.
#'
#' @export
#' @examples
#' # Create example scores and groups
#' scores_df <- data.frame(
#'   Strength_score  = c(rnorm(5, 55, 5), rnorm(5, 40, 5)),
#'   Mass_score      = c(rnorm(5, 50, 5), rnorm(5, 45, 5)),
#'   LeanMuscle_score = c(rnorm(5, 48, 5), rnorm(5, 38, 5)),
#'   Youth_score     = c(rnorm(5, 52, 5), rnorm(5, 35, 5)),
#'   Resilience_score = c(rnorm(5, 50, 5), rnorm(5, 45, 5)),
#'   MyoScore        = c(rnorm(5, 50, 3), rnorm(5, 40, 3))
#' )
#' groups <- rep(c("Healthy", "Disease"), each = 5)
#' \donttest{
#' myoscore_plot_boxplot(scores_df, groups = groups)
#' }
myoscore_plot_boxplot <- function(scores,
                                  groups,
                                  which = "MyoScore",
                                  colors = NULL,
                                  use_ggplot = TRUE,
                                  title = NULL,
                                  ...) {
  scores_df <- as.data.frame(scores)
  groups <- as.factor(groups)

  if (length(groups) != nrow(scores_df)) {
    stop("groups must have the same length as nrow(scores).")
  }

  # Resolve column names
  dims <- myoscore_dimensions()
  col_map <- stats::setNames(
    c(paste0(dims, "_score"), "MyoScore"),
    c(dims, "MyoScore")
  )

  if (which == "all") {
    plot_cols <- col_map
  } else {
    if (!which %in% names(col_map)) {
      stop(sprintf("which must be one of: %s",
                   paste(names(col_map), collapse = ", ")))
    }
    plot_cols <- col_map[which]
  }

  # Default colors
  if (is.null(colors)) {
    n_groups <- nlevels(groups)
    spec <- myoscore_colors("spectrum")
    if (n_groups <= length(spec)) {
      colors <- unname(spec[seq_len(n_groups)])
    } else {
      colors <- grDevices::rainbow(n_groups)
    }
  }

  # ggplot2 path
  if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    return(.boxplot_ggplot(scores_df, groups, plot_cols, colors, title, ...))
  }

  # Base R fallback
  .boxplot_base(scores_df, groups, plot_cols, colors, title, ...)
}


#' ggplot2 boxplot implementation
#' @noRd
.boxplot_ggplot <- function(scores_df, groups, plot_cols, colors, title, ...) {
  # Build long-form data
  long_list <- list()
  for (nm in names(plot_cols)) {
    col <- plot_cols[nm]
    if (!col %in% colnames(scores_df)) next
    long_list[[nm]] <- data.frame(
      Group = groups,
      Score = scores_df[[col]],
      Dimension = nm,
      stringsAsFactors = FALSE
    )
  }

  if (length(long_list) == 0) {
    stop("No matching score columns found in scores data.frame. ",
         "Expected columns: ", paste(plot_cols, collapse = ", "))
  }

  long <- do.call(rbind, long_list)
  long$Dimension <- factor(long$Dimension, levels = names(plot_cols))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = Group, y = Score,
                                           fill = Group)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.7, ...) +
    ggplot2::geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(y = "Score (0-100)", x = NULL)

  if (length(plot_cols) > 1) {
    p <- p + ggplot2::facet_wrap(~Dimension, scales = "free_y")
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}


#' Base R boxplot fallback
#' @noRd
.boxplot_base <- function(scores_df, groups, plot_cols, colors, title, ...) {
  n_panels <- length(plot_cols)

  if (n_panels > 1) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    nc <- min(n_panels, 3)
    nr <- ceiling(n_panels / nc)
    graphics::par(mfrow = c(nr, nc), mar = c(5, 4, 3, 1))
  }

  for (nm in names(plot_cols)) {
    col <- plot_cols[nm]
    if (!col %in% colnames(scores_df)) next

    graphics::boxplot(scores_df[[col]] ~ groups,
                      col = colors, main = nm,
                      ylab = "Score (0-100)", xlab = "",
                      las = 2, ...)
    graphics::stripchart(scores_df[[col]] ~ groups,
                         add = TRUE, vertical = TRUE, method = "jitter",
                         pch = 16, cex = 0.4, col = "grey30")
  }

  if (!is.null(title) && n_panels == 1) {
    graphics::title(main = title)
  }

  invisible(NULL)
}
