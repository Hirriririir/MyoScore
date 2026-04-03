#' Radar Chart of MyoScore Dimensions
#'
#' @description
#' Plot a radar (spider) chart showing the five MyoScore dimensions.
#' Supports plotting one or more groups (e.g., disease stages) as
#' overlaid or faceted panels.
#'
#' @param scores A data.frame from [myoscore_score()], or a named numeric
#'   vector of length 5 (one value per dimension), or a matrix/data.frame
#'   where each row is a group and columns are dimension scores.
#' @param groups Optional. A factor or character vector assigning each row
#'   of \code{scores} to a group. When provided, group means are plotted.
#' @param colors Optional. Character vector of colors (one per group).
#'   Default uses [myoscore_colors()] with type `"spectrum"`.
#' @param facet Logical. If \code{TRUE} and multiple groups exist, plot
#'   each group in a separate panel. Default \code{TRUE}.
#' @param title Optional main title.
#' @param show_values Logical. Show score values at vertices. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{fmsb::radarchart()}.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#'
#' @details
#' Requires the \pkg{fmsb} package (in Suggests).
#'
#' @export
#' @examples
#' # Radar chart from a named vector of dimension scores
#' dim_scores <- c(Strength = 55, Mass = 48, LeanMuscle = 42,
#'                 Youth = 60, Resilience = 50)
#' \donttest{
#' myoscore_plot_radar(dim_scores)
#' }
myoscore_plot_radar <- function(scores,
                                groups = NULL,
                                colors = NULL,
                                facet = TRUE,
                                title = NULL,
                                show_values = TRUE,
                                ...) {
  .check_suggests("fmsb", "for radar charts")

  dim_cols <- paste0(myoscore_dimensions(), "_score")
  dim_labels <- myoscore_dimensions()

  # Parse input into a matrix of group means
  if (is.numeric(scores) && is.null(dim(scores))) {
    # Named vector of 5 values — single sample radar with dimension colors
    if (length(scores) != 5) stop("Numeric vector must have 5 elements.")
    means_mat <- matrix(scores, nrow = 1)
    colnames(means_mat) <- dim_labels
    rownames(means_mat) <- if (!is.null(title)) title else "Sample"
    group_names <- rownames(means_mat)
    group_n <- NA
    .single_vector <- TRUE
  } else {
    .single_vector <- FALSE
    # Data.frame from myoscore_score()
    scores_df <- as.data.frame(scores)

    # Check for dimension score columns
    if (all(dim_cols %in% colnames(scores_df))) {
      mat <- as.matrix(scores_df[, dim_cols])
    } else if (all(dim_labels %in% colnames(scores_df))) {
      mat <- as.matrix(scores_df[, dim_labels])
    } else {
      stop("scores must contain dimension score columns.")
    }
    colnames(mat) <- dim_labels

    if (is.null(groups)) {
      # Overall mean
      means_mat <- matrix(colMeans(mat, na.rm = TRUE), nrow = 1)
      colnames(means_mat) <- dim_labels
      rownames(means_mat) <- "Overall"
      group_names <- "Overall"
      group_n <- nrow(mat)
    } else {
      groups <- as.character(groups)
      if (length(groups) != nrow(mat)) {
        stop("groups must have the same length as nrow(scores).")
      }
      group_names <- unique(groups)
      means_mat <- t(sapply(group_names, function(g) {
        colMeans(mat[groups == g, , drop = FALSE], na.rm = TRUE)
      }))
      colnames(means_mat) <- dim_labels
      group_n <- sapply(group_names, function(g) sum(groups == g))
    }
  }

  n_groups <- nrow(means_mat)

  # Colors
  if (is.null(colors)) {
    if (n_groups <= 5) {
      palette <- myoscore_colors("spectrum")
      colors <- unname(palette[seq_len(n_groups)])
    } else {
      colors <- grDevices::rainbow(n_groups)
    }
  }
  fill_colors <- grDevices::adjustcolor(colors, alpha.f = 0.25)

  # fmsb format: row1=max, row2=min, then data
  radar_data <- rbind(
    rep(100, 5),
    rep(0, 5),
    means_mat
  )
  colnames(radar_data) <- dim_labels
  radar_data <- as.data.frame(radar_data)

  if (facet && n_groups > 1) {
    # Faceted: one panel per group
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1, n_groups), mar = c(1, 1, 3, 1))

    for (i in seq_len(n_groups)) {
      single <- rbind(radar_data[1, ], radar_data[2, ], radar_data[i + 2, ])
      fmsb::radarchart(single,
                       axistype = 1,
                       pcol = colors[i],
                       pfcol = fill_colors[i],
                       plwd = 3, plty = 1,
                       cglcol = "grey80", cglty = 1, cglwd = 0.8,
                       axislabcol = "grey50",
                       caxislabels = c("0", "25", "50", "75", "100"),
                       calcex = 0.6, vlcex = 1.0,
                       vlabels = dim_labels,
                       ...)

      if (show_values) {
        .add_radar_values(means_mat[i, ], colors[i])
      }

      label <- group_names[i]
      if (!is.na(group_n[i])) {
        label <- sprintf("%s (n=%d)", label, group_n[i])
      }
      graphics::title(main = label, cex.main = 0.95,
                      col.main = colors[i], font.main = 2)
    }

    if (!is.null(title)) {
      graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2, line = -1)
    }
  } else if (.single_vector && n_groups == 1) {
    # Single-sample radar with per-dimension colored segments
    dim_colors <- unname(myoscore_colors("dimensions"))

    fmsb::radarchart(radar_data,
                     axistype = 1,
                     pcol = "grey40", pfcol = grDevices::adjustcolor("grey80", 0.15),
                     plwd = 2, plty = 1,
                     cglcol = "grey80", cglty = 1, cglwd = 0.8,
                     axislabcol = "grey50",
                     caxislabels = c("0", "25", "50", "75", "100"),
                     calcex = 0.7, vlcex = 1.0,
                     vlabels = dim_labels,
                     ...)

    # Draw colored segments per dimension
    angles <- seq(pi / 2, pi / 2 - 2 * pi * (4 / 5), length.out = 5)
    vals <- means_mat[1, ] / 100
    for (k in seq_along(vals)) {
      k2 <- if (k < 5) k + 1 else 1
      x_pts <- c(0, vals[k] * cos(angles[k]), vals[k2] * cos(angles[k2]))
      y_pts <- c(0, vals[k] * sin(angles[k]), vals[k2] * sin(angles[k2]))
      graphics::polygon(x_pts, y_pts,
                        col = grDevices::adjustcolor(dim_colors[k], 0.35),
                        border = dim_colors[k], lwd = 2)
    }

    if (show_values) {
      for (k in seq_along(vals)) {
        r <- vals[k]
        x <- r * cos(angles[k])
        y <- r * sin(angles[k])
        graphics::text(x * 1.15, y * 1.15, sprintf("%.1f", means_mat[1, k]),
                       cex = 0.85, font = 2, col = dim_colors[k])
      }
    }

    graphics::legend("bottomright", legend = dim_labels,
                     fill = dim_colors, cex = 0.75, bty = "n")

    if (!is.null(title)) {
      graphics::title(main = title, cex.main = 1.2, font.main = 2)
    }
  } else {
    # Overlaid: all groups in one chart
    fmsb::radarchart(radar_data,
                     axistype = 1,
                     pcol = colors, pfcol = fill_colors,
                     plwd = 2, plty = 1,
                     cglcol = "grey80", cglty = 1, cglwd = 0.8,
                     axislabcol = "grey50",
                     caxislabels = c("0", "25", "50", "75", "100"),
                     calcex = 0.7, vlcex = 1.0,
                     vlabels = dim_labels,
                     ...)

    if (n_groups > 1) {
      graphics::legend("topright", legend = group_names,
                       col = colors, lwd = 2, cex = 0.8, bty = "n")
    }

    if (!is.null(title)) {
      graphics::title(main = title, cex.main = 1.2, font.main = 2)
    }
  }

  invisible(NULL)
}


#' Add score values at radar chart vertices
#' @noRd
.add_radar_values <- function(values, color) {
  angles <- seq(pi / 2, pi / 2 - 2 * pi * (4 / 5), length.out = 5)
  for (k in seq_along(values)) {
    r <- values[k] / 100
    x <- r * cos(angles[k])
    y <- r * sin(angles[k])
    graphics::text(x * 1.15, y * 1.15, sprintf("%.1f", values[k]),
                   cex = 0.75, font = 2, col = color)
  }
}
