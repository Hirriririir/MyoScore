#' MyoScore Dimension Weights and Constants
#'
#' @description
#' Data-driven weights for the five MyoScore dimensions, derived from
#' GWAS-TWAS integration of 28 muscle-related phenotypes.
#'
#' @details
#' Weights represent the relative contribution of each dimension to overall
#' muscle health, determined by variance explained in the 1,722-sample
#' training cohort.
#'
#' @name myoscore-constants
NULL

#' Get MyoScore dimension weights
#'
#' @return Named numeric vector of dimension weights (sum to 1.0).
#' @export
#' @examples
#' myoscore_weights()
myoscore_weights <- function() {
  c(
    Strength    = 0.252,
    Mass        = 0.177,
    LeanMuscle  = 0.243,
    Youth       = 0.242,
    Resilience  = 0.087
  )
}

#' Get MyoScore dimension names
#'
#' @return Character vector of the five dimension names.
#' @export
#' @examples
#' myoscore_dimensions()
myoscore_dimensions <- function() {
  c("Strength", "Mass", "LeanMuscle", "Youth", "Resilience")
}

#' Get MyoScore color palette
#'
#' @param type One of \code{"dimensions"} (5 dimension colors),
#'   \code{"spectrum"} (unhealthy-to-healthy gradient), or \code{"all"}.
#' @return Named character vector of hex color codes.
#' @export
#' @examples
#' myoscore_colors("dimensions")
#' myoscore_colors("spectrum")
myoscore_colors <- function(type = c("dimensions", "spectrum", "all")) {
  type <- match.arg(type)

  dims <- c(
    Strength   = "#50327b",
    Mass       = "#46508b",
    LeanMuscle = "#f4e030",
    Youth      = "#72c95e",
    Resilience = "#31848f"
  )

  spectrum <- c(
    unhealthy = "#50327b",
    mild      = "#46508b",
    neutral   = "#31848f",
    positive  = "#72c95e",
    healthy   = "#f4e030"
  )

  switch(type,
    dimensions = dims,
    spectrum   = spectrum,
    all        = c(dims, spectrum)
  )
}
