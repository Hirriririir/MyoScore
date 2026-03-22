#' MyoScore Gene Weights
#'
#' A data frame containing 591 gene-dimension entries (417 unique genes)
#' used in MyoScore calculation, filtered to genes detectable in bulk
#' RNA-seq datasets.
#'
#' @format A data frame with 591 rows and 4 columns:
#' \describe{
#'   \item{ID}{Gene symbol (HGNC).}
#'   \item{weight}{Gene weight derived from TWAS Z-scores
#'     (|mean_Z| / n_phenotypes).}
#'   \item{direction_v3}{Direction of effect: +1 means high expression
#'     indicates health; -1 means high expression indicates disease.}
#'   \item{dimension}{One of five dimensions: Strength, Mass, LeanMuscle,
#'     Youth, Resilience.}
#' }
#'
#' @details
#' Genes were identified through TWAS (Transcriptome-Wide Association Study)
#' using FUSION with GTEx v8 skeletal muscle eQTL weights and 28 GWAS
#' phenotypes covering grip strength, body composition, MRI fat infiltration,
#' telomere length, and myopathy diagnoses.
#'
#' @source Myopathy Spectrum Project, GWAS-TWAS integration pipeline.
#' @examples
#' data(myoscore_genes)
#' table(myoscore_genes$dimension)
"myoscore_genes"
