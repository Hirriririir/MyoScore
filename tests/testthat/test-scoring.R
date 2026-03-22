test_that("myoscore_weights sums to ~1", {
  w <- myoscore_weights()
  expect_equal(sum(w), 1.001, tolerance = 0.01)
  expect_named(w, myoscore_dimensions())
})

test_that("myoscore_dimensions returns 5 dimensions", {
  dims <- myoscore_dimensions()
  expect_length(dims, 5)
  expect_true("Strength" %in% dims)
  expect_true("Youth" %in% dims)
})

test_that("myoscore_score_dimension returns NA for missing genes", {
  log2cpm <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(log2cpm) <- paste0("FakeGene", 1:5)
  colnames(log2cpm) <- paste0("S", 1:10)

  gene_weights <- data.frame(
    ID = paste0("RealGene", 1:3),
    weight = c(1, 2, 3),
    direction_v3 = c(1, -1, 1),
    dimension = "Strength",
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(
    myoscore_score_dimension(log2cpm, gene_weights, "Strength", verbose = FALSE)
  )
  expect_true(all(is.na(result)))
})

test_that("myoscore_score_dimension returns 0-100 range", {
  set.seed(42)
  n_genes <- 10
  n_samples <- 30

  log2cpm <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  gene_ids <- paste0("G", 1:n_genes)
  rownames(log2cpm) <- gene_ids
  colnames(log2cpm) <- paste0("S", 1:n_samples)

  gene_weights <- data.frame(
    ID = gene_ids,
    weight = runif(n_genes, 1, 10),
    direction_v3 = sample(c(-1, 1), n_genes, replace = TRUE),
    dimension = "Strength",
    stringsAsFactors = FALSE
  )

  result <- myoscore_score_dimension(log2cpm, gene_weights, "Strength",
                                     verbose = FALSE)
  expect_length(result, n_samples)
  expect_true(min(result) >= 0)
  expect_true(max(result) <= 100)
  expect_equal(min(result), 0)
  expect_equal(max(result), 100)
})

test_that("myoscore_score works with matrix input", {
  set.seed(123)
  n_genes <- 20
  n_samples <- 15

  # Create fake data using some real gene names from the package
  gene_weights <- data.frame(
    ID = paste0("TestGene", 1:n_genes),
    weight = runif(n_genes, 1, 10),
    direction_v3 = sample(c(-1, 1), n_genes, replace = TRUE),
    dimension = rep(myoscore_dimensions(), each = 4),
    stringsAsFactors = FALSE
  )

  counts <- matrix(rpois(n_genes * n_samples, lambda = 500),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- gene_weights$ID
  colnames(counts) <- paste0("S", 1:n_samples)

  result <- myoscore_score(counts, gene_weights = gene_weights, verbose = FALSE)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), n_samples)
  expect_true("MyoScore" %in% colnames(result))
  expect_true(all(c("Strength_score", "Mass_score", "LeanMuscle_score",
                     "Youth_score", "Resilience_score") %in% colnames(result)))
})

test_that("myoscore_score rejects single-sample input", {
  counts <- matrix(rpois(20, 100), nrow = 20, ncol = 1)
  rownames(counts) <- paste0("G", 1:20)
  colnames(counts) <- "S1"

  expect_error(myoscore_score(counts, verbose = FALSE), "At least 2 samples")
})

test_that(".validate_gene_weights catches missing columns", {
  bad_df <- data.frame(gene = "A", w = 1)
  expect_error(MyoScore:::.validate_gene_weights(bad_df), "missing columns")
})
