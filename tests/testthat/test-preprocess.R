test_that("myoscore_preprocess returns correct dimensions", {
  counts <- matrix(rpois(200, lambda = 100), nrow = 20, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:20)
  colnames(counts) <- paste0("Sample", 1:10)

  result <- myoscore_preprocess(counts, verbose = FALSE)

  expect_equal(dim(result), dim(counts))
  expect_equal(rownames(result), rownames(counts))
  expect_equal(colnames(result), colnames(counts))
})

test_that("myoscore_preprocess produces log2(CPM+1) values", {
  counts <- matrix(c(100, 200, 300, 400), nrow = 2, ncol = 2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")

  result <- myoscore_preprocess(counts, verbose = FALSE)

  # Manual CPM for S1: total = 300; CPM = c(100/300*1e6, 200/300*1e6)
  expected_cpm_s1 <- c(100 / 300 * 1e6, 200 / 300 * 1e6)
  expected_log2_s1 <- log2(expected_cpm_s1 + 1)

  expect_equal(unname(result[, 1]), expected_log2_s1, tolerance = 1e-6)
})

test_that("myoscore_preprocess rejects negative values", {
  counts <- matrix(c(-1, 2, 3, 4), nrow = 2)
  rownames(counts) <- c("A", "B")
  colnames(counts) <- c("S1", "S2")

  expect_error(myoscore_preprocess(counts, verbose = FALSE), "negative")
})

test_that("myoscore_preprocess handles data.frame input", {
  df <- data.frame(
    S1 = c(100, 200),
    S2 = c(300, 400),
    row.names = c("GeneA", "GeneB")
  )

  result <- myoscore_preprocess(df, verbose = FALSE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
})
