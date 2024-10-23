
test_that("PSO algorithm works", {
  deg_list <- list(c("GeneA", "GeneB"), c("GeneC", "GeneD"))
  result <- run_pso(deg_list)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("GMM clustering works", {
  data <- matrix(rnorm(100), nrow = 10)
  result <- run_gmm(data, n_clusters = 3)
  expect_true("model" %in% names(result))
  expect_true("clusters" %in% names(result))
})
