test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_density(pbmc_small)})[2], "ggplot")
})

test_that("error no hexbin", {
  expect_error(plot_hexbin_density(pbmc_small))
})

test_that("correct class sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class({plot_hexbin_density(pbmc_small)})[2], "ggplot")
})

test_that("error no hexbin sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(plot_hexbin_density(pbmc_small))
})
