test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_density(pbmc_small)})[2], "ggplot")
})

test_that("error no hexbin", {
  expect_error(plot_hexbin_density(pbmc_small))
})
