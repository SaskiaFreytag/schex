test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_gene(pbmc_small, type="counts", gene="TALDO1", action="prop_0")})[2], "ggplot")
})
