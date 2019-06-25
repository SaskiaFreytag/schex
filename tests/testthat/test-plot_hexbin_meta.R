test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1)})[2], "ggplot")
})
