test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
    make_hexbin_label(pbmc_small, col="RNA_snn_res.1")}), "data.frame")
})
