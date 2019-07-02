test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
    make_hexbin_label(pbmc_small, col="RNA_snn_res.1")}), "data.frame")
})

test_that("error no hexbin", {
  expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.1"))
})

test_that("error wrong col", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.2"))
})
