test_that("correct plot_hexbin_feature_plus Seurat", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_feature_plus(pbmc_small, col="RNA_snn_res.1", 
        type="counts", feature="NRBP1",action="mean"))[2], "ggplot")
})

test_that("correct plot_hexbin_feature_plus SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_feature_plus(pbmc_small, col="RNA_snn_res.1", 
      type="counts", feature="NRBP1",action="mean"))[2], "ggplot")
})

test_that("error plot_hexbin_feature_plus Seurat", {
  expect_error(plot_hexbin_feature_plus(pbmc_small, col="RNA_snn_res.1", 
     type="counts", feature="NRBP1",action="mean"))
})

test_that("error plot_hexbin_feature_plus SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(plot_hexbin_feature_plus(pbmc_small, col="RNA_snn_res.1", 
     type="counts", feature="NRBP1",action="mean"))
})


