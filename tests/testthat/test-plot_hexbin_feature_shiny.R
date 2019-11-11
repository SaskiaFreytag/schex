test_that("correct plot_hexbin_feature_shiny Seurat", {
  expect_equal(class(plot_hexbin_feature_shiny(pbmc_small, 
        type="counts", feature="TALDO1", 
        action="prop_0", min_nbins=2, max_nbins=10, dimension_reduction="PCA",
        mod="RNA")), "shiny.appobj")
})

test_that("correct plot_hexbin_feature_shiny SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(plot_hexbin_feature_shiny(pbmc_small, 
      type="counts", feature="TALDO1", 
      action="prop_0", min_nbins=2, max_nbins=10, dimension_reduction="PCA",
      mod="RNA")), "shiny.appobj")
})



