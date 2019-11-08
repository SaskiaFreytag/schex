test_that("correct plot_hexbin_meta_shiny Seurat", {
    expect_equal(class(plot_hexbin_meta_shiny(pbmc_small, col="RNA_snn_res.1", 
        action="majority", 
        min_nbins=2, max_nbins=10, dimension_reduction="PCA")), "shiny.appobj")
})

test_that("correct plot_hexbin_meta_shiny SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(plot_hexbin_meta_shiny(pbmc_small, col="RNA_snn_res.1", 
      action="majority", 
      min_nbins=2, max_nbins=10, dimension_reduction="PCA")), "shiny.appobj")
})
