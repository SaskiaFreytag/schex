test_that("correct plot_hexbin_meta_shiny Seurat", {
    expect_equal(class(plot_hexbin_meta_shiny(pbmc_small, col="RNA_snn_res.1", 
        action="majority", 
        min_nbins=2, max_nbins=10, dimension_reduction="PCA")), "shiny.appobj")
})

test_that("correct plot_hexbin_meta_shiny SingleCellExperiment", {
  pbmc_small <- mockSCE() 
  pbmc_small <- logNormCounts(pbmc_small) 
  pbmc_small <- runPCA(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace=TRUE))
  expect_equal(class(plot_hexbin_meta_shiny(pbmc_small, col="random", 
      action="majority", 
      min_nbins=2, max_nbins=10, dimension_reduction="PCA")), "shiny.appobj")
})
