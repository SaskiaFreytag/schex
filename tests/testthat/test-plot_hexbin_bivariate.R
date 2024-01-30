test_that("error feature not found plot_hexbin_bivariate SingleCellExperiment", {
  pbmc_small <- mockSCE() 
  pbmc_small <- logNormCounts(pbmc_small) 
  pbmc_small <- runPCA(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_bivariate(pbmc_small, type="counts", feature="TALDO0", fan = TRUE))
})

test_that("error mod plot_hexbin_bivariate SingleCellExperiment", {
  pbmc_small <- mockSCE() 
  pbmc_small <- logNormCounts(pbmc_small) 
  pbmc_small <- runPCA(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_bivariate(pbmc_small, type="counts", 
                                     mod="ADT", feature="ENSG00000188976", fan = TRUE))
})

test_that("error no hexbin plot_hexbin_bivariate SingleCellExperiment", {
  pbmc_small <- mockSCE() 
  pbmc_small <- logNormCounts(pbmc_small) 
  pbmc_small <- runPCA(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_bivariate(pbmc_small, type="counts", feature="ENSG00000188976", fan = TRUE))
})