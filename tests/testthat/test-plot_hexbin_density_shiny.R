test_that("correct plot_hexbin_density_shiny Seurat", {
  expect_equal(class(plot_hexbin_density_shiny(pbmc_small,3, 10, 
      dimension_reduction = "PCA")), "shiny.appobj")
})

test_that("correct plot_hexbin_density_shiny SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(plot_hexbin_density_shiny(pbmc_small,3, 10, 
      dimension_reduction = "PCA")), "shiny.appobj")
})



