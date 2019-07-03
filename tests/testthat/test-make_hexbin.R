test_that("correct class", {
  expect_equal(class(make_hexbin(pbmc_small, 10,
      dimension_reduction = "PCA"))[1], "Seurat")
})

test_that("done computing", {
  expect_equal(length(make_hexbin(pbmc_small, 10,
            dimension_reduction = "PCA")@misc$hexbin[[1]]), 80)
})

test_that("error dimension reduction", {
  expect_error(make_hexbin(pbmc_small, 10,
           dimension_reduction = "UMAP"))
})

test_that("error class", {
  sce <- c(1,2,3)
  expect_error(make_hexbin(sce, 10,
          dimension_reduction = "PCA"))
})

test_that("correct class sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(make_hexbin(pbmc_small, 10,
      dimension_reduction = "PCA"))[1], "SingleCellExperiment")
})

test_that("done computing sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(length(make_hexbin(pbmc_small, 10,
      dimension_reduction = "PCA")@metadata$hexbin[[1]]), 80)
})

test_that("error dimension reduction sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(make_hexbin(pbmc_small, 10,
                           dimension_reduction = "UMAP"))
})
