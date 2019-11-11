test_that("correct .extract_cID Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(.extract_cID(pbmc_small)), "integer")
})

test_that("correct .extract_cID SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(.extract_cID(pbmc_small)), "integer")
})

test_that("correct .extract_hexbin Seurat", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(.extract_hexbin(pbmc_small)), "matrix")
})

test_that("correct .extract_hexbin SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(.extract_hexbin(pbmc_small)), "matrix")
})
