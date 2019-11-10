test_that("correct class", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO1", action = "prop_0"
  ))[2], "ggplot")
})

test_that("error gene not found", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO13", action = "prop_0"
  ))
})

test_that("error assay not found", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_gene(pbmc_small,
    type = "countstt",
    gene = "TALDO1", action = "prop_0"
  ))
})

test_that("error assay no hexbin", {
  expect_error(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO1", action = "prop_0"
  ))
})

test_that("correct class sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO1", action = "prop_0"
  ))[2], "ggplot")
})

test_that("error gene not found sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(class(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO13", action = "prop_0"
  )))
})

test_that("error assay not found sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(class(plot_hexbin_gene(pbmc_small,
    type = "countsttt",
    gene = "TALDO1", action = "prop_0"
  )))
})

test_that("error assay no hexbin sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(plot_hexbin_gene(pbmc_small,
    type = "counts",
    gene = "TALDO1", action = "prop_0"
  ))
})
