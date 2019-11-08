test_that("correct plot_hexbin_density Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(
        plot_hexbin_density(pbmc_small))[2], "ggplot")
})

test_that("error no hexbin plot_hexbin_Seurat", {
    expect_error(plot_hexbin_density(pbmc_small))
})

test_that("correct plot_hexbin_density SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class({plot_hexbin_density(pbmc_small)})[2], "ggplot")
})

test_that("error no hexbin plot_hexbin_density SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    expect_error(plot_hexbin_density(pbmc_small))
})
