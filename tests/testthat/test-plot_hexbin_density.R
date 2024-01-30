test_that("correct plot_hexbin_density SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class({
        plot_hexbin_density(pbmc_small)
    })[2], "ggplot")
})

test_that("error no hexbin plot_hexbin_density SingleCellExperiment", {
    pbmc_small <- mockSCE()
    expect_error(plot_hexbin_density(pbmc_small))
})
