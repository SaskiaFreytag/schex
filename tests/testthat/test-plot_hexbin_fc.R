test_that("correct plot_hexbin_fc SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$test <- as.factor(sample(1:2, dim(pbmc_small)[2], replace = TRUE))
    expect_equal(class(plot_hexbin_fc(pbmc_small,
        col = "test", feature = "Gene_0001",
        type = "counts"
    ))[2], "ggplot")
})

test_that("error plot_hexbin_fc SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$test <- as.factor(sample(1:2, dim(pbmc_small)[2], replace = TRUE))
    expect_error(plot_hexbin_fc(pbmc_small,
        col = "test", feature = "CA2",
        type = "counts"
    ))
})
