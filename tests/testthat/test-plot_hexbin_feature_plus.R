test_that("correct plot_hexbin_feature_plus SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_equal(class(plot_hexbin_feature_plus(pbmc_small,
        col = "random",
        type = "counts", feature = "Gene_0001", action = "mean"
    ))[2], "ggplot")
})

test_that("error plot_hexbin_feature_plus SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    expect_error(plot_hexbin_feature_plus(pbmc_small,
        col = "RNA_snn_res.1",
        type = "counts", feature = "NRBP1", action = "mean"
    ))
})
