test_that("correct plot_hexbin_meta_plus SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    pbmc_small$random1 <- rnorm(ncol(pbmc_small))
    expect_equal(class(plot_hexbin_meta_plus(pbmc_small,
        col1 = "random",
        col2 = "random1", action = "mean"
    ))[2], "ggplot")
})

test_that("error plot_hexbin_meta_plus SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    pbmc_small$random1 <- rnorm(ncol(pbmc_small))
    expect_error(plot_hexbin_meta_plus(pbmc_small,
        col1 = "RNA_snn_res.1",
        col2 = "nCount_RNA", action = "mean"
    ))
})
