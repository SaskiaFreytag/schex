test_that("correct make_hexbin_label SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_equal(
        class(make_hexbin_label(pbmc_small, col = "random")),
        "data.frame"
    )
})

test_that("error make_hexbin_label SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    expect_error(make_hexbin_label(pbmc_small, col = "random"))
})

test_that("error wrong col make_hexbin_label SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(make_hexbin_label(pbmc_small, col = "random"))
})
