test_that("correct plot_hexbin_meta class prop SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_equal(class(plot_hexbin_meta(pbmc_small,
        col = "random",
        action = "prop", no = 1
    ))[2], "ggplot")
})

test_that("correct plot_hexbin_meta class majority SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_equal(class(plot_hexbin_meta(pbmc_small,
        col = "random",
        action = "majority"
    ))[2], "ggplot")
})

test_that("error plot_hexbin_meta no hexbin SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    expect_error(plot_hexbin_meta(pbmc_small,
        col = "RNA_snn_res.1",
        action = "prop", no = 1
    ))
})

test_that("error plot_hexbin_meta wrong col SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_meta(pbmc_small,
        col = "RNA_snn_res.2",
        action = "prop", no = 1
    ))
})
