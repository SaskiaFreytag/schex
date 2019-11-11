test_that("correct make_hexbin_label Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(make_hexbin_label(pbmc_small, col="RNA_snn_res.1")),
        "data.frame")
})

test_that("error make_hexbin_label Seurat", {
    expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.1"))
})

test_that("error wrong col make_hexbin_label Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.2"))
})

test_that("correct make_hexbin_label SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(make_hexbin_label(pbmc_small, col="RNA_snn_res.1")),
        "data.frame")
})

test_that("error make_hexbin_label SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.1"))
})

test_that("error wrong col make_hexbin_label SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(make_hexbin_label(pbmc_small, col="RNA_snn_res.2"))
})

