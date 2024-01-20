test_that("correct class SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    expect_equal(class(make_hexbin(pbmc_small, 10,
        dimension_reduction = "PCA"
    ))[1], "SingleCellExperiment")
})

test_that("correct done SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    expect_equal(length(make_hexbin(pbmc_small, 10,
        dimension_reduction = "PCA"
    )@metadata$hexbin[[1]]), 200)
})

test_that("error dimension reduction SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    expect_error(make_hexbin(pbmc_small, 10,
        dimension_reduction = "UMAP"
    ))
})
