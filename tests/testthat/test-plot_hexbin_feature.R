test_that("correct class plot_hexbin_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_feature(pbmc_small,
        type = "counts", mod = "ADT",
        feature = "A1", action = "prop_0"
    ))[2], "ggplot")
})

test_that("error feature not found plot_hexbin_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_feature(pbmc_small,
        type = "counts", mod = "ADT",
        feature = "A11", action = "prop_0"
    ))
})

test_that("error assay not found plot_hexbin_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_feature(pbmc_small,
        type = "counts", mod = "ADA",
        feature = "A1", action = "prop_0"
    ))
})

test_that("error no hexbin plot_hexbin_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_error(plot_hexbin_feature(pbmc_small,
        type = "counts", mod = "ADT",
        feature = "A1", action = "prop_0"
    ))
})

test_that("correct gene plot_hexbin_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_feature(pbmc_small,
        type = "counts",
        feature = "Gene_0001", action = "prop_0"
    ))[2], "ggplot")
})
