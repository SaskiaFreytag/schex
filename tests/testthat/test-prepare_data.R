test_that("correct .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_equal(class(schex:::.prepare_data_feature(pbmc_small,
        mod = "ADT",
        type = "counts", feature = "A1"
    )), "numeric")
})

test_that("error no mod .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_error(schex:::.prepare_data_feature(pbmc_small,
        mod = "ADF",
        type = "counts", feature = "A1"
    ))
})


test_that("error type .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_error(schex:::.prepare_data_feature(pbmc_small,
        mod = "ADT",
        type = "logcounts", feature = "A1"
    ))
})

test_that("error feature .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_error(schex:::.prepare_data_feature(pbmc_small,
        mod = "ADT",
        type = "counts", feature = "A11"
    ))
})


test_that("correct .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(schex:::.prepare_data_feature(pbmc_small,
        mod = "RNA",
        type = "counts", feature = "Gene_0001"
    )), "numeric")
})


test_that("error type .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(schex:::.prepare_data_feature(pbmc_small,
        mod = "RNA",
        type = "klcounts", feature = "Gene_0001"
    ))
})

test_that("error feature RNA .prepare_data_feature SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(schex:::.prepare_data_feature(pbmc_small,
        mod = "RNA",
        type = "counts", feature = "MS4A12"
    ))
})


test_that("correct .prepare_data_meta SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_equal(class(schex:::.prepare_data_meta(pbmc_small,
        col = "random"
    )), "factor")
})


test_that("error .prepare_data_meta SingleCellExperiment", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    pbmc_small$random <- factor(sample(1:3, ncol(pbmc_small), replace = TRUE))
    expect_error(schex:::.prepare_data_meta(pbmc_small, col = "ran"))
})
