test_that("correct class sce", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_interact(pbmc_small,
        type = c("counts", "counts"), mod = c("RNA", "ADT"),
        feature = c("Gene_0001", "A1"),
        interact = "mi"
    ))[2], "ggplot")
})

test_that("correct class fc sce", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_interact(pbmc_small,
        type = c("counts", "counts"), mod = c("RNA", "RNA"),
        feature = c("Gene_0001", "Gene_0002"), interact = "fc"
    ))[2], "ggplot")
})

test_that("error feature not found sce", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_interact(pbmc_small,
        type = c("counts", "counts"),
        mod = c("RNA", "ADT"), feature = c("Gene_001", "A1"), interact = "mi"
    ))
})

test_that("error assay not found sce", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_interact(pbmc_small,
        type = c("counts", "counts"),
        mod = c("RNA", "ADA"), feature = c("CD7", "A1"), interact = "mi"
    ))
})

test_that("error no hexbin sce", {
    pbmc_small <- mockSCE()
    pbmc_small <- logNormCounts(pbmc_small)
    pbmc_small <- runPCA(pbmc_small)
    protein <- matrix(rnorm(10 * ncol(pbmc_small)), ncol = ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1, 10, 1))
    colnames(protein) <- colnames(pbmc_small)
    alt_adt <- SummarizedExperiment(assays = list(counts = protein))
    altExp(pbmc_small, "ADT") <- alt_adt
    expect_error(plot_hexbin_interact(pbmc_small,
        type = c("counts", "counts"),
        mod = c("RNA", "ADT"), feature = c("CD7", "A1"), interact = "mi"
    ))
})
