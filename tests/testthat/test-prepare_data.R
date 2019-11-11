test_that("correct .prepare_data_feature Seurat", {
    protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
    rownames(protein) <- paste0("A", seq(1,10,1))
    colnames(protein) <- colnames(pbmc_small)
    pbmc_small[["ADT"]] <- CreateAssayObject(counts = protein)
    expect_equal(class(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
        type="counts", feature="A1")), "numeric")
})

test_that("correct .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  alt_adt <- SummarizedExperiment(assays=list(counts=protein))
  altExp(pbmc_small, "ADT") <- alt_adt
  expect_equal(class(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
        type="counts", feature="A1")), "numeric")
})

test_that("error no mod .prepare_data_feature Seurat", {
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  pbmc_small[["ADT"]] <- CreateAssayObject(counts = protein)
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADF", 
        type="counts", feature="A1"))
})

test_that("error no mod .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  alt_adt <- SummarizedExperiment(assays=list(counts=protein))
  altExp(pbmc_small, "ADT") <- alt_adt
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADF", 
       type="counts", feature="A1"))
})

test_that("error type .prepare_data_feature Seurat", {
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  pbmc_small[["ADT"]] <- CreateAssayObject(counts = protein)
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
        type="logcounts", feature="A1"))
})

test_that("error type .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  alt_adt <- SummarizedExperiment(assays=list(counts=protein))
  altExp(pbmc_small, "ADT") <- alt_adt
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
        type="logcounts", feature="A1"))
})

test_that("error feature .prepare_data_feature Seurat", {
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  pbmc_small[["ADT"]] <- CreateAssayObject(counts = protein)
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
      type="counts", feature="A11"))
})

test_that("error feature .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  protein <- matrix(rnorm(10* ncol(pbmc_small)), ncol=ncol(pbmc_small))
  rownames(protein) <- paste0("A", seq(1,10,1))
  colnames(protein) <- colnames(pbmc_small)
  alt_adt <- SummarizedExperiment(assays=list(counts=protein))
  altExp(pbmc_small, "ADT") <- alt_adt
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="ADT", 
        type="counts", feature="A11"))
})

test_that("correct .prepare_data_feature Seurat", {
  expect_equal(class(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
        type="counts", feature="MS4A1")), "numeric")
})

test_that("correct .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
          type="counts", feature="MS4A1")), "numeric")
})

test_that("error type .prepare_data_feature Seurat", {
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
      type="klcounts", feature="MS4A1"))
})

test_that("error type .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
      type="klcounts", feature="MS4A1"))
})

test_that("error feature RNA .prepare_data_feature Seurat", {
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
      type="counts", feature="MS4A12"))
})

test_that("error feature RNA .prepare_data_feature SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(schex:::.prepare_data_feature(pbmc_small, mod="RNA", 
       type="counts", feature="MS4A12"))
})

test_that("correct .prepare_data_meta Seurat", {
  expect_equal(class(schex:::.prepare_data_meta(pbmc_small, 
      col="RNA_snn_res.0.8")), "factor")
})

test_that("correct .prepare_data_meta SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_equal(class(schex:::.prepare_data_meta(pbmc_small, 
      col="RNA_snn_res.0.8")), "factor")
})

test_that("error .prepare_data_meta Seurat", {
  expect_error(schex:::.prepare_data_meta(pbmc_small, col="RNA_snn_"))
})

test_that("error .prepare_data_meta SingleCellExperiment", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(schex:::.prepare_data_meta(pbmc_small, col="RNA_snn_"))
})

