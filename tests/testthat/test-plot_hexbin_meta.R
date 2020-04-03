test_that("correct plot_hexbin_meta prop Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="prop", no=1))[2], "ggplot")
})

test_that("correct plot_hexbin_meta majority Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="majority"))[2], "ggplot")
})

test_that("error plot_hexbin_meta no hexbin Seurat", {
    expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="prop", no=1))
})

test_that("error plot_hexbin_meta wrong col Seurat", {
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_meta(pbmc_small,
        col="RNA_snn_res.2", action="prop", no=1))
})

test_that("error plot_hexbin", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  expect_error(schex:::.plot_hexbin(drhex, colour_by = "Cluster_majority"))
})

test_that("correct plot_hexbin majority", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_majority <- c(rep(c("A", "B"), dim(drhex)[1] / 2))
  expect_equal(
    class(.plot_hexbin(drhex, colour_by = "lab_majority",
                               action="majority"))[2],
    "ggplot"
  )
})

test_that("correct plot_hexbin prop", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_prop_2 <- c(rep(c("A", "B"), dim(drhex)[1] / 2))
  expect_equal(
    class(.plot_hexbin(drhex, colour_by = "lab_prop_2",
                       action="prop"))[2],
    "ggplot"
  )
})

test_that("correct plot_hexbin numeric", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_mean <- c(rep(c(1, 2), dim(drhex)[1] / 2))
  expect_equal(
    class(.plot_hexbin(drhex, colour_by = "lab_mean",
        action="mean"))[2],
    "ggplot"
  )
})

test_that("correct plot_hexbin_meta class prop SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="prop", no=1))[2], "ggplot")
})

test_that("correct plot_hexbin_meta class majority SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="majority"))[2], "ggplot")
})

test_that("error plot_hexbin_meta no hexbin SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
        action="prop", no=1))
})

test_that("error plot_hexbin_meta wrong col SingleCellExperiment", {
    pbmc_small <- as.SingleCellExperiment(pbmc_small)
    pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
    expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.2",
        action="prop", no=1))
})
