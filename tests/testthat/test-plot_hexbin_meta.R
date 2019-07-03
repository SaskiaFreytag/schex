test_that("correct class prop", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1)})[2], "ggplot")
})

test_that("correct class majority", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="majority")})[2], "ggplot")
})

test_that("error no hexbin", {
  expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1))
})

test_that("error wrong col", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.2", action="prop", no=1))
})

test_that("error plot_hexbin", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  expect_error(schex:::.plot_hexbin(drhex, colour_by="Cluster_majority"))
})

test_that("correct plot_hexbin majority", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_majority <- c(rep(c("A", "B"), dim(drhex)[1]/2))
  expect_equal(class(schex:::.plot_hexbin(drhex, colour_by="lab_majority"))[2],
               "ggplot")
})

test_that("correct plot_hexbin prop", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_prop_2 <- c(rep(c("A", "B"), dim(drhex)[1]/2))
  expect_equal(class(schex:::.plot_hexbin(drhex, colour_by="lab_prop_2"))[2],
               "ggplot")
})

test_that("correct plot_hexbin numeric", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_mean <- c(rep(c(1, 2), dim(drhex)[1]/2))
  expect_equal(class(schex:::.plot_hexbin(drhex, colour_by="lab_mean"))[2],
               "ggplot")
})



test_that("correct class prop sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
      action="prop", no=1))[2], "ggplot")
})

test_that("correct class majority sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_equal(class(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
    action="majority"))[2], "ggplot")
})

test_that("error no hexbin sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1",
    action="prop", no=1))
})

test_that("error wrong col sce", {
  pbmc_small <- as.SingleCellExperiment(pbmc_small)
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  expect_error(plot_hexbin_meta(pbmc_small, col="RNA_snn_res.2",
          action="prop", no=1))
})


