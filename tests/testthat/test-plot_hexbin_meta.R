test_that("correct class", {
  expect_equal(class({pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA");
  plot_hexbin_meta(pbmc_small, col="RNA_snn_res.1", action="prop", no=1)})[2], "ggplot")
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

test_that("correct plot_hexbin", {
  pbmc_small <- make_hexbin(pbmc_small, 10, dimension_reduction = "PCA")
  drhex <- pbmc_small@misc$hexbin$hexbin.matrix
  drhex <- as_tibble(drhex)
  drhex$lab_majority <- c(rep(c("A", "B"), dim(drhex)[1]/2))
  expect_equal(class(schex:::.plot_hexbin(drhex, colour_by="lab_majority"))[2],
               "ggplot")
})
