test_that("correct class", {
  expect_equal(class(make_hexbin(pbmc_small, 10,
      dimension_reduction = "PCA"))[1], "Seurat")
})

test_that("done computing", {
  expect_equal(length(make_hexbin(pbmc_small, 10,
            dimension_reduction = "PCA")@misc$hexbin[[1]]), 80)
})


