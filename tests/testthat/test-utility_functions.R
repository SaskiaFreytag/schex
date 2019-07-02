test_that("get_mode works", {
  expect_equal(schex:::.get_mode(c(2,2,2,2,1,4,5)), 2)
})

test_that("make_hexbin_colnames works", {
  expect_equal(schex:::.make_hexbin_colnames(c("A", "B"), "col"),
   c("col_prop_A", "col_prop_B"))
})

test_that("make_hexbin_function majority", {
  x <- as.factor(rep(c("A", "A", "B"), 3))
  action <- "majority"
  cID <- c(1,1,1,2,2,2,3,3,3)
  expect_equal(levels(schex:::.make_hexbin_function(x, action, cID)),
               "A")
})

test_that("make_hexbin_function mode", {
  x <- rep(c(2,2,1), 3)
  action <- "mode"
  cID <- c(1,1,1,2,2,2,3,3,3)
  expect_equal(schex:::.make_hexbin_function(x, action, cID),
               c(2,2,2))
})

test_that("make_hexbin_function prop", {
  x <-  as.factor(rep(c("A", "A", "B"), 3))
  action <- "prop"
  cID <- c(1,1,1,2,2,2,3,3,3)
  expect_equal(class(schex:::.make_hexbin_function(x, action, cID)),
               "matrix")
})

test_that("make_hexbin_function mean", {
  x <- rep(c(2,3,1), 3)
  action <- "mean"
  cID <- c(1,1,1,2,2,2,3,3,3)
  expect_equal(schex:::.make_hexbin_function(x, action, cID),
               c(2,2,2))
})

test_that("make_hexbin_function median", {
  x <- rep(c(2,3,1), 3)
  action <- "mean"
  cID <- c(1,1,1,2,2,2,3,3,3)
  expect_equal(schex:::.make_hexbin_function(x, action, cID),
               c(2,2,2))
  })
