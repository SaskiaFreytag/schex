test_that("correct .make_hexbin_fc_function", {
    x <- as.factor(rep(1:2, 16))
    x_gene <- runif(32)
    cID <- rep(1:4, 8)
    expect_equal(
        class(schex:::.make_hexbin_fc_function(x, x_gene, cID)),
        "list"
    )
})

test_that("correct1 .make_hexbin_fc_function", {
    x <- as.factor(rep(1:2, 32))
    x_gene <- runif(64)
    cID <- sample(1:2, 64, replace = TRUE)
    expect_equal(class(schex:::.make_hexbin_fc_function(x, x_gene, cID)), "list")
})

test_that("error .make_hexbin_fc_function", {
    x <- rep(1:2, 16)
    x_gene <- runif(32)
    cID <- rep(1:4, 8)
    expect_error(schex:::.make_hexbin_fc_function(x, x_gene, cID))
})
