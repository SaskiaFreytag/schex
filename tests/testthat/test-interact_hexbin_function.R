test_that("correct spearman .interact_hexbin_function", {
    first_x <- runif(64)
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_equal(class(schex:::.interact_hexbin_function(
        first_x, second_x,
        "corr_spearman", cID
    )), "numeric")
})

test_that("correct mi .interact_hexbin_function", {
    first_x <- runif(64)
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_equal(class(schex:::.interact_hexbin_function(
        first_x, second_x,
        "mi", cID
    )), "numeric")
})

test_that("correct fc .interact_hexbin_function", {
    first_x <- runif(64)
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_equal(class(schex:::.interact_hexbin_function(
        first_x, second_x,
        "fc", cID
    )), "numeric")
})

test_that("error spearman .interact_hexbin_function", {
    first_x <- as.character(runif(64))
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_error(schex:::.interact_hexbin_function(
        first_x, second_x,
        "corr_spearman", cID
    ))
})

test_that("error mi .interact_hexbin_function", {
    first_x <- as.character(runif(64))
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_error(schex::.interact_hexbin_function(
        first_x, second_x,
        "mi", cID
    ))
})

test_that("error fc .interact_hexbin_function", {
    first_x <- as.character(runif(64))
    second_x <- runif(64)
    cID <- sample(1:4, 64, replace = TRUE)
    expect_error(schex::.interact_hexbin_function(
        first_x, second_x,
        "fc", cID
    ))
})
