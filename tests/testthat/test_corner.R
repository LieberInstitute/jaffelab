context("corner")

lis <- list(iris, mtcars, matrix(rnorm(1000), ncol = 10))

test_that("corner", {
    expect_equal(corner(lis), lapply(lis, function(x) head(x[, seq_len(min(6, ncol(x)))])))
})
