context('Smaller functions')

test_that('getOR', {
    expect_equal(getOR(matrix(c(75, 75, 100, 150), ncol = 2)), 1.5)
})