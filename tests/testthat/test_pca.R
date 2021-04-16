context("PCA")

pca <- prcomp(USArrests)
test_that("getPcaVars", {
    expect_equal(getPcaVars(pca)[1], 96.6)
    expect_equal(getPcaVars(pca, digits = 5)[2], 2.7817)
})
