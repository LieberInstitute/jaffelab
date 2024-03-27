context("String wrappers")

info <- c("DLPFC_polyA", "DLPFC_RiboZero")
test_that("ss", {
    expect_equal(unique(ss(info, "_", 1)), "DLPFC")
    expect_equal(ss(info, "_", 2)[1], "polyA")
})


test_that("splitit and split0", {
    expect_equal(splitit(letters[seq_len(2)]), list(a = 1, b = 2))
    expect_equal(splitit(letters), split0(letters))
})

test_that("ucsc_to_grances", {
    expect_equal(ucsc_to_granges(c("chr1:1000-2000", "chrY:1-100")), GRanges(c("chr1", "chrY"), IRanges(c(1000, 1), c(2000, 100))))
    expect_equal(c("chr1:1000-2000", "chrY:1-100"), granges_to_ucsc(ucsc_to_granges(c("chr1:1000-2000", "chrY:1-100"))))
})
