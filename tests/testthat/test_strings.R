context('String wrappers')

info <- c('DLPFC_polyA', 'DLPFC_RiboZero')
test_that('ss', {
    expect_equal(unique(ss(info, '_', 1)), 'DLPFC')
    expect_equal(ss(info, '_', 2)[1], 'polyA')
})


test_that('splitit and split0', {
    expect_equal(splitit(letters[1:2]), list(a = 1, b = 2))
    expect_equal(splitit(letters), split0(letters))
})