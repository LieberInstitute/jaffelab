context('String wrappers')

info <- c('DLPFC_polyA', 'DLPFC_RiboZero')
test_that('ss', {
    expect_equal(unique(ss(info, '_', 1)), 'DLPFC')
    expect_equal(ss(info, '_', 2)[1], 'polyA')
})
