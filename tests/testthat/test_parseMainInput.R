context("test parsing main input file into long dataframe")

test_that("test fasta input", {
    a <- phyloprofile::create_long_matrix("main_fasta_test.txt")
    expect_true(ncol(a) == 5)
})

test_that("test wide input", {
    a <- phyloprofile::create_long_matrix("main_wide_test.txt")
    expect_true(nrow(a) == 8)
})

test_that("test xml input", {
    a <- phyloprofile::create_long_matrix("main_orthoxml_test.xml")
    expect_true(ncol(a) == 6)
})
