context("test parsing domain input file(s) into dataframe")

test_that("test with single file", {
    a <- parse_domain_input(NULL, "domains/OG_1029.domains", "file")
    expect_true(ncol(a) == 5)
})

test_that("test with domain folder", {
    a <- parse_domain_input("OG_1029", "domains", "folder")
    expect_true(ncol(a) == 5)

    b <- parse_domain_input("OG_1020", "domains", "folder")
    expect_true(b == "noFileInFolder")
})
