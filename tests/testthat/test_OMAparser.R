context("test functions for getting data from OMA browser")

test_that("test for a wrong OMA ID", {
    id <- "RATNOBLABLA"
    expect_error(check_oma_id(id))
})

test_that("test for getting HOG orthologs", {
    id <- "HUMAN29397"
    orthologs <- get_oma_members(id, "HOG")
    expect_true(length(orthologs) > 1)
})

test_that("test for getting protein annotations and fasta for one OMA id", {
    id <- "HUMAN29397"
    oma_df <- get_data_for_one_ortholog(id)
    expect_true(ncol(oma_df) == 5)
})

test_that("test for parsing annotation", {
    id <- "HUMAN29397"
    oma_df <- get_data_for_one_oma(id, "HOG")
    domain_df <- get_all_domains_oma(oma_df)
    expect_true(ncol(domain_df) == 6)
})


