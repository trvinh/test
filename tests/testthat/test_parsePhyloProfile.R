context("test parsing and pre-processing phyloprofile input")

test_that("test connection to taxonomy files and getting input taxa", {
    rank_name <- "species"
    input_df <- phyloprofile::create_long_matrix("main_wide_test.txt")

    input_taxonID <- phyloprofile::get_input_taxa_id(input_df)
    expect_true(length(input_taxonID) == 4)

    input_taxonName <- phyloprofile::get_input_taxa_name(rank_name,
                                                         input_taxonID)
    expect_true(nrow(input_taxonName) == 4)
})
