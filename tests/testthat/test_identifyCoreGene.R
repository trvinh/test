context("test identification of core genes for a selected set of taxa")

test_that("test core gene estimation", {
    # load full processed data
    data("full_processed_profile", package="phyloprofile")

    # calculate gene ages
    rank_name <- "class"
    taxa_core <- c("Mammalia", "Mucorales", "Alphaproteobacteria")
    var1_cutoff <- c(0,1)
    var2_cutoff <- c(0,1)
    percent_cutoff <- c(0,1)
    core_coverage <- 100
    core_gene <- get_core_gene(rank_name,
                               taxa_core,
                               full_processed_profile,
                               var1_cutoff, var2_cutoff,
                               percent_cutoff, core_coverage)

    expect_true(length(core_gene) == 1)
})
