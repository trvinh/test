context("test estimation of gene ages")

test_that("test estimation of gene ages", {
    # load full processed data
    data("full_processed_profile", package="phyloprofile")

    # calculate gene ages
    rank_name <- "class"
    ref_taxon <- "Mammalia"
    var1_cutoff <- c(0,1)
    var2_cutoff <- c(0,1)
    percent_cutoff <- c(0,1)
    gene_age <- estimate_gene_age(full_processed_profile,
                                  rank_name, ref_taxon,
                                  var1_cutoff, var2_cutoff, percent_cutoff)

    expect_true(gene_age$age[gene_age$geneID == "OG_1019"] == "07_LUCA")
})

test_that("test plotting gene age plot", {
    gene_age_df <- data.frame(
        geneID = c("OG_1017", "OG_1019"),
        cat = c("0000001", "0000001"),
        age = c("07_LUCA", "07_LUCA"),
        stringsAsFactors = FALSE
    )
    plot_df <- gene_age_plotDf(gene_age_df)
    gene_age_text <- 1

    p <- create_gene_age_plot(plot_df, gene_age_text)
    expect_true(p$guides$fill$nrow == 1)
})
