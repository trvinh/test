context("test creating data and plots for distritbution analysis function")

test_that("test distritbution analyze functions", {
    # raw input
    input_data <- phyloprofile::create_long_matrix("main_wide_test.txt")

    # selected rank name
    rank_name <- "species"

    # variable thresholds
    var1_cutoff_min <- 0.0
    var1_cutoff_max <- 1.0
    var2_cutoff_min <- 0.0
    var2_cutoff_max <- 1.0

    # distribution data for percentage of species present in supertaxa
    percent_distribution_data <- create_percentage_distribution_data(
        input_data, rank_name
    )
    expect_true(nrow(percent_distribution_data) == 6)

    # distribution data for 2 additional variables
    distribution_data <- create_variable_distribution_data(
        input_data,
        var1_cutoff_min,
        var1_cutoff_max,
        var2_cutoff_min,
        var2_cutoff_max
    )
    expect_true(nrow(distribution_data) == 6)

    # distribution data for 2 additional variables of a subset of taxa
    input_taxonID <- phyloprofile::get_input_taxa_id(input_data)
    input_taxonName <- phyloprofile::get_input_taxa_name(rank_name,
                                                         input_taxonID)
    ref_taxon <- input_taxonName$fullName[1]
    taxa_tree <- NULL

    sorted_taxa <- phyloprofile::sort_input_taxa(input_taxonID,
                                                 input_taxonName,
                                                 rank_name,
                                                 ref_taxon,
                                                 taxa_tree)

    full_profile_data <- phyloprofile::parse_info_profile(input_data,
                                                          sorted_taxa,
                                                          "max",
                                                          "mean")
    selected_genes <- c("OG_1017", "OG_1019")
    selected_taxa <- c("Arabidopsis thaliana", "Encephalitozoon intestinalis")
    subset_distribution_data <- create_variable_distribution_data_subset(
        full_profile_data,
        distribution_data,
        selected_genes,
        selected_taxa
    )
    expect_true(ncol(subset_distribution_data) == 5)

    # plot distribution of var1
    p <- create_var_dist_plot(distribution_data,
                              "variable 1 name",
                              "var1",
                              NULL,
                              12)
    expect_true(nrow(p$data) == 6)
})
