context("test creating data for plotting profiles")

test_that("test profile plot data generation", {
    data("full_processed_profile", package="phyloprofile")
    
    plot_df <- data_main_plot(full_processed_profile)
    expect_true(nrow(plot_df) == 20)
    
    selected_taxa <- c("Mammalia", "Echinoidea", "Gunneridae")
    selected_seq <- "all"
    customized_plot_df <- data_customized_plot(
        full_processed_profile, selected_taxa, selected_seq
    )
    expect_true(nrow(customized_plot_df) == 10)
})
