context("test creating domain plot for a pair seed and ortholog proteins")

test_that("test domain plot", {
    info <- c("OG_1029", "E.intestinalis@876142@Eint_030020")
    domain_df <- phyloprofile::parse_domain_input(
        "OG_1029","domains/OG_1029.domains","file"
    )
    expect_true(nrow(domain_df) == 8)

    label_size <- 12
    title_size <- 12
    p <- phyloprofile::create_archi_plot(info, domain_df,
                                         label_size, title_size)
    expect_true(nrow(p$layout) == 2)
})
