#' Get the dataframe with the significant genes
#' @export
#' @usage get_significant_genes(in_group, selected_genes_list, rank, var,
#'     use_common_ancestor, reference_taxon, parameters, ncbi_id_list, 
#'     data_full, right_format_features, domains)
#' @param in_group list of taxa
#' @param selected_genes_list list of genes
#' @param rank selected taxonamy rank
#' @param var variable for which to  calculate the significance
#' @param use_common_ancestor boolean if the common anchestor should be used
#' @param reference_taxon taxon which is used as reference
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param ncbi_id_list list of ncbi ids
#' @param data_full full processed main data
#' @param right_format_features boolean if the features have the right format
#' @param domains dataframe holding the domain input
#' @return dataframe with the significant genes
#' @importFrom stats p.adjust
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' data("full_processed_profile", package="phyloprofile")
#' in_group <- c("Aeropyrum pernix", "Agrobacterium fabrum")
#' selected_genes_list <- "OG_1017"
#' rank <- "species"
#' var <- ""
#' use_common_ancestor <- FALSE
#' reference_taxon <- "Homo sapiens"
#' ncbi_id_list <- get_input_taxa_id(main_long_raw)
#' data_full <- full_processed_profile
#' right_format_features <- TRUE
#' domain_file <- system.file(
#'     "extdata", "domain_files/OG_1017.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' domains <- parse_domain_input("OG_1017", domain_file, "file")
#' parameters <- list(
#'     "show_p_value" = TRUE,
#'     "highlight_significant" = FALSE,
#'     "significance" = 0.05,
#'     "var1_id" = "1st variable",
#'     "var2_id" = "2nd variable",
#'     "x_size_gc" = 9,
#'     "y_size_gc" = 9,
#'     "interesting_features" = NULL,
#'     "angle_gc" = 90,
#'     "legend_gc" = "none",
#'     "legend_size_gc" = 9,
#'     "p_values_size" = 9
#' )
#' get_significant_genes(
#'     in_group,
#'     selected_genes_list,
#'     rank,
#'     var,
#'     use_common_ancestor,
#'     reference_taxon,
#'     parameters,
#'     ncbi_id_list,
#'     data_full,
#'     right_format_features,
#'     domains
#' )

get_significant_genes <- function(
    in_group,
    selected_genes_list,
    rank,
    var,
    use_common_ancestor,
    reference_taxon,
    parameters,
    ncbi_id_list,
    data_full,
    right_format_features,
    domains
){
    if (is.null(in_group) | length(selected_genes_list) == 0) return()
    significance_level <- parameters$significance

    name_list <- get_name_list() # load name List
    taxa_list <- get_taxonomy_matrix(FALSE, ncbi_id_list) # load unsorted taxa

    # Get the rank and the in-group
    # if there is more than one element in the in_group -> use common anchestor
    if (use_common_ancestor == TRUE) {
        ancestor <- get_common_ancestor(in_group, rank,
                                        name_list, taxa_list, ncbi_id_list)
        if (is.null(ancestor)) return("No common ancestor found")
        in_group <- ancestor[1]
        rank <- ancestor[2]
    } else {
        # rank <- substring(rank, 4)
        rank <- rank
    }
    if (is.na(rank)) return("No common ancestor found")

    # provide the empty data frame
    if (var == "Both") {
        significant_genes_df <- data.frame(
            geneID = character(),
            in_group = I(list()),
            out_group = I(list()),
            pvalues_1 = I(list()),
            pvalues_2 = I(list()),
            features = I(list()),
            databases = I(list()))
    } else {
        significant_genes_df <- data.frame(
            geneID = character(),
            in_group = I(list()),
            out_group = I(list()),
            pvalues = I(list()),
            features = I(list()),
            databases = I(list()))
    }

    # Get the list of genes to look at
    if (is.element("all", selected_genes_list)) {
        genes <- data_full$geneID
        genes <- genes[!duplicated(genes)]
    } else {
        genes <- selected_genes_list
    }
    genes <- sort(genes)

    # Subset depending on the rank and the in_group
    selected_subset <- get_selected_subset(rank, in_group, name_list, taxa_list)
    selected_subset <- subset(
        selected_subset, !selected_subset$fullName == reference_taxon
    )

    # Check for each gene if it is significant
    for (gene in genes) {
        print(paste("Analyzing the distribution of", gene, "..."))
        # Processing the dataframes for in- and out-group
        selected_gene_df <- subset(data_full, data_full$geneID == gene)

        in_group_df <- {
            subset(
                selected_gene_df,
                selected_gene_df$abbrName %in% selected_subset$abbrName
            )
        }
        out_group_df <- {
            subset(
                selected_gene_df,
                !(selected_gene_df$abbrName %in% selected_subset$abbrName)
            )
        }
        out_group_df <- {
            subset(out_group_df, !out_group_df$fullName == reference_taxon)
        }

        # Generate and check the p_values for the gene
        pvalue <- get_p_values(in_group_df, out_group_df, var, gene, parameters)
        new_row <- data.frame(
            geneID = gene,
            in_group = NA,
            out_group = NA,
            pvalues = NA,
            features = NA
        )
        new_row$in_group <- list(in_group_df)
        new_row$out_group <- list(out_group_df)

        if (var == "Both") {
            new_row$pvalues_1 <- pvalue[1]
            new_row$pvalues_2 <- pvalue[2]
        } else {
            new_row$pvalues <- pvalue
        }

        features  <- get_features(gene, domains)
        new_row$features <- list(features)
        if (right_format_features) {
            new_row$databases <- list(get_prefix_features(features))
        }
        significant_genes_df <- rbind(significant_genes_df, new_row)
    }

    if (var == "Both") {
        significant_genes_df$pvalues_1 <- {
            p.adjust(
                significant_genes_df$pvalues_1, method = "holm",
                n = length(significant_genes_df$pvalues_1)
            )
        }

        significant_genes_df$pvalues_2 <- {
            p.adjust(
                significant_genes_df$pvalues_2, method = "holm",
                n = length(significant_genes_df$pvalues_2)
            )
        }

        significant_genes_df <-
            significant_genes_df[
                significant_genes_df$pvalues_1 <= significance_level
                | significant_genes_df$pvalues_2 <= significance_level ,
            ]

        significant_genes_df <-
            significant_genes_df[
                !is.na(significant_genes_df$pvalues_1) |
                !is.na(significant_genes_df$pvalues_1),
            ]
    } else {
        significant_genes_df$pvalues <- {
            p.adjust(
                significant_genes_df$pvalues, method = "holm",
                n = length(significant_genes_df$pvalues)
            )
        }

        significant_genes_df <- {
            significant_genes_df[
                significant_genes_df$pvalues <= significance_level,
            ]
        }

        significant_genes_df <- {
            significant_genes_df[!is.na(significant_genes_df$pvalues) ,]
        }
    }

    # return the significant genes
    if (nrow(significant_genes_df) != 0) {
        significant_genes_df$var <- var
        significant_genes_df$rank <- rank
        return(significant_genes_df)
    } else {
        return("No candidate genes found")
    }
}

#' Generate the in_group
#' @export
#' @usage get_common_ancestor(in_group, rank, name_list, selected_in_group, 
#'     ncbi_id_list)
#' @param in_group list of taxa
#' @param rank selected rank
#' @param name_list contains "ncbiID", "fullName", "rank", "parentID"
#' @param selected_in_group contains "abbrName, "ncbiID", fullName", "strain",
#' "genus"..
#' @param ncbi_id_list list of input taxon IDs
#' @return common anchestor
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' in_group <- c("Aeropyrum pernix", "Agrobacterium fabrum")
#' rank <- "species"
#' name_list_file <- system.file(
#'     "extdata", "data/taxonNamesFull.txt",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' name_list <- as.data.frame(data.table::fread(name_list_file))
#' selected_in_group <- get_taxonomy_matrix(FALSE, NULL)
#' ncbi_id_list <- c(56636, 1176649)
#' get_common_ancestor(
#'     in_group,
#'     rank,
#'     name_list,
#'     selected_in_group,
#'     ncbi_id_list
#' )

get_common_ancestor <- function(in_group,
                                rank,
                                name_list,
                                selected_in_group,
                                ncbi_id_list){

    all_ranks <- get_taxonomy_ranks()

    selected_in_group <- {
        selected_in_group[!duplicated(selected_in_group), ]
    }

    # ranks were all elements of the in_group might be in the same taxon
    # possible_ranks <- all_ranks[all_ranks >= rank]
    possible_ranks <- all_ranks[match(rank, all_ranks):(length(all_ranks)-1)]
    position <-  1
    if (length(in_group) == 1) rank <- rank #substring(rank, 4)

    # find the common ancestor of all taxa in the in_group
    while (length(in_group) > 1 & position < length(possible_ranks)) {

        current_rank <- as.character(possible_ranks[position][1])
        next_rank <- as.character(possible_ranks[position + 1][1])

        # dataframe with all elements with fitting rank
        df_in_group <- subset(name_list, name_list$rank == current_rank)

        # subset of df_in_group  with elements that belong to the in_group
        df_in_group <- subset(df_in_group, df_in_group$fullName %in% in_group)

        # get all elements which could belong to the in-group
        possible_in_group <- subset(selected_in_group,
                                    select = c(current_rank, next_rank))
        possible_in_group <- {
            possible_in_group[
                possible_in_group[,current_rank] %in% df_in_group$ncbiID,
            ]
        }
        possible_in_group <- possible_in_group[!duplicated(possible_in_group), ]

        # only consider elements that have the next higher rank
        subset_next_rank <- taxa_select_gc(next_rank, ncbi_id_list)
        subset_next_rank <- subset_next_rank[!duplicated(subset_next_rank), ]
        subset_next_rank <- {
            subset(
                subset_next_rank,
                subset_next_rank$ncbiID %in% possible_in_group[, next_rank]
            )
        }
        in_group <- subset_next_rank$fullName
        position <- position + 1
        rank <- next_rank
    }

    # Return the in-group and the rank
    if (position > length(possible_ranks)) return()
    return(c(in_group, rank))
}


#' Get the subset depending on the choosen rank
#' @param rank selected rank
#' @param in_group list of taxa
#' @param name_list contains "ncbiID", "fullName", "rank", "parentID"
#' @param taxa_list contains "abbrName, "ncbiID", fullName", "strain", "genus"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_selected_subset <- function(rank, in_group, name_list, taxa_list){
    # Look if the fullName is in the in_group
    name_list$fullName <- as.character(name_list$fullName)
    name_list_rank <- subset(name_list, name_list$rank == rank)
    in_group_subset <- subset(name_list, name_list$fullName %in% in_group)

    # Look if it has the right rank
    selected_subset <- taxa_list[taxa_list[, rank] %in% in_group_subset$ncbiID,]

    return(selected_subset)
}

#' Decide if the gene is significant
#' @param in_group contains "supertaxon", "geneID", "ncbiID", "orthoID",
#' "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param out_group  as in-group but with information containing the out-group
#' @param variable variable(s) to claculate the plots for
#' @param gene gene to calculate the p-values for
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return return the pvalues
#' @author Carla Mölbert (carla.moelbert@gmx.de)

get_p_values <- function(in_group, out_group, variable, gene, parameters){
    significance_level <- parameters$significance

    # get the p-values for both variables
    if (variable == "Both") {
        var1 <- parameters$var1
        var2 <- parameters$var2

        pvalues1 <- calculate_p_value(
            in_group$var1,
            out_group$var1,
            significance_level
        )

        pvalues2 <- calculate_p_value(
            in_group$var2,
            out_group$var2,
            significance_level
        )
        pvalues <- list(pvalues1, pvalues2)
        return(pvalues)
    }
    # get the p-values for one variable
    else {
        # Check which variable is selected and get the p_values
        if (variable == parameters$var1_id) {
            pvalues <- calculate_p_value(
                in_group$var1,
                out_group$var1,
                significance_level
            )
        } else {
            pvalues <- calculate_p_value(
                in_group$var2,
                out_group$var2,
                significance_level
            )
        }

        return(pvalues)
    }
}

#' calculate the p_values
#' @param var_in list of values for the variable concerning the in-group
#' @param var_out list of values for the variable concerning the out-group
#' @param significance_level significant cutoff for statistical test
#' @return return the pvalues
#' @importFrom stats ks.test
#' @importFrom stats wilcox.test
#' @author Carla Mölbert (carla.moelbert@gmx.de)

calculate_p_value <- function(var_in, var_out, significance_level){
    # delete all entrys that are NA
    var_in <- var_in[!is.na(var_in)]
    var_out <- var_out[!is.na(var_out)]

    # if there is no data in one of the groups the p-value is NULL
    if (length(var_in) == 0) return(NA)
    else if (length(var_out) == 0) return(NA)
    else {
        # * Kolmogorov-Smirnov Test
        # H0 : The two samples have the same distribution
        ks <- suppressWarnings(
            ks.test(unique(var_in), unique(var_out), exact = FALSE)
        )

        p_value <- ks$p.value # probabilitiy to recet H0, if it is correct

        if (p_value < significance_level) pvalue <- c(p_value)

        else {
            # * Wilcoxon-Mann-Whitney Test
            # H0: the samples have the same location parameters

            wilcox <- suppressWarnings(
                wilcox.test(var_in,
                            var_out,
                            alternative = "two.sided",
                            #exact = FALSE,
                            paired = FALSE)
            )
            p_value_wilcox <- wilcox$p.value
            # pvalue <- c(p_value, p_value_wilcox)
            pvalue <- c(p_value_wilcox)
        }

        # perm <- jmuOutlier:: perm.test(
        #     unique(var_in), unique(var_out),
        #     alternative = c("two.sided"),
        #     mu = 0, # Hypothesis: Samples do not differ
        #     paired = FALSE,
        #     all.perms = TRUE, # Tries to get a exact p value
        #     num.sim = 1000000,
        #     plot = FALSE, # does not plot
        #     stat = mean
        # ) # compaires the means of the distributions
        # pvalue <- perm$p.value

        # return the calculated pvalues ----------------------------------------
        return(pvalue)
    }
}

#' get the list with all the features in the gene
#' @param selected_gene gene to get the feartures for
#' @param domains contains "seedID", "orthoID", "feature", "start",   "end"
#' @return dataframe for the specific gene containing "seedID",  "orthoID",
#' "feature", "start",   "end"
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_features <- function(selected_gene, domains){
    subset_domains <- {
        subset(
            domains,
            substr(
                domains$seedID,
                1,
                nchar(as.character(selected_gene))
            ) == selected_gene
        )
    }
    subset_domains <- subset_domains[!duplicated(subset_domains), ]
    return(subset_domains)
}

#' Get the database for each feature in a specific gene
#' @param data contains "seedID", "orthoID", "feature", "start", "end"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_prefix_features <- function(data){
    features <- data$feature
    choices <- gsub("_.*", "", features)
    choices <- choices[!duplicated(choices)]
    return(choices)
}

#' print list of available taxa
#' @param rank_select_gc rank selected for group compariosn
#' @param input_taxonID contains "seedID",  "orthoID", "feature", "start", "end"
#' @return avilable taxa containing "ncbiID", "fullName", "rank", "parentID"
#' @author Carla Mölbert (carla.moelbert@gmx.de)
taxa_select_gc <- function(rank_select_gc, input_taxonID){
    # if there is no rank set, there can not be any available taxa
    if (length(rank_select_gc) == 0) return()
    else{
        # load list of unsorted taxa
        if (is.null(input_taxonID)) dt <- get_taxonomy_matrix(
            FALSE, input_taxonID
        )
        else dt <- get_taxonomy_matrix(TRUE, input_taxonID)

        # load list of taxon name
        name_list <- phyloprofile::get_name_list()

        # get rank name from rank_select
        if (substr(rank_select_gc,3,3) == "_") {
            # rank_name <- substr(rank_select_gc, 4, nchar(rank_select_gc))
            rank_name <- rank_select_gc
        }
        else rank_name <- rank_select_gc

        choice <- as.data.frame
        choice <- rbind(dt[rank_name])
        colnames(choice) <- "ncbiID"
        choice <- merge(choice, name_list, by = "ncbiID", all = FALSE)
        return(choice)
    }
}

#' Create a list with all main taxanomy ranks
#' @export
#' @return A list of all main ranks (from strain to superkingdom)
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' get_taxonomy_ranks()

get_taxonomy_ranks <- function(){
    all_ranks <- list(
        "Strain " = "strain",
        "Species" = "species",
        "Genus" = "genus",
        "Family" = "family",
        "Order" = "order",
        "Class" = "class",
        "Phylum" = "phylum",
        "Kingdom" = "kingdom",
        "Superkingdom" = "superkingdom",
        "unselected" = ""
    )

    return(all_ranks)
}
