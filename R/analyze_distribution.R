#' Create data for percentage present taxa distribution
#' @param input_data dataframe contains raw input data in long format
#' @param rank_name name of the working taxonomy rank (e.g. "species", "family")
#' @return A dataframe ready for analysing the distribution of the percentage of
#' species in the selected supertaxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{main_long_raw}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' create_percentage_distribution_data(main_long_raw, "class")

create_percentage_distribution_data <- function(input_data,
                                                rank_name) {
    mdData <- input_data
    if (ncol(mdData) < 4) {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID")
    } else if (ncol(mdData) < 5) {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID", "var1")
    } else {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    }

    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"

    # (3) GET SORTED TAXONOMY LIST (3)
    input_taxonID <- phyloprofile::get_input_taxa_id(input_data)
    input_taxonName <- phyloprofile::get_input_taxa_name(
        rank_name, input_taxonID
    )
    ref_taxon <- input_taxonName$fullName[1]
    taxa_tree <- NULL

    taxa_list <- phyloprofile::sort_input_taxa(
        input_taxonID, input_taxonName, rank_name, ref_taxon, taxa_tree
    )

    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list, "supertaxon")

    # merge mdData, mdDatavar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
    if ("var1" %in% colnames(taxaMdData)) {
        taxaMdData$var1 <-
            suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    }
    if ("var2" %in% colnames(taxaMdData)) {
        taxaMdData$var2 <-
            suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    }
    # calculate % present species
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)

    finalPresSpecDt[!is.na(finalPresSpecDt$geneID),]
    return(finalPresSpecDt)
}

#' Create data for additional variable distribution
#' @usage create_variable_distribution_data(input_data, var1_cutoff_min,
#'     var1_cutoff_max, var2_cutoff_min, var2_cutoff_max)
#' @param input_data dataframe contains raw input data in long format
#' @param var1_cutoff_min min cutoff for var1
#' @param var1_cutoff_max max cutoff for var1
#' @param var2_cutoff_min min cutoff for var2
#' @param var2_cutoff_max max cutoff for var2
#' @return A dataframe ready for analysing the distribution of the additional
#' variable(s)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{main_long_raw}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' create_variable_distribution_data(
#'     main_long_raw, 0, 1, 0.5, 1
#' )

create_variable_distribution_data <- function(
    input_data,
    var1_cutoff_min,
    var1_cutoff_max,
    var2_cutoff_min,
    var2_cutoff_max
) {
    dataOrig <- input_data
    if (ncol(dataOrig) < 4) {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID")
        splitDt <- dataOrig[, c("orthoID")]
    } else if (ncol(dataOrig) < 5) {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID",
                                "var1")
        splitDt <- dataOrig[, c("orthoID",
                                "var1")]
    } else {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID",
                                "var1",
                                "var2")
        splitDt <- dataOrig[, c("orthoID", "var1", "var2")]
    }

    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[complete.cases(splitDt), ]

    if (length(levels(as.factor(splitDt$var2))) == 1) {
        if (levels(as.factor(splitDt$var2)) == "") {
            splitDt$var2 <- 0
        }
    }

    # convert factor into numeric for "var1" & "var2" column
    if ("var1" %in% colnames(splitDt)) {
        splitDt$var1 <- suppressWarnings(as.numeric(as.character(splitDt$var1)))
        # filter splitDt based on selected var1 cutoff
        splitDt <- splitDt[splitDt$var1 >= var1_cutoff_min
                            & splitDt$var1 <= var1_cutoff_max, ]
    }
    if ("var2" %in% colnames(splitDt)) {
        splitDt$var2 <- suppressWarnings(as.numeric(as.character(splitDt$var2)))
        # filter splitDt based on selected var2 cutoff
        splitDt <- splitDt[splitDt$var2 >= var2_cutoff_min
                            & splitDt$var2 <= var2_cutoff_max, ]
    }

    return(splitDt)
}

#' Create data for additional variable distribution (for a subset data)
#' @usage create_variable_distribution_data_subset(full_profile_data,
#'     distribution_data, selected_genes, selected_taxa)
#' @param full_profile_data dataframe contains the full processed profiles
#' @param distribution_data dataframe contains the full distribution data
#' @param selected_genes list of genes of interst
#' @param selected_taxa list of taxa of interest
#' @return A dataframe ready for analysing the distribution of the additional
#' variable(s) for a subset of genes and/or taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parse_info_profile}},
#' \code{\link{create_variable_distribution_data}},
#' \code{\link{full_processed_profile}}, \code{\link{main_long_raw}}
#' @examples
#' data("full_processed_profile", package="phyloprofile")
#' data("main_long_raw", package="phyloprofile")
#' distribution_data <- create_variable_distribution_data(
#'     main_long_raw, 0, 1, 0.5, 1
#' )
#' selected_genes <- "OG_1019"
#' selected_taxa <- c("Mammalia", "Echinoidea", "Gunneridae", "Mucorales",
#' "Alphaproteobacteria")
#' create_variable_distribution_data_subset(
#'     full_processed_profile,
#'     distribution_data,
#'     selected_genes,
#'     selected_taxa
#' )

create_variable_distribution_data_subset <- function(
    full_profile_data,
    distribution_data,
    selected_genes,
    selected_taxa
) {
    geneID <- NULL
    orthoID <- NULL
    var1.x <- NULL
    var2.y <- NULL
    supertaxonMod <- NULL

    allData <- full_profile_data
    splitDt <- distribution_data

    # get geneID and supertaxon name for splitDt
    splitDtName <- merge(
        splitDt, allData,
        by = "orthoID",
        all.x = TRUE
    )
    splitDtName$supertaxonMod <-
        substr(
            splitDtName$supertaxon,
            6,
            nchar(as.character(splitDtName$supertaxon))
        )
    splitDtName <- subset(
        splitDtName,
        select = c(
            orthoID,
            var1.x,
            var2.y,
            supertaxonMod,
            geneID
        )
    )
    colnames(splitDtName) <- c(
        "orthoID",
        "var1",
        "var2",
        "supertaxonMod",
        "geneID"
    )

    # filter
    if (selected_taxa[1] == "all" & selected_genes[1] != "all") {
        # select data from dataHeat for selected sequences only
        splitDt <- subset(splitDtName, geneID %in% selected_genes)
    } else if (selected_genes[1] == "all" & selected_taxa[1] != "all") {
        # select data from dataHeat for selected taxa only
        splitDt <- subset(splitDtName, supertaxonMod %in% selected_taxa)
    } else {
        # select data from dataHeat for selected sequences and taxa
        splitDt <- subset(
            splitDtName,
            geneID %in% selected_genes
            & supertaxonMod %in% selected_taxa
        )
    }

    return(splitDt)
}

#' Create distribution plot
#' @description Create distribution plot for one of the additional variable or
#' the percentage of the species present in the supertaxa.
#' @param data dataframe contains data for plotting
#' @param var_name name of the variable that need to be analyzed (either name of
#' variable 1 or variable 2 or "percentage of present taxa")
#' @param var_type type of variable (either "var1", "var2" or "presSpec")
#' @param percent range of percentage cutoff
#' @param dist_text_size text size of the distribution plot
#' @return A distribution plot as a ggplot object
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{main_long_raw}},
#' \code{\link{create_variable_distribution_data}},
#' \code{\link{create_variable_distribution_data_subset}},
#' \code{\link{create_percentage_distribution_data}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' data <- create_variable_distribution_data(
#'     main_long_raw, 0, 1, 0.5, 1
#' )
#' var_name <- "Variable abc"
#' var_type <- "var1"
#' percent <- c(0,1)
#' dist_text_size <- 12
#' create_var_dist_plot(
#'     data,
#'     var_name,
#'     var_type,
#'     percent,
#'     dist_text_size
#' )

create_var_dist_plot <- function(
    data,
    var_name,
    var_type,
    percent,
    dist_text_size
) {
    if (var_type == "presSpec") {
        # remove presSpec < cutoff_min or > cutoff_max
        if (percent[1] > 0) {
            data <- data[data$presSpec >= percent[1]
                        & data$presSpec <= percent[2], ]
        } else {
            if (percent[2] > 0) {
                data <- data[data$presSpec > 0 & data$presSpec <= percent[2], ]
            } else {
                data <- data[data$presSpec > 0, ]
            }
        }
    } else {
        data <- data[!is.na(data[,var_type]), ]
    }

    data.mean <- mean(data[,var_type])

    p <- ggplot(data, aes(x = data[,var_type])) +
        geom_histogram(binwidth = .01, alpha = .5, position = "identity") +
        geom_vline(
            data = data,
            aes(xintercept = data.mean, colour = "red"),
            linetype = "dashed",
            size = 1
        ) +
        theme_minimal()
    p <- p +
        theme(
            legend.position = "none",
            axis.title = element_text(size = dist_text_size),
            axis.text = element_text(size = dist_text_size)
        ) +
        labs(
            x = paste0(
                var_name,
                " (mean = ",
                round(mean(data[,var_type]), 3),
                ")"
            ),
            y = "Frequency"
        )

    return(p)
}
