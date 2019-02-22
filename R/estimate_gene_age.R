#' Calculate the phylogenetic gene age from profiles
#' @export
#' @usage estimate_gene_age(processed_profile_data, rank_name, ref_taxon,
#'     var1_cutoff, var2_cutoff, percent_cutoff)
#' @param processed_profile_data dataframe contains the full processed
#' phylogenetic profiles
#' @param rank_name taxonomy rank (e.g. "species", "genus", "family")
#' @param ref_taxon reference taxon (e.g. "Homo sapiens", "Homo", "Hominidae")
#' @param var1_cutoff cutoff for var1
#' @param var2_cutoff cutoff for var2
#' @param percent_cutoff cutoff for percentage of species present in each
#' supertaxon
#' @return A dataframe contains estimated gene ages
#' @importFrom data.table setDT
#' @importFrom data.table setnames
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{parse_info_profile}} for creating a full processed
#' profile dataframe; \code{\link{get_name_list}} and
#' \code{\link{get_taxonomy_matrix}} for getting taxonomy info,
#' \code{\link{full_processed_profile}} for a demo input dataframe
#' @examples
#' data("full_processed_profile", package="phyloprofile")
#' rank_name <- "class"
#' ref_taxon <- "Mammalia"
#' processed_profile_data <- full_processed_profile
#' var1_cutoff <- c(0, 1)
#' var2_cutoff <- c(0, 1)
#' percent_cutoff <- c(0, 1)
#' estimate_gene_age(
#'     processed_profile_data,
#'     rank_name,
#'     ref_taxon,
#'     var1_cutoff, var2_cutoff, percent_cutoff
#' )

estimate_gene_age <- function(
    processed_profile_data,
    rank_name, ref_taxon,
    var1_cutoff, var2_cutoff, percent_cutoff
){
    rankList <- c(
        "family", "class", "phylum", "kingdom", "superkingdom", "root"
    )

    # get selected (super)taxon ID
    taxa_list <- phyloprofile::get_name_list()
    superID <- {
        as.numeric(
            taxa_list[
                taxa_list$fullName == ref_taxon & taxa_list$rank == rank_name,
            ]$ncbiID
        )
    }

    # full non-duplicated taxonomy data
    Dt <- get_taxonomy_matrix(FALSE, NULL)

    # subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[, c("abbrName", rankList)]

    # get (super)taxa IDs for one of representative species
    # get all taxon info for 1 representative
    first_line <- Dt[Dt[, rank_name] == superID, ][1, ]
    sub_first_line <- first_line[, c("abbrName", rankList)]

    # compare each taxon ncbi IDs with selected taxon
    # and create a "category" data frame
    catDf <- data.frame("ncbiID" = character(),
                        "cat" = character(),
                        stringsAsFactors = FALSE)
    for (i in seq_len(nrow(subDt))) {
        cat <- subDt[i, ] %in% sub_first_line
        cat[cat == FALSE] <- 0
        cat[cat == TRUE] <- 1
        cat <- paste0(cat, collapse = "")
        catDf[i, ] <- c(as.character(subDt[i, ]$abbrName), cat)
    }

    # get main input data
    mdData <- droplevels(processed_profile_data)
    mdData <- mdData[, c(
        "geneID", "ncbiID", "orthoID", "var1", "var2", "presSpec"
    )]

    ### add "category" into mdData
    mdDataExtended <- merge(mdData,
                            catDf,
                            by = "ncbiID",
                            all.x = TRUE)

    mdDataExtended$var1[mdDataExtended$var1 == "NA"
                        | is.na(mdDataExtended$var1)] <- 0
    mdDataExtended$var2[mdDataExtended$var2 == "NA"
                        | is.na(mdDataExtended$var2)] <- 0

    # remove cat for "NA" orthologs
    # and also for orthologs that do not fit cutoffs
    if (nrow(mdDataExtended[mdDataExtended$orthoID == "NA"
        | is.na(mdDataExtended$orthoID), ]) > 0) {
        mdDataExtended[mdDataExtended$orthoID == "NA"
        | is.na(mdDataExtended$orthoID), ]$cat <- NA
    }

    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

    # filter by %specpres, var1, var2 ..
    mdDataExtended$cat[mdDataExtended$var1 < var1_cutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var1 > var1_cutoff[2]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 < var2_cutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 > var2_cutoff[2]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec < percent_cutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec > percent_cutoff[2]] <- NA

    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

    ### get the furthest common taxon with selected taxon for each gene
    gene_ageDf <- as.data.frame(
        tapply(mdDataExtended$cat, mdDataExtended$geneID, min)
    )

    setDT(gene_ageDf, keep.rownames = TRUE)[]
    setnames(gene_ageDf, seq_len(2), c("geneID", "cat"))  # rename columns
    row.names(gene_ageDf) <- NULL   # remove row names

    ### convert cat into gene_age
    gene_ageDf$age[gene_ageDf$cat == "0000001"] <- "07_LUCA"
    gene_ageDf$age[gene_ageDf$cat == "0000011" | gene_ageDf$cat == "0000010"] <-
        paste0(
            "06_",
            as.character(
                taxa_list$fullName[
                    taxa_list$ncbiID == sub_first_line$superkingdom
                    & taxa_list$rank == "superkingdom"
                ]
            )
        )

    gene_ageDf$age[gene_ageDf$cat == "0000111"] <-
        paste0(
            "05_",
            as.character(
                taxa_list$fullName[
                    taxa_list$ncbiID == sub_first_line$kingdom
                    & taxa_list$rank == "kingdom"
                ]
            )
        )

    gene_ageDf$age[gene_ageDf$cat == "0001111"] <-
        paste0(
            "04_",
            as.character(
                taxa_list$fullName[
                    taxa_list$ncbiID == sub_first_line$phylum
                    & taxa_list$rank == "phylum"
                ]
            )
        )

    gene_ageDf$age[gene_ageDf$cat == "0011111"] <-
        paste0(
            "03_",
            as.character(
                taxa_list$fullName[
                    taxa_list$ncbiID == sub_first_line$class
                    & taxa_list$rank == "class"
                ]
            )
        )

    gene_ageDf$age[gene_ageDf$cat == "0111111"] <-
        paste0(
            "02_",
            as.character(
                taxa_list$fullName[
                    taxa_list$ncbiID == sub_first_line$family
                    & taxa_list$rank == "family"
                ]
            )
        )

    gene_ageDf$age[gene_ageDf$cat == "1111111"] <-
        paste0(
            "01_",
            as.character(
                taxa_list$fullName[
                    taxa_list$fullName == ref_taxon
                    & taxa_list$rank == rank_name
                ]
            )
        )

    # return gene_age data frame
    gene_ageDf <- gene_ageDf[, c("geneID", "cat", "age")]
    gene_ageDf$age[is.na(gene_ageDf$age)] <- "Undef"

    return(gene_ageDf)
}

#' Create data for plotting gene ages
#' @param gene_ageDf data of estimated gene age
#' @return A dataframe for plotting gene age plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{estimate_gene_age}}
#' @examples
#' gene_ageDf <- data.frame(
#' geneID = c("OG_1017", "OG_1019"),
#' cat = c("0000001", "0000001"),
#' age = c("07_LUCA", "07_LUCA")
#' )
#' gene_age_plotDf(gene_ageDf)

gene_age_plotDf <- function(gene_ageDf){
    countDf <- plyr::count(gene_ageDf, c("age"))
    countDf$percentage <- round(countDf$freq / sum(countDf$freq) * 100)
    countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
    return(countDf)
}

#' Create gene age plot
#' @param count_df data for plotting gene age
#' @param gene_age_text text size
#' @return A gene age distribution plot as a ggplot2 object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @export
#' @seealso \code{\link{estimate_gene_age}} and \code{\link{gene_age_plotDf}}
#' @examples
#' count_df <- data.frame(
#'     age = "OG_1017",
#'     frea = 2,
#'     percentage = 100,
#'     pos = 50
#' )
#' gene_age_text <- 1
#' create_gene_age_plot(count_df, gene_age_text)

create_gene_age_plot <- function(count_df, gene_age_text){
    age <- NULL
    percentage <- NULL
    pos <- NULL
    freq <- NULL

    p <- ggplot(count_df, aes(fill = age, y = percentage, x = 1)) +
        geom_bar(stat = "identity") +
        scale_y_reverse() +
        coord_flip() +
        theme_minimal()
    p <- p + geom_text(
        data = count_df,
        aes(x = 1, y = 100 - pos, label = paste0(freq, "\n", percentage, "%")),
        size = 4 * gene_age_text
    )
    p <- p + theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 12 * gene_age_text),
            axis.title = element_blank(), axis.text = element_blank()
        ) +
        scale_fill_brewer(palette = "Spectral") +
        guides(fill = guide_legend(
            nrow = max(round(nrow(count_df) / 3, 0), 1),
            byrow = TRUE
        ))
    return(p)
}
