#' Identify core genes for a list of selected taxa
#' @export
#' @usage get_core_gene(rank_name, taxa_core, processed_profile_data,
#'     var1_cutoff, var2_cutoff, percent_cutoff, core_coverage)
#' @param rank_name taxonomy rank (e.g. "species", "genus", "family")
#' @param taxa_core name list of selected taxa
#' @param processed_profile_data dataframe contains the full processed
#' phylogenetic profiles
#' @param var1_cutoff cutoff for var1
#' @param var2_cutoff cutoff for var2
#' @param percent_cutoff cutoff for percentage of species present in each
#' supertaxon
#' @param core_coverage the least number of selected taxa should be considered
#' @return A list of core genes
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{parse_info_profile}} for creating a full processed
#' profile dataframe
#' @examples
#' data("full_processed_profile", package="phyloprofile")
#' rank_name <- "class"
#' taxa_core <- c("Mammalia", "Mucorales", "Alphaproteobacteria")
#' processed_profile_data <- full_processed_profile
#' var1_cutoff <- c(0.75, 1.0)
#' var2_cutoff <- c(0.75, 1.0)
#' percent_cutoff <- c(0.0, 1.0)
#' core_coverage <- 1
#' get_core_gene(
#'     rank_name,
#'     taxa_core,
#'     processed_profile_data,
#'     var1_cutoff, var2_cutoff,
#'     percent_cutoff, core_coverage
#' )

get_core_gene <- function(
    rank_name,
    taxa_core,
    processed_profile_data,
    var1_cutoff, var2_cutoff,
    percent_cutoff, core_coverage
) {
    supertaxonID <- NULL
    mVar1 <- NULL
    mVar2 <- NULL
    presSpec <- NULL
    Freq <- NULL

    # get ID list of chosen taxa
    taxa_list <- phyloprofile::get_name_list()

    if ("none" %in% taxa_core) {
        superID <- NA
    } else {
        superID <- taxa_list$ncbiID[
            taxa_list$fullName %in% taxa_core
            & taxa_list$rank %in% c(rank_name, "norank")
        ]
    }

    # get main input data
    mdData <- processed_profile_data
    if (is.null(mdData)) return()
    mdData <- mdData[, c(
        "geneID",
        "ncbiID",
        "fullName",
        "supertaxon",
        "supertaxonID",
        "rank",
        "presSpec",
        "mVar1",
        "mVar2"
    )]

    # filter by var1 and var2 cutoffs
    var1_cutoff_min <- var1_cutoff[1]
    var1_cutoff_max <- var1_cutoff[2]
    var2_cutoff_min <- var2_cutoff[1]
    var2_cutoff_max <- var2_cutoff[2]

    if (!is.null(var1_cutoff_max)) {
        if (!is.na(var1_cutoff_max)) {
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar1 >= var1_cutoff_min
            )
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar1 <= var1_cutoff_max
            )
        }
    }

    if (!is.null(var2_cutoff_max)) {
        if (!is.na(var2_cutoff_max)) {
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar2 >= var2_cutoff_min
            )
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar2 <= var2_cutoff_max
            )
        }
    }

    # filter by selecting taxa
    if (is.na(superID[1])) return(NULL) #mdData <- NULL
    else {
        data <- subset(
            mdData, supertaxonID %in% superID
            & presSpec >= percent_cutoff[1]
        )
        data <- subset(
            data, supertaxonID %in% superID
            & presSpec <= percent_cutoff[2]
        )

        # get supertaxa present in each geneID
        supertaxonCount <- as.data.frame(
            plyr::count(data, c("geneID", "supertaxonID"))
        )

        # count number of supertaxa present in each geneID
        # and get min number of supertaxa muss be taken into account
        count <- as.data.frame(table(supertaxonCount$geneID))
        require_coverage <- length(superID) * (core_coverage / 100)

        # get only gene that contains orthologs in that coverage # of taxa
        core_gene <- subset(count, Freq >= require_coverage)
        core_gene$Var1 <- factor(core_gene$Var1)

        return(levels(core_gene$Var1))
    }
}
