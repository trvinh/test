# Functions for parsing and pre-processing phyloprofile input

#' Get list of pre-installed NCBI taxon names
#' @description Get all NCBI taxon names from "data/taxonNamesReduced.txt"
#' @export
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' get_name_list()

get_name_list <- function() {
    nameReduced_file <- paste(
        system.file(package="phyloprofile"),
        "phyloprofile/data/taxonNamesReduced.txt",
        sep="/"
    )

    if (!file.exists(nameReduced_file)) {
        fileURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "taxonNamesReduced.txt"
        )
        res <- tryCatch(
            utils::download.file(
                fileURL,
                destfile = nameReduced_file,
                method="auto"
            ),
            error=function(e) 1
        )
    }

    name_list <- as.data.frame(read.table(
        nameReduced_file,
        sep = "\t",
        header = TRUE,
        fill = TRUE
    ))
    name_list$fullName <- as.character(name_list$fullName)
    name_list$rank <- as.character(name_list$rank)
    name_list <- name_list[!duplicated(name_list), ]

    return(name_list)
}

#' Get taxonomy matrix
#' @description Get full taxonomy matrix from "data/taxonomyMatrix.txt" or
#' only a subset of matrix based on an input taxon list
#' @export
#' @param subset_taxa_check subset taxonomy matrix based on input taxon IDs
#' (TRUE/FALSE)
#' @param taxon_IDs list of input taxon IDs
#' @return Data frame contains the (subset of) taxonomy matrix
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' # get full pre-installed taxonomy matrix
#' get_taxonomy_matrix(FALSE, NULL)
#' # get taxonomy matrix for a list of taxon IDs
#' taxon_ids <- c("ncbi10020", "ncbi10090")
#' get_taxonomy_matrix(TRUE, taxon_ids)

get_taxonomy_matrix <- function(subset_taxa_check, taxon_IDs){
    taxonomyMatrix_file <- paste(
        system.file(package="phyloprofile"),
        "phyloprofile/data/taxonomyMatrix.txt",
        sep="/"
    )

    file.exists(taxonomyMatrix_file)
    if (!file.exists(taxonomyMatrix_file)) {
        fileURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "taxonomyMatrix.txt"
        )
        res <- tryCatch(
            utils::download.file(
                fileURL,
                destfile = taxonomyMatrix_file,
                method="auto"
            ),
            error=function(e) 1
        )
    }

    dt <- as.data.frame(read.table(
        taxonomyMatrix_file,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = TRUE
    ))
    if (subset_taxa_check) {
        dt <- dt[dt$abbrName  %in% taxon_IDs, ]
    }
    return(dt)
}

#' Get ID list of input taxa from the main input
#' @param raw_profile long dataframe of input profile
#' @return List of all input taxon IDs
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{create_long_matrix}}, \code{\link{main_long_raw}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' get_input_taxa_id(main_long_raw)

get_input_taxa_id <- function(raw_profile){
    if (is.null(raw_profile)) return()
    inputTaxa <- levels(as.factor(raw_profile$ncbiID))

    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    if (inputTaxa[1] == "geneID") {
        # remove "geneID" element from vector inputTaxa
        inputTaxa <- inputTaxa[-1]
    }
    # return input taxa
    return(inputTaxa)
}

#' Get NCBI taxon names for a selected list of taxa
#' @description Get NCBI taxon names from "data/taxonNamesReduced.txt" for
#' a selected list of taxon
#' @param rank_name taxonomy rank (e.g. "species","phylum",...)
#' @param taxon_IDs list of taxon IDs (check get_input_taxa_id())
#' @return List of full names, taxonomy ranks and parent IDs for the input taxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{get_input_taxa_id}} for getting input taxon IDs,
#' \code{\link{get_name_list}} for getting the full taxon name list
#' @examples
#' taxon_ids <- c("ncbi10020", "ncbi10090")
#' get_input_taxa_name("species", taxon_ids)

get_input_taxa_name <- function(rank_name, taxon_IDs){
    # load list of unsorted taxa
    Dt <- get_taxonomy_matrix(TRUE, taxon_IDs)

    # load list of taxon name
    nameList <- get_name_list()

    choice <- as.data.frame
    choice <- rbind(Dt[rank_name])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,
                    nameList,
                    by = "ncbiID",
                    all = FALSE)
    return(choice)
}

#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @param taxon_IDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...)
#' @param taxon_names list of taxon names and their corresponding taxonomy
#' ranks, IDs,...
#' @param rank_name working taxonomy rank (e.g. "species", "phylum",...)
#' @param ref_taxon selected reference taxon
#' @param taxa_tree input taxonomy tree (optional)
#' @return Taxonomy matrix for the input taxa ordered by the selected reference
#' taxon. This matrix is sorted either based on the NCBI taxonomy info, or
#' based on an user-defined taxonomy tree (if provided).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{get_name_list}}, \code{\link{get_taxonomy_matrix}},
#' \code{\link{create_rooted_tree}}, \code{\link{sort_taxa_from_tree}},
#' \code{\link{get_input_taxa_name}}, \code{\link{get_input_taxa_id}},
#' \code{\link{create_long_matrix}}
#' @examples
#' taxon_IDs <- c(
#'     "ncbi272557", "ncbi176299", "ncbi3702", "ncbi876142", "ncbi9606"
#' )
#' taxon_names <- get_input_taxa_name("species", taxon_IDs)
#' sort_input_taxa(taxon_IDs, taxon_names, "species", "Homo sapiens", NULL)

sort_input_taxa <- function(taxon_IDs,
                            taxon_names,
                            rank_name,
                            ref_taxon,
                            taxa_tree){
    ncbiID <- NULL
    fullName <- NULL
    abbrName <- NULL

    fullname_list <- get_name_list()

    # get selected supertaxon ID(s)
    rankNameTMP <- taxon_names$rank[taxon_names$fullName == ref_taxon]
    if (rank_name == "strain") {
        superID <-
            as.numeric(
                fullname_list$ncbiID[
                    fullname_list$fullName == ref_taxon
                    & fullname_list$rank == "norank"
                ]
            )
    } else {
        superID <-
            as.numeric(
                fullname_list$ncbiID[
                    fullname_list$fullName == ref_taxon
                    & fullname_list$rank == rankNameTMP[1]
                ]
            )
    }

    # get full taxonomy data
    Dt <- get_taxonomy_matrix(FALSE, NULL)

    # representative taxon
    repTaxon <- Dt[Dt[, rank_name] == superID, ][1, ]

    # THEN, SORT TAXON LIST BASED ON TAXONOMY TREE
    if (is.null(taxa_tree)) {
        # prepare Df for calculating distance matrix
        distDf <- subset(Dt, select = -c(ncbiID, fullName))
        row.names(distDf) <- distDf$abbrName
        distDf <- distDf[, -1]
        # create taxonomy tree rooted by ref_taxon
        taxa_tree <- create_rooted_tree(distDf, as.character(repTaxon$abbrName))
    } else {
        taxa_tree <- ape::root(
            taxa_tree,
            outgroup = as.character(repTaxon$abbrName),
            resolve.root = TRUE
        )
    }

    taxonList <- sort_taxa_from_tree(taxa_tree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]

    # subset to get list of input taxa only
    sortedDt <- subset(sortedDt, abbrName %in% taxon_IDs)

    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(
        sortedDt,
        select = c(
            "abbrName", "ncbiID", "fullName", as.character(rank_name)
        )
    )
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")

    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullname_list,
        by = "ncbiID",
        all.x = TRUE,
        sort = FALSE
    )
    sortedOut$species <- as.character(sortedOut$species)

    # add order_prefix to supertaxon name
    # and add prefix "ncbi" to taxon_ncbiID (column "species")
    prefix <- 1001

    ## create new column for sorted supertaxon
    sortedOut$sortedSupertaxon <- 0
    sortedOut$sortedSupertaxon[1] <- paste0(prefix,
                                            "_",
                                            sortedOut$fullName.y[1])
    sortedOut$species[1] <- paste0("ncbi", sortedOut$species[1])

    if (nrow(sortedOut) > 1) {
        for (i in 2:nrow(sortedOut)) {
            ## increase prefix if changing to another supertaxon
            if (sortedOut$fullName.y[i] != sortedOut$fullName.y[i - 1]) {
                prefix <- prefix + 1
            }
            sortedOut$sortedSupertaxon[i] <- paste0(prefix,
                                                    "_",
                                                    sortedOut$fullName.y[i])
            sortedOut$species[i] <- paste0("ncbi", sortedOut$species[i])
        }
    }

    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c(
        "abbrName",
        "taxonID",
        "fullName.x",
        "species",
        "ncbiID",
        "sortedSupertaxon",
        "rank",
        "category"
    )]
    colnames(sortedOut) <- c(
        "abbrName",
        "taxonID",
        "fullName",
        "ncbiID",
        "supertaxonID",
        "supertaxon",
        "rank",
        "category"
    )

    sortedOut$taxonID <- as.numeric(sortedOut$taxonID)
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)

    return(sortedOut)
}

#' Calculate percentage of present species in each ortholog group
#' @export
#' @param profile_with_tax long data frame of main phyloprofile input together
#' with their taxonomy info
#' @param taxa_count number of species occur in each supertaxon
#' @return A data frame with % of present species for each seed protein in
#' each supertaxon
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{profile_with_taxonomy}} for a demo input data
#' @examples
#' # NOTE: for internal testing only - not recommended for outside using
#' data("profile_with_taxonomy", package="phyloprofile")
#' taxa_count <- plyr::count(profile_with_taxonomy, "supertaxon")
#' taxa_count$freq <- 1
#' calc_pres_spec(profile_with_taxonomy, taxa_count)

calc_pres_spec <- function(profile_with_tax, taxa_count){
    paralog <- NULL
    abbrName <- NULL
    profile_with_tax <- profile_with_tax[profile_with_tax$orthoID != "NA", ]

    # get geneID and supertaxon
    gene_id_supertaxon <- subset(
        profile_with_tax,
        select = c("geneID", "supertaxon", "paralog", "abbrName")
    )
    # remove duplicated rows
    gene_id_supertaxon <- gene_id_supertaxon[!duplicated(gene_id_supertaxon), ]

    # remove NA rows from profile_with_tax
    profile_with_tax_no_na <-
        profile_with_tax[profile_with_tax$orthoID != "NA", ]

    # count present frequency of supertaxon for each gene
    gene_supertaxon_count <- plyr::count(
        profile_with_tax_no_na,
        c("geneID", "supertaxon")
    )

    # merge with taxa_count to get total number of species of each supertaxon
    # and calculate presSpec
    pres_spec_dt <- merge(
        gene_supertaxon_count, taxa_count, by = "supertaxon", all.x = TRUE
    )

    spec_count <- plyr::count(gene_id_supertaxon, c("geneID", "supertaxon"))
    pres_spec_dt <- merge(
        pres_spec_dt, spec_count, by = c("geneID", "supertaxon")
    )

    pres_spec_dt$presSpec <- pres_spec_dt$freq / pres_spec_dt$freq.y

    pres_spec_dt <- pres_spec_dt[pres_spec_dt$presSpec <= 1, ]
    pres_spec_dt <- pres_spec_dt[order(pres_spec_dt$geneID), ]
    pres_spec_dt <- pres_spec_dt[, c("geneID", "supertaxon", "presSpec")]

    # add absent supertaxon into pres_spec_dt
    gene_id_supertaxon <- subset(
        gene_id_supertaxon, select = -c(paralog, abbrName)
    )
    final_pres_spec_dt <- merge(pres_spec_dt,
                                gene_id_supertaxon,
                                by = c("geneID", "supertaxon"),
                                all.y = TRUE)
    final_pres_spec_dt$presSpec[is.na(final_pres_spec_dt$presSpec)] <- 0

    # remove duplicated rows
    final_pres_spec_dt <- final_pres_spec_dt[!duplicated(final_pres_spec_dt), ]

    # return final_pres_spec_dt
    return(final_pres_spec_dt)
}

#' Parsing info for phylogenetic profiles
#' @description Creating main dataframe for the input phylogenetic profiles with
#' the selected input taxonomy level (e.g. strain, species) and reference taxon.
#' The output contains the number of paralogs, percentage of species presence
#' in each supertaxon, and the max/min/mean/median of VAR1 and VAR2.
#' @usage parse_info_profile( input_df, sorted_input_taxa, var1_aggregate_by, 
#'     var2_aggregate_by)
#' @param input_df input profiles in long format
#' @param sorted_input_taxa sorted taxonomy data for the input taxa
#' (check sort_input_taxa())
#' @param var1_aggregate_by aggregate method for VAR1 (min, max, mean or median)
#' @param var2_aggregate_by aggregate method for VAR2 (min, max, mean or median)
#' @importFrom stats aggregate
#' @return Dataframe contains all info for input profiles (a full processed
#' profile)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{create_long_matrix}}, \code{\link{sort_input_taxa}},
#' \code{\link{calc_pres_spec}}, \code{\link{main_long_raw}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' input_df <- main_long_raw
#' taxon_IDs <- get_input_taxa_id(input_df)
#' taxon_names <- get_input_taxa_name("class", taxon_IDs)
#' sorted_input_taxa <- sort_input_taxa(
#'     taxon_IDs, taxon_names, "class", "Mammalia", NULL
#' )
#' var1_aggregate_by <- "max"
#' var2_aggregate_by <- "mean"
#' parse_info_profile(
#'     input_df, sorted_input_taxa, var1_aggregate_by, var2_aggregate_by
#' )

parse_info_profile <- function(
    input_df,
    sorted_input_taxa,
    var1_aggregate_by,
    var2_aggregate_by
) {
    mdData <- input_df
    # rename columns of 2 additional variables
    if (ncol(mdData) > 3) {
        if (ncol(mdData) < 5) {
            colnames(mdData)[4] <- "var1"
        } else {
            colnames(mdData)[c(4,5)] <- c("var1", "var2")
        }
    }

    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"

    # (1) GET SORTED TAXONOMY LIST
    taxa_list <- sorted_input_taxa

    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list, "supertaxon")

    # merge mdData, mdDataVar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
    taxaMdData$var1 <-
        suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    taxaMdData$var2 <-
        suppressWarnings(as.numeric(as.character(taxaMdData$var2)))

    # (2) calculate PERCENTAGE of PRESENT SPECIES
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)

    # (3) calculate max/min/mean/median VAR1 for every supertaxon of each gene
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1), ]
    # calculate m var1
    mVar1Dt <- stats::aggregate(
        taxaMdDataNoNA[, "var1"],
        list(taxaMdDataNoNA$supertaxon, taxaMdDataNoNA$geneID),
        FUN = var1_aggregate_by
    )
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")

    # (4) calculate max/min/mean/median VAR2 for each super taxon
    # remove NA rows from taxaMdData
    taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2), ]
    # calculate max/min/mean/median VAR2
    if (nrow(taxaMdDataNoNA_var2) > 0) {
        mVar2Dt <- aggregate(
            taxaMdDataNoNA_var2[, "var2"],
            list(taxaMdDataNoNA_var2$supertaxon, taxaMdDataNoNA_var2$geneID),
            FUN = var2_aggregate_by
        )
        colnames(mVar2Dt) <- c("supertaxon", "geneID", "mVar2")
    } else {
        mVar2Dt <- taxaMdData[, c("supertaxon", "geneID")]
        mVar2Dt$mVar2 <- 0
    }

    # (3+4) & join mVar2 together with mVar1 scores into one df
    scoreDf <- merge(
        mVar1Dt, mVar2Dt, by = c("supertaxon", "geneID"), all = TRUE
    )

    # (2+3+4) add presSpec and mVar1 into taxaMdData
    presMdData <- merge(taxaMdData,
                        finalPresSpecDt,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(presMdData,
                        scoreDf,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(fullMdData,
                        taxaCount, by = ("supertaxon"),
                        all.x = TRUE)
    # rename "freq" into "numberSpec"
    names(fullMdData)[names(fullMdData) == "freq"] <- "numberSpec"

    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData) == "orthoID.x"] <- "orthoID"

    # parsed input data frame and return
    fullMdData <- fullMdData[!duplicated(fullMdData), ]
    return(fullMdData)
}

#' Reduce the full processed profile into supertaxon level
#' @description Reduce data of the phylogenetic profiles from input taxonomy
#' rank into supertaxon level (e.g. from species to phylum)
#' @param full_profile dataframe contains the full processed profiles
#' @return A reduced dataframe contains only profiles for the selected rank.
#' This dataframe contains only supertaxa and their value (%present, mVar1 &
#' mVar2) for each gene.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parse_info_profile}} for creating a full processed
#' profile dataframe, \code{\link{full_processed_profile}} for a demo full
#' processed profile dataframe
#' @examples
#' data("full_processed_profile", package="phyloprofile")
#' reduce_profile(full_processed_profile)

reduce_profile <- function(full_profile) {
    fullMdData <- full_profile

    # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1) {
        if (unique(levels(as.factor(fullMdData$numberSpec))) == 1) {
            superDfExt <- fullMdData[, c(
                "geneID",
                "supertaxon",
                "supertaxonID",
                "var1",
                "presSpec",
                "category",
                "orthoID",
                "var2",
                "paralog"
            )]
            flag <- 0
        }
    }

    if (flag == 1) {
        # get representative orthoID that has m VAR1 for each supertaxon
        mOrthoID <- fullMdData[, c(
            "geneID",
            "supertaxon",
            "var1",
            "mVar1",
            "orthoID"
        )]
        mOrthoID <- subset(mOrthoID, mOrthoID$var1 == mOrthoID$mVar1)
        colnames(mOrthoID) <- c(
            "geneID",
            "supertaxon",
            "var1",
            "mVar1",
            "orthoID"
        )
        mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
        mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID")]
        mOrthoID <- mOrthoID[!duplicated(mOrthoID[, seq_len(2)]), ]

        # get data set for phyloprofile plotting (contains only supertaxa info)
        superDf <- subset(fullMdData, select = c(
            "geneID",
            "supertaxon",
            "supertaxonID",
            "mVar1",
            "presSpec",
            "category",
            "mVar2",
            "paralog"
        ))
        superDf$paralog <- 1
        superDf <- superDf[!duplicated(superDf), ]

        superDfExt <- merge(superDf, mOrthoID, by = c("geneID", "supertaxon"),
                            all.x = TRUE)
        superDfExt <- superDfExt[, c(
            "geneID",
            "supertaxon",
            "supertaxonID",
            "mVar1",
            "presSpec",
            "category",
            "orthoID",
            "mVar2",
            "paralog"
        )]

        # rename mVar to var
        names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
        names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }

    return(superDfExt)
}

#' Create data for plotting the phylogentic profiles
#' @usage create_profile_data(superTaxon_data, ref_taxon, percent_cutoff,
#'     coortholog_cutoff_max, var1_cutoff, var2_cutoff, var1_relation,
#'     var2_relation, group_by_cat, cat_dt)
#' @param superTaxon_data a reduced dataframe contains info for all profiles in
#' the selected taxonomy rank.
#' @param ref_taxon selected reference taxon
#' @param percent_cutoff min and max cutoffs for percentage of species present
#' in a supertaxon
#' @param coortholog_cutoff_max maximum number of co-orthologs allowed
#' @param var1_cutoff min and max cutoffs for var1
#' @param var2_cutoff min anc max cutoffs for var2
#' @param var1_relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param var2_relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param group_by_cat group genes by their categories (TRUE or FALSE)
#' @param cat_dt dataframe contains gene categories
#' (optional, NULL if group_by_cat = FALSE or no info provided)
#' @return A dataframe ready for generating profile plot.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parse_info_profile}} and \code{\link{reduce_profile}}
#' for generating input dataframe, \code{\link{full_processed_profile}} for a
#' demo full processed profile dataframe
#' @examples
#' data("full_processed_profile", package="phyloprofile")
#' superTaxon_data <- reduce_profile(full_processed_profile)
#' ref_taxon <- "Mammalia"
#' percent_cutoff <- c(0.0, 1.0)
#' coortholog_cutoff_max <- 10
#' var1_cutoff <- c(0.75, 1.0)
#' var2_cutoff <- c(0.5, 1.0)
#' var1_relation <- "protein"
#' var2_relation <- "species"
#' group_by_cat <- FALSE
#' cat_dt <- NULL
#' create_profile_data(
#'     superTaxon_data,
#'     ref_taxon,
#'     percent_cutoff,
#'     coortholog_cutoff_max,
#'     var1_cutoff,
#'     var2_cutoff,
#'     var1_relation,
#'     var2_relation,
#'     group_by_cat,
#'     cat_dt
#' )

create_profile_data <- function(
    superTaxon_data,
    ref_taxon,
    percent_cutoff,
    coortholog_cutoff_max,
    var1_cutoff,
    var2_cutoff,
    var1_relation,
    var2_relation,
    group_by_cat,
    cat_dt
) {
    dataHeat <- superTaxon_data
    
    # cutoffs
    percent_cutoff_min <- percent_cutoff[1]
    percent_cutoff_max <- percent_cutoff[2]
    var1_cutoff_min <- var1_cutoff[1]
    var1_cutoff_max <- var1_cutoff[2]
    var2_cutoff_min <- var2_cutoff[1]
    var2_cutoff_max <- var2_cutoff[2]

    # get selected supertaxon name
    in_select <- ref_taxon

    ### replace insufficient values according to the thresholds by NA or 0
    # based on presSpec or # of co-orthologs
    number_coortholog <- levels(as.factor(dataHeat$paralog))
    if (length(number_coortholog) > 1) {
        dataHeat$presSpec[
            dataHeat$supertaxon != in_select
            & dataHeat$paralog > coortholog_cutoff_max
        ] <- 0
        if (var2_relation == "protein") {
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$paralog > coortholog_cutoff_max
            ] <- NA
        }
    } else {
        if (length(levels(as.factor(dataHeat$presSpec))) > 1)
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$presSpec < percent_cutoff_min
            ] <- 0
        dataHeat$presSpec[
            dataHeat$supertaxon != in_select
            & dataHeat$presSpec > percent_cutoff_max
        ] <- 0
        if (var2_relation == "protein") {
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$presSpec < percent_cutoff_min
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$presSpec > percent_cutoff_max
            ] <- NA
        }
    }

    dataHeat$presSpec[
        dataHeat$supertaxon != in_select
        & dataHeat$var1 < var1_cutoff_min
    ] <- 0
    dataHeat$presSpec[
        dataHeat$supertaxon != in_select
        & dataHeat$var1 > var1_cutoff_max
    ] <- 0
    if (var1_relation == "protein") {
        if (var2_relation == "protein") {
            # prot-prot: remove complete cell if one variable not sufficient
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 < var2_cutoff_min
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 > var2_cutoff_max
            ] <- 0
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 < var1_cutoff_min
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 > var1_cutoff_max
            ] <- NA
            dataHeat$var1[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 < var2_cutoff_min
            ] <- NA
            dataHeat$var1[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 > var2_cutoff_max
            ] <- NA
        } else {
            # prot-spec: var1 depend on var2
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 < var2_cutoff_min
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var2 > var2_cutoff_max
            ] <- 0
        }
    } else {
        if (var2_relation == "species") {
            # spec-spec: remove var1 and var2 independently
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 < var1_cutoff_min
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 > var1_cutoff_max
            ] <- 0
        } else {
            # spec-prot: var2 depend on var1
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 < var1_cutoff_min
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != in_select
                & dataHeat$var1 > var1_cutoff_max
            ] <- NA
        }
    }

    dataHeat$var1[
        dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min
    ] <- NA
    dataHeat$var1[
        dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max
    ] <- NA
    dataHeat$var2[
        dataHeat$supertaxon != in_select & dataHeat$var2 < var2_cutoff_min
    ] <- NA
    dataHeat$var2[
        dataHeat$supertaxon != in_select & dataHeat$var2 > var2_cutoff_max
    ] <- NA

    dataHeat <- droplevels(dataHeat)  # delete unused levels
    dataHeat$geneID <- as.factor(dataHeat$geneID)
    dataHeat$supertaxon <- as.factor(dataHeat$supertaxon)

    ### add gene categories (if provided)
    if (group_by_cat == TRUE) {
        if (is.null(cat_dt)) {
            cat_dt <- data.frame( geneID = levels(dataHeat$geneID))
            cat_dt$group <- "no_category"
        }

        # create a dataframe that contain all genes and all taxa
        dataHeat_cat <- data.frame(
            supertaxon = rep(
                levels(dataHeat$supertaxon), nlevels(dataHeat$geneID)
            ),
            geneID = rep(
                levels(dataHeat$geneID), each = nlevels(dataHeat$supertaxon)
            )
        )

        dataHeat_cat <- merge(dataHeat_cat, cat_dt, by = "geneID")

        # add categories into dataHeat
        dataHeat <- merge(
            dataHeat_cat, dataHeat, by = c("geneID","supertaxon"), all.x = TRUE
        )
    }

    return(dataHeat)
}

#' Create data for plotting profiles (from raw input to final dataframe)
#' @description Create data needed for plotting phylogenetic profile
#' from raw input file.
#' @usage from_input_to_profile(raw_input, rank_name, ref_taxon, taxa_tree,
#'     var1_aggregate_by, var2_aggregate_by, percent_cutoff,
#'     coortholog_cutoff_max, var1_cutoff, var2_cutoff, var1_relation,
#'     var2_relation, group_by_cat, cat_dt)
#' @param raw_input input file (in long, wide, multi-fasta or orthoxml format)
#' @param rank_name taxonomy rank (e.g. "species","phylum",...)
#' @param ref_taxon selected reference taxon
#' @param taxa_tree input taxonomy tree (optional)
#' @param var1_aggregate_by aggregate method for VAR1 (min, max, mean or median)
#' @param var2_aggregate_by aggregate method for VAR2 (min, max, mean or median)
#' @param percent_cutoff min and max cutoffs for percentage of species present
#' in a supertaxon
#' @param coortholog_cutoff_max maximum number of co-orthologs allowed
#' @param var1_cutoff min and max cutoffs for var1
#' @param var2_cutoff min and max cutoffs for var2
#' @param var1_relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param var2_relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param group_by_cat group genes by their categories (TRUE or FALSE)
#' @param cat_dt dataframe contains gene categories
#' @return dataframe for generating profile plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{create_long_matrix}}, \code{\link{get_input_taxa_id}},
#' \code{\link{get_input_taxa_name}}, \code{\link{sort_input_taxa}},
#' \code{\link{parse_info_profile}}, \code{\link{reduce_profile}},
#' \code{\link{create_profile_data}}
#' @examples
#' raw_input <- system.file(
#'     "extdata", "test.main.long", package = "phyloprofile", mustWork = TRUE
#' )
#' rank_name <- "class"
#' ref_taxon <- "Mammalia"
#' taxa_tree <- NULL
#' var1_aggregate_by <- "max"
#' var2_aggregate_by <- "mean"
#' percent_cutoff <- c(0.0, 1.0)
#' coortholog_cutoff_max <- 10
#' var1_cutoff <- c(0.75, 1.0)
#' var2_cutoff <- c(0.5, 1.0)
#' var1_relation <- "protein"
#' var2_relation <- "species"
#' group_by_cat <- FALSE
#' cat_dt <- NULL
#' from_input_to_profile(
#'     raw_input,
#'     rank_name,
#'     ref_taxon,
#'     taxa_tree,
#'     var1_aggregate_by,
#'     var2_aggregate_by,
#'     percent_cutoff,
#'     coortholog_cutoff_max,
#'     var1_cutoff,
#'     var2_cutoff,
#'     var1_relation,
#'     var2_relation,
#'     group_by_cat,
#'     cat_dt
#' )

from_input_to_profile <- function(
    raw_input,
    rank_name,
    ref_taxon,
    taxa_tree,
    var1_aggregate_by,
    var2_aggregate_by,
    percent_cutoff,
    coortholog_cutoff_max,
    var1_cutoff,
    var2_cutoff,
    var1_relation,
    var2_relation,
    group_by_cat,
    cat_dt
) {
    # convert raw input into long format
    input_df <- phyloprofile::create_long_matrix(raw_input)

    # get input taxon IDs and names
    input_taxonID <- get_input_taxa_id(input_df)
    input_taxonName <- get_input_taxa_name(rank_name, input_taxonID)

    # sort input taxa based on selected reference taxon or input taxonomy tree
    sortedtaxa_list <- sort_input_taxa(
        input_taxonID, input_taxonName, rank_name, ref_taxon, taxa_tree
    )

    # parse info (additional values...) into profile df
    fullMdData <- parse_info_profile(
        input_df,
        sorted_input_taxa = sortedtaxa_list,
        var1_aggregate_by, var2_aggregate_by
    )

    # reduce profile df into supertaxon level
    dataSupertaxa <- reduce_profile(fullMdData)
    
    # cutoffs
    percent_cutoff_min <- percent_cutoff[1]
    percent_cutoff_max <- percent_cutoff[2]
    var1_cutoff_min <- var1_cutoff[1]
    var1_cutoff_max <- var1_cutoff[2]
    var2_cutoff_min <- var2_cutoff[1]
    var2_cutoff_max <- var2_cutoff[2]

    # create final df
    dataHeat <- create_profile_data(superTaxon_data = dataSupertaxa,
                                    ref_taxon = ref_taxon,
                                    percent_cutoff,
                                    coortholog_cutoff_max,
                                    var1_cutoff,
                                    var2_cutoff,
                                    var1_relation,
                                    var2_relation,
                                    group_by_cat = FALSE,
                                    cat_dt = NULL)

    return(dataHeat)
}
