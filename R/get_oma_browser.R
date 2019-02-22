# Get OMA information functions
# using OmaDB package (https://github.com/klarakaleb/OmaDB)

#' Check OMA IDs
#' @description Check if input IDs are valid IDs for OMA Browser
#' (either OMA IDs or UniProt IDs)
#' @export
#' @param ids list of ids needs to be checked
#' @return list of invalid IDs (not readable for OMA)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' check_oma_id("HUMAN29398")

check_oma_id <- function(ids) {
    invalid <- list()
    for (id in ids) {
        id <- as.character(id)
        data <- OmaDB::getData("protein", id)
        if (is.null(data$entry_nr)) invalid <- c(invalid, id)
    }
    return(invalid)
}

#' Get OMA members
#' @description Get OMA orthologs for a seed protein from OMA Browser
#' @export
#' @param id ID of the seed protein (OMA or UniProt ID)
#' @param ortho_type type of OMA orthologs: either "HOG", "OG"
#' (orthologous group) or "PAIR" (orthologous pair - CURRENTLY NOT WORKING)
#' @return list of ortholog members
#' @author Carla Mölbert {carla.moelbert@gmx.de}
#' @examples
#' get_oma_members("HUMAN29397", "OG")

get_oma_members <- function(id, ortho_type) {
    # get the members of the Hierarchical Orthologous Group
    if (ortho_type == "HOG") {
        members <- suppressWarnings(
            OmaDB::getHOG(id = id, level = "root", members = TRUE)$members$omaid
        )
    }
    # get the members of the Ortholoug group
    else if (ortho_type == "OG") {
        members <- suppressWarnings(
            OmaDB::getData(type = "group", id = id)$members$omaid
        )
    }
    # # get the members of the Orthologous Pair
    # else if (ortho_type == "PAIR") {
    #   members <- OmaDB::resolveURL(OmaDB::getData(type = "protein",
    #                                               id = id)$orthologs)$omaid
    #   # add query ID into output list
    #   seed <- OmaDB::getData("protein",id)$omaid
    #   members <- c(seed,members)
    # }

    return(members)
}

#' Get domain annotation from OMA Browser
#' @description Get domain annotation from OMA Browser based on a URL or a
#' raw data frame contains annotation info from OMA
#' @param domainURL URL address for domain annotation of ONE OMA id
#' @return data frame contains feature names, start and end positions
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' get_domain_from_url("https://omabrowser.org/api/protein/7916808/domains/")
#' }

get_domain_from_url <- function(domainURL) {
    if (grepl("https://", domainURL)) {
        domains <- OmaDB::resolveURL(domainURL)$regions
    } else {
        domains <- domainURL$regions
    }

    domains$feature <- NA
    domains$start <- NA
    domains$end <- NA
    for (i in seq_len(nrow(domains))) {
        pos <- unlist(strsplit(domains$location[i], ":"))
        domains[i,]$start <- pos[1]
        domains[i,]$end <- pos[2]

        if (nchar(domains$name[i]) > 0) {
            domains[i,]$feature <- paste0(
                domains$source[i],"_",domains$domainid[i],
                " (",domains$name[i],")"
            )
        } else {
            domains[i,]$feature <- paste0(
                domains$source[i],"_",domains$domainid[i]
            )
        }
        domains[i,]$feature <- gsub("#", "-", domains[i,]$feature)
    }
    return(domains[, c("feature","start","end")])
}

#' Get taxonomy ID, sequence and annotation for one OMA sequence
#' @param id oma ID of a ortholog
#' @return data frame
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' get_data_for_one_ortholog(HUMAN29397)
#' }

get_data_for_one_ortholog <- function(id) {
    omaDf <- data.frame(
        "ortho_id" = character(),
        "taxon_id" = character(),
        "seq" = character(),
        "length" = numeric(),
        "domains" = character(),
        stringsAsFactors = FALSE
    )

    # get ncbi taxonomy id
    spec_name <- substr(id, 1, 5)
    taxon_id <- paste0(
        "ncbi", OmaDB::getTaxonomy(members = spec_name, newick = FALSE)$id
    )
    # get raw data
    raw <- OmaDB::getData("protein",id)

    # get sequence
    seq <- as.character(raw$sequence)

    # get sequence length
    length <- raw$sequence_length

    # get annotation
    raw_domains <- suppressWarnings(raw$domains)
    domainDf <- suppressWarnings(get_domain_from_url(raw_domains))
    domainDf_join <- c(domainDf, sep = "#")
    domains <- paste(unlist(do.call(paste, domainDf_join)), collapse = "\t")

    # return data frame contains all info
    omaDf[1,] <- c(id, taxon_id, seq, length, domains)
    return(omaDf)
}

#' Get OMA info for a query protein and its orthologs
#' @description Get taxonomy IDs, sequences and annotations for an OMA
#' orthologous group (or OMA HOG).
#' @export
#' @param seed_id protein query ID in OMA or UniProt format
#' @param ortho_type type of OMA orthologs
#' @return data frame contains info for all sequences of that OMA group (or HOG)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' get_data_for_one_oma("HUMAN29397", "OG")

get_data_for_one_oma <- function(seed_id, ortho_type){
    final_omaDf <- data.frame()

    # get members
    members <- get_oma_members(seed_id, ortho_type)
    oma_seed_id <- OmaDB::getData("protein",seed_id)$omaid

    # get all data
    j <- 1
    for (ortho in members) {
        orthoDf <- get_data_for_one_ortholog(ortho)
        orthoDf$seed <- seed_id
        if (ortho == oma_seed_id) {
            orthoDf$ortho_id <- seed_id
        }
        final_omaDf <- rbind(final_omaDf, orthoDf)

        p <- j / length(members) * 100
        svMisc::progress(p)
        j <- j + 1
    }

    return(final_omaDf)
}

#' Create phylogenetic profile from a raw OMA dataframe
#' @export
#' @param final_oma_df raw OMA data for a list of input IDs
#' @return phylogenetic profiles in long format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' oma_data <- get_data_for_one_oma("HUMAN29397", "OG")
#' create_profile_from_oma(oma_data)
#' @seealso \code{\link{get_data_for_one_oma}}

create_profile_from_oma <- function(final_oma_df) {
    profile_df <- final_oma_df[, c("seed", "taxon_id", "ortho_id")]
    colnames(profile_df) <- c("geneID", "ncbiID", "orthoID")
    return(profile_df[!duplicated(profile_df), ])
}

#' Create domain annotation dataframe for one OMA protein
#' @param domain_id protein domain ID
#' @param ortho_id protein ID
#' @param length protein length
#' @param domain_list list of all domains and their positions for this protein
#' @return domain annotation in a dataframe
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_domain_df <- function(domain_id, ortho_id, length, domain_list) {
    domain_df = data.frame(
        seedID = character(),
        orthoID = character(),
        length = numeric(),
        feature = character(),
        start = numeric(),
        end = numeric(),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(length(domain_list))) {
        anno_info <- strsplit(domain_list[i], "#")[[1]]
        domain_df[i,] <- c( domain_id, ortho_id, length, anno_info)
    }

    return(domain_df)
}

#' Create domain annotation dataframe from a raw OMA dataframe
#' @export
#' @param final_oma_df raw OMA data for a list of input IDs
#' @return domain annotation in a dataframe to input into PhyloProfile, which
#' contains seed ID, ortholog ID, ortholog length, annotated feature, start
#' and end position of that feature.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' oma_data <- get_data_for_one_oma("HUMAN29397", "OG")
#' get_all_domains_oma(oma_data)
#' @seealso \code{\link{get_data_for_one_oma}}

get_all_domains_oma <- function(final_oma_df) {
    oma_domain_df <- data.frame()
    for (i in seq_len(nrow(final_oma_df))) {
        domainID <- paste0(
            final_oma_df[i,]$seed, "#", final_oma_df[i,]$ortho_id
        )

        seed_line <-
            final_oma_df[final_oma_df$ortho_id == final_oma_df[i,]$seed, ]
        seed_domains <- strsplit(as.character(seed_line$domains), "\t")[[1]]
        seed_domainDf <- create_domain_df(
            domainID, seed_line$ortho_id, seed_line$length, seed_domains
        )
        oma_domain_df <- rbind(oma_domain_df, seed_domainDf)

        ortho_domains <-
            strsplit(as.character(final_oma_df[i,]$domains), "\t")[[1]]
        ortho_domainDf <- create_domain_df(
            domainID,
            final_oma_df[i,]$ortho_id,
            final_oma_df[i,]$length,
            ortho_domains
        )
        oma_domain_df <- rbind(oma_domain_df, ortho_domainDf)
    }

    oma_domain_df$length <- as.numeric(oma_domain_df$length)
    oma_domain_df$start <- as.numeric(oma_domain_df$start)
    oma_domain_df$end <- as.numeric(oma_domain_df$end)
    return(oma_domain_df)
}

#' Get all fasta sequences from a raw OMA dataframe
#' @export
#' @param final_oma_df raw OMA data for a list of input IDs
#' @return dataframe contains all fasta sequences
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' oma_data <- get_data_for_one_oma("HUMAN29397", "OG")
#' get_all_fasta_oma(oma_data)
#' @seealso \code{\link{get_data_for_one_oma}}

get_all_fasta_oma <- function(final_oma_df) {
    oma_fasta_df <- data.frame()

    fasta_df <- final_oma_df[, c("ortho_id", "seq")]
    for (i in seq_len(nrow(fasta_df))) {
        seq_id <- as.character(fasta_df$ortho_id[i])
        seq <- as.character(fasta_df$seq[i])
        fasta_out <- paste(paste0(">", seq_id), seq, sep = "\n")
        oma_fasta_df <- rbind(oma_fasta_df, as.data.frame(fasta_out))
    }

    return(oma_fasta_df[!duplicated(oma_fasta_df), ])
}

#' Get selected fasta sequences from a raw OMA dataframe
#' @export
#' @param final_oma_df raw OMA data for a list of input IDs
#' @param seq_id sequence need to be returned
#' @return required sequence in fasta format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' oma_data <- get_data_for_one_oma("HUMAN29397", "OG")
#' get_selected_fasta_oma(oma_data, "HUMAN29397")
#' @seealso \code{\link{get_data_for_one_oma}}

get_selected_fasta_oma <- function(final_oma_df, seq_id) {
    selected_df <- subset(
        final_oma_df[, c("ortho_id", "seq")],
        final_oma_df$ortho_id == seq_id
    )
    header <- paste0(">", selected_df$ortho_id[1])
    return(paste(header, selected_df$seq[1], sep = "\n"))
}
