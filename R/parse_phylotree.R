# Functions for working with phylogenetic trees

#' Check the validity of input newick tree
#' @export
#' @param tree input newick tree
#' @param input_taxonID list of all input taxon ID
#' @return checking result (1 = missing parenthesis; 2 = missing comma;
#' 3 = tree has singleton; or list of missing taxa in profile)
#' @importFrom stringr regex
#' @importFrom stringr str_count
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{get_input_taxa_id}} for getting input taxon IDs,
#' \code{\link{pp_tree}} for an example of input tree
#' @examples
#' data("pp_tree", package="phyloprofile")
#' check_newick(pp_tree, c("ncbi3702", "ncbi3711", "ncbi7029"))

check_newick <- function(tree, input_taxonID){
    # get tree structure
    tree_struc <- gsub(stringr::regex("\\w"), "", as.character(tree$V1))

    open <- stringr::str_count(tree_struc, "\\(")
    close <- stringr::str_count(tree_struc, "\\)")
    comma <- stringr::str_count(tree_struc, "\\,")
    singleton <- stringr::str_count(tree_struc, "\\(\\)")

    if (singleton > 0) {
        return(3) # tree contains singleton
    }
    if (open != close) {
        return(1) # missing parenthesis
    } else {
        if ((comma - open) > 1 | (comma - open) < 0) {
            # return(2) # missing comma
        } else {
            # get list of tips
            node_string <- gsub(regex("\\W+"), "#", as.character(tree$V1))
            node_list <- unlist(strsplit(node_string, "#"))
            # list of input taxa
            input_taxa <- input_taxonID

            missing_taxa <- list()
            j <- 1
            for (i in seq_len(length(node_list))) {
                if (nchar(node_list[i]) > 0 & !(node_list[i] %in% input_taxa)) {
                    missing_taxa[[j]] <- node_list[i]
                    j <- j + 1
                }
            }
            if (length(missing_taxa) > 0) {
                # contains taxa that not exist in main input
                return(paste(missing_taxa, collapse = "; "))
            } else {
                return(0)
            }
        }
    }

    return(0)
}

#' Create rooted tree from a taxonomy matrix
#' @export
#' @param df data frame contains taxonomy matrix used for creating tree
#' @param root_taxon taxon used for rooting tree
#' @importFrom stats hclust
#' @importFrom ape as.phylo
#' @importFrom ape root
#' @return rooted tree based on root_taxon
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{taxa2dist}} for distance matrix generation from a
#' taxonomy matrix, \code{\link{get_taxonomy_matrix}} for getting taxonomy
#' matrix, \code{\link{pp_taxonomy_matrix}} for a demo taxonomy matrix data
#' @examples
#' data("pp_taxonomy_matrix", package = "phyloprofile")
#' # prepare matrix for calculating distances
#' distDf <- subset(pp_taxonomy_matrix, select = -c(ncbiID, fullName))
#' row.names(distDf) <- distDf$abbrName
#' distDf <- distDf[, -1]
#' # create taxonomy tree rooted by ncbi10090
#' create_rooted_tree(distDf, "ncbi10090")

create_rooted_tree <- function(df, root_taxon){
    # calculate distance matrix
    taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
    # create tree
    tree <- ape::as.phylo(stats::hclust(taxdis))
    # root tree
    tree <- ape::root(tree, outgroup = root_taxon, resolve.root = TRUE)
    # return
    return(tree)
}

#' Get sorted supertaxon list based on a rooted taxonomy tree
#' @export
#' @param tree rooted taxonomy tree
#' @return list of taxa sorted according to the taxonomy tree
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{pp_taxonomy_matrix}} for a demo taxonomy matrix data
#' @examples
#' data("pp_taxonomy_matrix", package = "phyloprofile")
#' # prepare matrix for calculating distances
#' distDf <- subset(pp_taxonomy_matrix, select = -c(ncbiID, fullName))
#' row.names(distDf) <- distDf$abbrName
#' distDf <- distDf[, -1]
#' # create taxonomy tree rooted by ncbi10090
#' rooted_tree <- create_rooted_tree(distDf, "ncbi10090")
#' # get taxon list sorted from tree
#' sort_taxa_from_tree(rooted_tree)

sort_taxa_from_tree <- function(tree){
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips <- tree$edge[is_tip, 2]
    taxon_list <- rev(tree$tip.label[ordered_tips])
    return(taxon_list)
}
