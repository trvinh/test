#' @title Search NCBI taxonomy IDs for a list of taxon names
#' @param taxa list of taxon names
#' @importFrom taxize get_uid
#' @importFrom taxize gnr_resolve
#' @return dataframe contains NCBI taxonomy IDs for all input taxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @examples
#' taxa <- c("homo sapiens")
#' search_taxonID_online(taxa)

search_taxonID_online <- function(taxa){
    id_df <- data.frame("name" = character(),
                        "new_name" = character(),
                        "id" = character(),
                        "type" = character(),
                        stringsAsFactors = FALSE)
    for (i in seq_len(length(taxa))) {
        id_df[i,] <- get_taxon_id(taxa[i])
    }
    return(id_df)
}

get_taxon_id <- function(tax_name){
    id <- taxize::get_uid(sciname = tax_name)[1]

    id_df <- data.frame("name" = character(),
                        "new_name" = character(),
                        "id" = character(),
                        "type" = character(),
                        stringsAsFactors = FALSE)

    if (is.na(id)) {
        temp <- taxize::gnr_resolve(names = as.character(tax_name))
        if (nrow(temp) > 0) {
            new_id <- taxize::get_uid(sciname = temp[1, 3])[1]
            if (is.na(new_id)) {
                id_df[1, ] <- c(as.character(tax_name),
                                as.character(temp[1, 3]),
                                paste0("NA"), "notfound")
            } else {
                id_df[1, ] <- c(as.character(tax_name),
                                as.character(temp[1, 3]),
                                paste0("ncbi", new_id),
                                "notfound")
            }
        } else {
            id_df[1, ] <- c(as.character(tax_name),
                            paste0("no alternative"),
                            paste0("NA"),
                            "notfound")
        }
    } else {
        id_df[1, ] <- c(as.character(tax_name),
                        "NA",
                        paste0("ncbi", id),
                        "retrieved")
    }
    return(id_df)
}
