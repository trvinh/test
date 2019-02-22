#' Parse domain input file
#' @export
#' @param seed seed ID(s)
#' @param input_file name of input data (demo data, file name or path to folder)
#' @param type type of data (demo, file or folder)
#' @importFrom utils read.csv
#' @return A dataframe containing protein domain info
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{get_domain_online}}, \code{\link{get_domain_file}},
#' \code{\link{get_domain_folder}}
#' @examples
#' seed <- "OG_1009"
#' input_file <- system.file(
#'     "extdata", "domain_files/OG_1009.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' type <- "file"
#' parse_domain_input(seed, input_file, type)

parse_domain_input <- function(seed, input_file, type) {
    domains <- data.frame()
    file <- NULL

    # get domain file(s) from online data
    if (type == "demo") {
        file <- get_domain_online(input_file)
    }
    # or from single file
    else if (type == "file") {
        file <- get_domain_file(input_file)
    }
    # or from a domain folder
    else {
        file <- get_domain_folder(seed, input_file)
        if (file == "noSelectHit") return("noSelectHit")
        else if (file == "noFileInFolder") return("noFileInFolder")
    }

    # parse domain file
    # for demo data
    if (type == "demo") {
        domain_df <- as.data.frame(read.csv(file,
                                            sep = "\t",
                                            header = FALSE,
                                            comment.char = "",
                                            stringsAsFactors = FALSE,
                                            quote = ""))
        domains <- rbind(domains, domain_df)

        if (ncol(domains) == 5) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end"
            )
        } else if (ncol(domains) == 6) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end", "weight"
            )
        } else if (ncol(domains) == 7) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end", "weight", "path"
            )
        }
    }
    # for user input data
    else {
        if (file != FALSE) {
            exeptions <- c("noFileInput", "noSelectHit", "noFileInFolder")
            if (!(file %in% exeptions)) {
                domain_df <- as.data.frame(read.table(
                    file,
                    sep = "\t",
                    header = FALSE,
                    comment.char = ""
                ))
                domains <- rbind(domains, domain_df)
            }
        }

        if (ncol(domains) == 5) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end"
            )
        } else if (ncol(domains) == 6) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end"
            )
        } else if (ncol(domains) == 7) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end",
                "weight"
            )
        } else if (ncol(domains) == 8) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end",
                "weight", "path"
            )
        } else {
            return("ERR")
        }
    }

    if (nrow(domains) == 0) return("ERR_0")

    domains$seedID <- as.character(domains$seedID)
    domains$orthoID <- as.character(domains$orthoID)
    domains$seedID <- gsub("\\|",":",domains$seedID)
    domains$orthoID <- gsub("\\|",":",domains$orthoID)

    return(domains)
}

#' Get domain file(s) for online data set
#' @param demo_data demo data name (either lca-micros or ampk-tor)
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' get_domain_online("lca-micors")
#' }

get_domain_online <- function(demo_data) {
    if (demo_data == "lca-micros") {
        file_domain <- {
            suppressWarnings(
                paste0(
                    "https://github.com/BIONF/phyloprofile-data/blob/master/",
                    "demo/domain_files/concatenate.domains?raw=true"
                )
            )
        }
    } else {
        file_domain <- {
            suppressWarnings(
                paste0(
                    "https://raw.githubusercontent.com/BIONF/phyloprofile-data",
                    "/master/expTestData/ampk-tor/ampk-tor.domains_F?raw=true"
                )
            )
        }
    }

    return(file_domain)
}

#' Get domain file from a single (concatenate) file
#' @param input_file concatenate domain file
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' input_file <- system.file(
#'     "extdata", "domain_files/OG_1009.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' get_domain_file("lca-micors")
#' }

get_domain_file <- function(input_file) {
    file_domain <- input_file
    return(file_domain)
}

#' Get domain files from a folder
#' @param seed seed ID
#' @param domain_path path to domain folder
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' domain_path <- paste0(
#'     path.package("phyloprofile", quiet = FALSE), "extdata/domain_files"
#' )
#' get_domain_online("OG_1009", domain_path)
#' }

get_domain_folder <- function(seed, domain_path){
    if (is.null(seed)) {
        file_domain <- "noSelectHit"
    } else {
        # check file extension
        all_extension <- c("txt", "csv", "list", "domains", "architecture")
        flag <- 0

        for (i in seq_len(length(all_extension))) {
            file_domain <- paste0(
                domain_path, "/", seed, ".", all_extension[i]
            )
            if (file.exists(file_domain) == TRUE) {
                flag <- 1
                break()
            }
        }

        if (flag == 0) {
            file_domain <- "noFileInFolder"
        }
    }

    return(file_domain)
}
