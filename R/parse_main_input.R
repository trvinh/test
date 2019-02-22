# Parse main input functions
# These functions check the validity of the main input file
# and convert different input format into long format.

#' Check the validity of the main input file
#' @description Check if input file has one of the following format: orthoXML,
#' multiple FASTA, wide or long matrix, or a list of OMA IDs.
#' @export
#' @param filein input file
#' @return The format of the input file format, or type of error
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{check_oma_id}}
#' @examples
#' filein <- system.file(
#'     "extdata", "test.main.wide", package = "phyloprofile", mustWork = TRUE
#' )
#' check_input_validity(filein)

check_input_validity <- function(filein) {
    input_dt <- as.data.frame(read.table(
        file = filein,
        sep = "\t",
        header = FALSE,
        check.names = FALSE,
        comment.char = "",
        fill = TRUE
    ))

    if (is.na(input_dt[1, ncol(input_dt)])) {
        return("moreCol")
    } else {
        names(input_dt) <- as.character(unlist(input_dt[1, ]))

        # XML format (starts with <?xml)
        if (grepl("<?xml", colnames(input_dt)[1])) {
            return("xml")
        }
        # FASTA format (starts with ">" )
        else if (grepl(">", colnames(input_dt)[1]) == TRUE) {
            return("fasta")
        }
        # LONG or WIDE format (starts with "geneID")
        else {
            if (grepl("geneID", colnames(input_dt)[1])) {
                # LONG format
                if (is.na(pmatch("ncbi", colnames(input_dt)[3])) ||
                    is.na(pmatch("ncbi", colnames(input_dt)[4])) ||
                    is.na(pmatch("ncbi", colnames(input_dt)[5]))) {
                    return("long")
                }
                # WIDE format
                else {
                    tmp <- input_dt[input_dt == ""][1]
                    if (!is.na(tmp) & tmp == "") {
                        return("emptyCell")
                    } else {
                        return("wide")
                    }
                }
            }
            # OMA ids
            else {
                invalid_oma <- check_oma_id(levels(input_dt[,1]))
                if (length(invalid_oma) == 0) {
                    return("oma")
                } else {
                    return(invalid_oma)
                }
            }
        }
    }
}

#' Parse orthoXML input file
#' @param input_file input file in xml format
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' input_file <- system.file(
#'     "extdata", "test.main.xml", package = "phyloprofile", mustWork = TRUE
#' )
#' xml_parser(input_file)
#' }

xml_parser <- function(input_file){
    path <- paste(
        system.file(package="phyloprofile"),
        "phyloprofile/scripts/orthoxmlParser.py",
        sep="/"
    )
    cmd <- paste("python ", path, " -i ", input_file, sep = "")

    df_in <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

    # the first row will be the header
    colnames(df_in) <- as.character(unlist(df_in[1, ]))

    df_in <- subset(df_in[df_in$geneID != "geneID", ])
    df_in <- droplevels(df_in)

    return(df_in)
}

#' Parse multi-fasta input file
#' @param input_file input multiple fasta file
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' input_file <- system.file(
#'     "extdata", "test.main.fasta", package = "phyloprofile", mustWork = TRUE
#' )
#' fasta_parser(input_file)
#' }

fasta_parser <- function(input_file){
    path <- paste(
        system.file(package="phyloprofile"),
        "phyloprofile/scripts/fastaParser.py",
        sep="/"
    )
    cmd <- paste("python ", path, " -i ", input_file, sep = "")

    var1 <- NULL
    var2 <- NULL

    df_in <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

    # the first row will be the header
    colnames(df_in) <- as.character(unlist(df_in[1, ]))
    df_in <- subset(df_in[df_in$geneID != "geneID", ])
    df_in <- droplevels(df_in)

    # remove var1 and var2 columns if they are all NAs
    if (all(is.na(df_in$var2))) {
        df_in <- subset(df_in, select = -c(var2) )
    }
    if (all(is.na(df_in$var1))) {
        df_in <- subset(df_in, select = -c(var1) )
    }

    return(df_in)
}

#' Transform input file in wide matrix into long matrix format
#' @param input_file input file in wide matrix format
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' input_file <- system.file(
#'     "extdata", "test.main.wide", package = "phyloprofile", mustWork = TRUE
#' )
#' wide_to_long(input_file)
#' }

wide_to_long <- function(input_file){
    wide_dataframe <- as.data.frame(read.table(
        file = input_file,
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        comment.char = ""
    ))
    long_dataframe <- data.frame()
    row_nr_long <- 0
    ncbi_ids <- colnames(wide_dataframe)

    for (row_nr in seq_len(nrow(wide_dataframe))) {
        geneID <- wide_dataframe[row_nr, 1]
        for (column_nr in 2:ncol(wide_dataframe)) {
            current_cell <- as.character(wide_dataframe[row_nr, column_nr])
            new_row_info <- unlist(strsplit(current_cell, "#"))
            row_nr_long <- row_nr_long + 1
            long_dataframe[row_nr_long, 1] <- geneID
            long_dataframe[row_nr_long, 2] <- ncbi_ids[column_nr]
            long_dataframe[row_nr_long, 3] <- new_row_info[1]
            long_dataframe[row_nr_long, 4] <- new_row_info[2]
            long_dataframe[row_nr_long, 5] <- new_row_info[3]
        }
    }

    colnames(long_dataframe) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    return(long_dataframe)
}

#' Create a long matrix format for all kinds of input file
#' @export
#' @param input_file input file in orthoXML, multiple FASTA, wide or long
#' matrix format.
#' @return A data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{xml_parser}}, \code{\link{fasta_parser}},
#' \code{\link{wide_to_long}}
#' @examples
#' input_file <- system.file(
#'     "extdata", "test.main.wide", package = "phyloprofile", mustWork = TRUE
#' )
#' create_long_matrix(input_file)

create_long_matrix <- function(input_file){
    if (input_file[1] == "lca-micros") {
        input_url <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "demo/test.main.long"
        )
        long_dataframe <-
            read.table(
                input_url,
                sep = "\t",
                header = TRUE,
                fill = TRUE,
                stringsAsFactors = FALSE
            )
    } else if (input_file[1] == "ampk-tor") {
        input_url <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "expTestData/ampk-tor/ampk-tor.phyloprofile"
        )
        long_dataframe <-
            read.table(
                input_url,
                sep = "\t",
                header = TRUE,
                fill = TRUE,
                stringsAsFactors = FALSE
            )
    } else {
        filein <- input_file
        if (is.null(filein)) return()
        input_type <- check_input_validity(filein)

        # XML
        if (input_type == "xml") {
            long_dataframe <- xml_parser(filein)
        }
        # FASTA
        else if (input_type == "fasta") {
            long_dataframe <- fasta_parser(filein)
        }
        # LONG
        else if (input_type == "long") {
            long_dataframe <- as.data.frame(read.table(
                file = filein,
                sep = "\t",
                header = TRUE,
                check.names = FALSE,
                comment.char = ""
            ))
        }
        # WIDE
        else if (input_type == "wide") {
            long_dataframe <- wide_to_long(filein)
        }
        else {
            return(NULL)
        }
    }

    # make sure all columns have the same type (factor)
    for (i in seq_len(ncol(long_dataframe))) {
        long_dataframe[, i] <- as.factor(long_dataframe[, i])
    }

    # long_dataframe$orthoID <- gsub("\\|",":",long_dataframe$orthoID)
    return(long_dataframe)
}
