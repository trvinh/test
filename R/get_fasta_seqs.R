#' Get fasta sequences
#' @export
#' @usage get_fasta_seqs(data_in, filein, demo_data, input_type_fasta,
#'     concat_fasta, path, dir_format, file_ext, id_format)
#' @param data_in a dataframe of the input phylogenetic profiles, that contains
#' at least 3 columns: geneIDs, orthoIDs and ncbiIDs
#' @param filein main input file (for checking input type)
#' @param demo_data name of demo data set (either "ampk-tor", "lca-micros",
#' or "none")
#' @param input_type_fasta source of fasta sequences ("Concatenated fasta file"
#' or "Fasta folder")
#' @param concat_fasta input concatenated fasta file
#' @param path path to fasta folder
#' @param dir_format directory format (either "path/speciesID.fa*" or
#' "path/speciesID/speciesID.fa*")
#' @param file_ext fasta file extension ("fa", "fasta", "fas" or "txt")
#' @param id_format fasta header format (">speciesID:seqID",
#' ">speciesID@seqID", ">speciesID|seqID" or only "seqID")
#' @return A dataframe contains fasta sequences
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{check_input_validity}}, \code{\link{str_reverse}},
#' \code{\link{main_long_raw}}
#' @examples
#' data("main_long_raw", package="phyloprofile")
#' data_in <- main_long_raw
#' filein <- system.file(
#'     "extdata", "test.main.long", package = "phyloprofile", mustWork = TRUE
#' )
#' demo_data <- "none"
#' input_type_fasta <- "Concatenated fasta file"
#' concat_fasta <- system.file(
#'     "extdata", "fasta_files/concatenatedFile.fa",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' # the following parameters are conly required if upload a fasta folder!
#' path <- NULL
#' dir_format <- NULL
#' file_ext <- NULL
#' id_format <- NULL
#' # get fasta sequences
#' get_fasta_seqs(
#'     data_in, filein, demo_data,
#'     input_type_fasta,
#'     concat_fasta,
#'     path,
#'     dir_format,
#'     file_ext,
#'     id_format
#' )

get_fasta_seqs <- function(
    data_in, filein, demo_data,
    input_type_fasta,
    concat_fasta,
    path,
    dir_format,
    file_ext,
    id_format
){

    fasta_out_df <- data.frame()

    # check main input
    if (!is.null(filein)) {
        input_type <- check_input_validity(filein)
    } else{
        input_type <- "NA"
    }

    # get seqs for AMPK-TOR and microsporidia ONLINE demo data -----------------
    if (demo_data == "ampk-tor" | demo_data == "lca-micros") {
        fasta_url <-
            paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/demo/fasta_file/concatenatedSeq.fa"
            )
        if (demo_data == "ampk-tor") {
            fasta_url <-
                paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/expTestData/ampk-tor/",
                    "ampk-tor.extended.fa"
                )
        }

        if (RCurl::url.exists(fasta_url)) {
            # load fasta file
            fa_file <- as.data.frame(read.table(fasta_url,
                                                sep = "\t",
                                                header = FALSE,
                                                fill = TRUE,
                                                stringsAsFactors = FALSE,
                                                quote = ""))
            fa_df <- data.frame("seqID" = fa_file$V1[grepl(">", fa_file$V1)],
                                "seq" = fa_file$V1[!grepl(">", fa_file$V1)],
                                stringsAsFactors = FALSE)
            # get sequences
            for (j in seq_len(nrow(data_in))) {
                seq_id <- as.character(data_in$orthoID[j])
                group_id <- as.character(data_in$geneID[j])

                seq <- as.character(
                    fa_df$seq[fa_df$seqID == paste0(">", seq_id)]
                )
                fasta_out <- paste(paste0(">", seq_id), seq, sep = "\n")
                fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
            }
        } else {
            fasta_out <- paste0(fasta_url, " not found!!!")
            fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
        }
    }

    # get seqs for fasta main input --------------------------------------------
    if (input_type == "fasta") {
        file <- filein
        fasta_file <- Biostrings::readAAStringSet(file)

        seq_name <- names(fasta_file)
        sequence <- paste(fasta_file)
        # data frame contains all sequences from input file
        fa <- data.frame(seq_name, sequence)

        for (j in seq_len(nrow(data_in))) {
            seq_id <- paste0(
                as.character(data_in$geneID[j]), "|ncbi",
                as.character(data_in$ncbiID[j]), "|",
                as.character(data_in$orthoID[j])
            )

            seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]

            if (length(seq[1]) < 1) {
                fasta_out <- paste0(seq_id,
                                    " not found in ",
                                    file,
                                    "! Please check again!")
            } else{
                fasta_out <- paste(paste0(">", seq_id), seq[1], sep = "\n")
            }
            fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
        }
    }

    # get seqs for main input in other formats ---------------------------------
    else {
        # * get seqs from concatenated fasta file ------------------------------
        if (demo_data == "none" &
            input_type_fasta == "Concatenated fasta file") {
            if (!is.null(concat_fasta)) {
                fas_in <- concat_fasta
                file <- toString(fas_in)

                # read fasta file and save sequences into dataframe
                fasta_file <- Biostrings::readAAStringSet(file)

                seq_name <- names(fasta_file)
                sequence <- paste(fasta_file)
                # data frame contains all sequences from input file
                fa <- data.frame(seq_name, sequence)

                # get selected sequences
                for (j in seq_len(nrow(data_in))) {
                    seq_id <- as.character(data_in$orthoID[j])
                    group_id <- as.character(data_in$geneID[j])
                    seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
                    flag <- 1
                    if (is.na(seq)) {
                        seq_id <- paste0(
                            as.character(data_in$geneID[j]), "|ncbi",
                            as.character(data_in$ncbiID[j]), "|",
                            as.character(data_in$orthoID[j])
                        )
                        seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
                        flag <- 0
                    }
                    if (length(seq[1]) < 1) {
                        fasta_out <-
                            paste0(
                                seq_id,
                                " not found in ",
                                file,
                                "! Please check the header format in
                                FASTA file!"
                            )
                    } else {
                        if (!is.na(seq[1])) {
                            if (flag == 1) {
                                fasta_out <- paste(
                                    paste0(">", seq_id), seq[1], sep = "\n"
                                )
                            } else {
                                fasta_out <- paste(
                                    paste0(">", seq_id), seq[1], sep = "\n"
                                )
                            }

                        } else {
                            fasta_out <-
                                paste0(
                                    seq_id,
                                    " not found in uploaded FASTA file!!! ",
                                    "Please check again!!!"
                                )
                        }
                    }
                    fasta_out_df <- rbind(
                        fasta_out_df, as.data.frame(fasta_out)
                    )
                }
            } else {
                if (input_type != "fasta") {
                    fasta_out <- {
                        paste0(
                            "Please provide FASTA file(s) in ",
                            "Input & settings page!"
                        )
                    }
                    fasta_out_df <- rbind(
                        fasta_out_df, as.data.frame(fasta_out)
                    )
                }
            }
        }

        # * get seqs from an folder --------------------------------------------
        if (demo_data == "none"
            & input_type_fasta == "Fasta folder"
            & input_type != "fasta") {
            if (path != "") {
                # get list of species IDs
                if (id_format == 1) {
                    spec_df <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                str_reverse(as.character(data_in$orthoID)),
                                ":", 2
                            )
                        )
                    spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
                } else if (id_format == 2) {
                    spec_df <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                str_reverse(as.character(data_in$orthoID)),
                                "@", 2
                            )
                        )
                    spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
                } else if (id_format == 3) {
                    spec_df <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                str_reverse(as.character(data_in$orthoID)),
                                "|", 2
                            )
                        )
                    spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
                }

                # read all specices FASTA files at once
                fa <- data.frame()
                if (id_format == 4) {
                    file_list <- list.files(path, pattern = file_ext)
                    for (file in file_list){
                        file_with_path = paste0(path, "/", file)
                        fasta_file <-
                            Biostrings::readAAStringSet(file_with_path)

                        seq_name <- names(fasta_file)
                        sequence <- paste(fasta_file)
                        # data frame contains all sequences from input file
                        fa <- rbind(fa, data.frame(seq_name, sequence))
                    }
                } else {
                    for (
                        i in
                        seq_len(length(levels(as.factor(spec_df$specID))))
                    ) {
                        spec_id <- as.character(
                            levels(as.factor(spec_df$specID))[i]
                        )

                        # full path fasta file
                        file <- paste0(path, "/", spec_id, ".", file_ext)
                        if (dir_format == 2) {
                            file <- paste0(
                                path, "/", spec_id, "/", spec_id, ".", file_ext
                            )
                        }

                        # read fasta file and save sequences into dataframe
                        if (file.exists(file)) {
                            fasta_file <- Biostrings::readAAStringSet(file)
                            seq_name <- names(fasta_file)
                            sequence <- paste(fasta_file)
                            # data frame contains all sequences from input file
                            fa <- rbind(fa, data.frame(seq_name, sequence))
                        }
                    }
                }

                # now get selected sequences
                if (nrow(fa) > 0) {
                    for (j in seq_len(nrow(data_in))) {
                        seq_id <- as.character(data_in$orthoID[j])
                        group_id <- as.character(data_in$geneID[j])

                        seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]

                        if (length(seq[1]) < 1) {
                            fasta_out <- {
                                paste0(
                                    seq_id,
                                    " not found in ",
                                    file,
                                    "! Please check id_format in FASTA config ",
                                    "again!"
                                )
                            }
                        } else {
                            fasta_out <- paste(
                                paste0(">", seq_id), seq[1], sep = "\n"
                            )
                        }
                        fasta_out_df <- rbind(
                            fasta_out_df, as.data.frame(fasta_out)
                        )
                    }
                } else {
                    fasta_out <-
                        paste0(
                            "No fasta file has been found in ",
                            path,
                            "!!! Please check the full path to FASTA folder ",
                            "and the id_format (header format) in FASTA config",
                            " again!!!"
                        )
                    fasta_out_df <- rbind(
                        fasta_out_df, as.data.frame(fasta_out)
                    )
                }
            } else {
                fasta_out <- {
                    paste0(
                        "Please provide FASTA files in Input & settings page!"
                    )
                }
                fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
            }
        }
    }
    # remove duplicated sequences
    fasta_out_df <- fasta_out_df[!duplicated(fasta_out_df), ]

    return(fasta_out_df)
}

#' Reverse string
#' @return a reversed string
#' @param x input string
#' @examples
#' \dontrun{
#' str_reverse("abc")
#' }

str_reverse <- function(x) {
    vapply(
        lapply(strsplit(x, NULL), rev), paste, collapse = "",
        FUN.VALUE = character(1)
    )
}
