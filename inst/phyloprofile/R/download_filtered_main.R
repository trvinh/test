#' Download filtered data from main profile
#'
#' @export
#' @param full processed main data (from reactive fn "get_data_filtered")
#' @param fasta fasta sequences (from reactive fn "main_fasta_download")
#' @param var1_id name of 1st variable (from input$var1_id)
#' @param var2_id name of 2nd variable (from input$var2_id)
#' @param var1 cutoff value of 1st variable (from input$var1)
#' @param var2 cutoff value of 2nd variable (from input$var2)
#' @param percent cutoff value of percentage species in each supertaxon
#' (from input$percent)
#' @return data of main profile for downloading
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

download_filtered_main_ui <- function(id) {
    ns <- NS(id)

    tabPanel(
        "Main data",
        column(
            4,
            checkboxInput(
                ns("get_representative_main"),
                strong(em("Download representative sequences")),
                value = FALSE,
                width = NULL
            )
        ),
        column(
            4,
            conditionalPanel(
                condition = {
                    sprintf("input['%s'] == true",
                            ns("get_representative_main"))
                },
                uiOutput(ns("ref_var_main.ui"))
            )
        ),
        column(
            4,
            conditionalPanel(
                condition = {
                    sprintf("input['%s'] == true",
                            ns("get_representative_main"))
                },
                radioButtons(
                    inputId = ns("ref_type_main"),
                    label = {
                        "Select representative by"
                    },
                    choices = list("max", "min"),
                    selected = "max",
                    inline = TRUE
                )
            )
        ),
        column(
            12,
            dataTableOutput(ns("filtered_main_data"))
        ),
        column(
            4,
            downloadButton(ns("download_data"),
                           "Download filtered data")
        ),
        column(
            4,
            downloadButton(ns("download_fasta"),
                           "Download FASTA sequences")
        ),
        column(
            4,
            downloadButton(ns("download_long"),
                           "Download data as PhyloProfile input format")
        )
    )
}

download_filtered_main <- function(input, output, session,
                                   data,
                                   fasta,
                                   var1_id, var2_id,
                                   var1, var2, percent){

    # render options for downloading -------------------------------------------
    output$ref_var_main.ui <- renderUI({
        ns <- session$ns
        if (nchar(var2_id()) < 1 & nchar(var1_id()) < 1) {
            radioButtons(
                inputId = ns("ref_var_main"), label = "Reference variable",
                choices = list(var1_id(), var2_id()),
                selected = var1_id()
            )
        } else if (nchar(var2_id()) < 1) {
            radioButtons(
                inputId = ns("ref_var_main"),
                label = "Reference variable",
                choices = list(var1_id()),
                selected = var1_id()
            )
        } else {
            radioButtons(
                inputId = ns("ref_var_main"),
                label = "Reference variable",
                choices = list(var1_id(), var2_id()),
                selected = var1_id()
            )
        }
    })

    # filtered data for downloading (Main Profile ) ----------------------------
    download_data <- reactive({
        if (is.null(data())) return()
        ### filtered data
        data_out <- data()

        data_out <- as.data.frame(data_out[data_out$presSpec > 0, ])
        data_out <- data_out[!is.na(data_out$geneID), ]

        data_out <- as.data.frame(data_out[data_out$presSpec >= percent()[1], ])
        if (length(var1()) == 1) {
            data_out <- as.data.frame(data_out[data_out$var1 >= var1()[1], ])
        } else {
            data_out <- as.data.frame(data_out[data_out$var1 >= var1()[1]
                                               & data_out$var1 <= var1()[2], ])
        }

        if (!all(is.na(data_out$var2))) {
            if (length(var2()) == 1) {
                data_out <-
                    as.data.frame(data_out[data_out$var2 >= var2()[1], ])
            } else {
                data_out <-
                    as.data.frame(data_out[data_out$var2 >= var2()[1]
                                           & data_out$var2 <= var2()[2], ])
            }
        } else {
            data_out$var2 <- 0
        }

        ### select only representative genes if chosen
        if (input$get_representative_main == TRUE) {
            if (is.null(input$ref_var_main)) return()
            else {
                if (input$ref_var_main == var1_id()) {
                    data_out_agg <- aggregate(
                        as.numeric(data_out$var1),
                        by = list(data_out$geneID, data_out$ncbiID),
                        FUN = input$ref_type_main
                    )
                } else if (input$ref_var_main == var2_id()) {
                    data_out_agg <- aggregate(
                        as.numeric(data_out$var2),
                        by = list(data_out$geneID, data_out$ncbiID),
                        FUN = input$ref_type_main
                    )
                } else {
                    data_out_agg <-
                        data_out[data_out, c("geneID", "ncbiID", "var1")]
                }
                colnames(data_out_agg) <- c("geneID", "ncbiID", "var_best")

                data_out_representative <- merge(data_out, data_out_agg,
                                                 by = c("geneID", "ncbiID"),
                                                 all.x = TRUE)

                if (input$ref_var_main == var1_id()) {
                    data_out <-
                        data_out_representative[
                            data_out_representative$var1 ==
                                data_out_representative$var_best,
                            ]
                } else if (input$ref_var_main == var2_id()) {
                    data_out <-
                        data_out_representative[
                            data_out_representative$var2 ==
                                data_out_representative$var_best,
                            ]
                } else {
                    data_out <- data_out
                }
                # used to select only one ortholog,
                # if there exist more than one "representative"
                data_out$dup <- paste0(data_out$geneID, "#", data_out$ncbiID)
                data_out <- data_out[!duplicated(c(data_out$dup)), ]
            }
        }

        # sub select columns of dataout
        data_out <- data_out[, c("geneID",
                                 "orthoID",
                                 "fullName",
                                 "ncbiID",
                                 "supertaxon",
                                 "var1",
                                 "var2",
                                 "presSpec")]
        data_out <- data_out[order(data_out$geneID, data_out$supertaxon), ]
        data_out <- data_out[complete.cases(data_out), ]

        data_out$geneID <- as.character(data_out$geneID)
        data_out$fullName <- as.character(data_out$fullName)
        data_out$ncbiID <- as.character(data_out$ncbiID)
        data_out$supertaxon <- substr(data_out$supertaxon,
                                      6,
                                      nchar(as.character(data_out$supertaxon)))
        data_out$var1 <- as.character(data_out$var1)
        data_out$var2 <- as.character(data_out$var2)
        data_out$presSpec <- as.numeric(data_out$presSpec)

        # rename columns
        names(data_out)[names(data_out) == "presSpec"] <- "%Spec"
        if (nchar(var1_id()) > 0) {
            names(data_out)[names(data_out) == "var1"] <- var1_id()
        } else {
            data_out <- subset(data_out, select = -c(var1) )
        }
        if (nchar(var2_id()) > 0) {
            names(data_out)[names(data_out) == "var2"] <- var2_id()
        } else {
            data_out <- subset(data_out, select = -c(var2) )
        }

        # return data for downloading
        data_out <- as.matrix(data_out)
        return(data_out)
    })

    # download data ------------------------------------------------------------
    output$download_data <- downloadHandler(
        filename = function(){
            c("filteredData.out")
        },
        content = function(file){
            data_out <- download_data()
            write.table(data_out, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # render download data table -----------------------------------------------
    output$filtered_main_data <- renderDataTable(rownames = FALSE, {
        data <- download_data()
        data
    })

    # download FASTA -----------------------------------------------------------
    output$download_fasta <- downloadHandler(
        filename = function(){
            c("filteredSeq.fa")
        },
        content = function(file){
            fasta_out_df <- fasta()
            write.table(fasta_out_df, file,
                        sep = "\t",
                        col.names = FALSE,
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # download data as long format ---------------------------------------------
    download_data_long <- reactive({
        download_data <- download_data()

        if (ncol(download_data) == 6) {
            download_data_long <- download_data[,c(1,4,2)]
        } else if (ncol(download_data) == 7) {
            download_data_long <- download_data[,c(1,4,2,6)]
        } else if (ncol(download_data) == 8) {
            download_data_long <- download_data[,c(1,4,2,6,7)]
        }

        return(download_data_long)
    })

    output$download_long <- downloadHandler(
        filename = function(){
            c("filteredData.phyloprofile")
        },
        content = function(file){
            data_out <- download_data_long()
            write.table(data_out, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(download_data)
}
