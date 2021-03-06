#' Download filtered data from customized profile
#'
#' @export
#' @param data MAIN data for downloading (from module "download_filtered_main.R"
#' @param fasta fasta sequences (from reactive fn "customized_fasta_download")
#' @param in_seq selected sequences in customized profile (from input$in_seq)
#' @param in_taxa selected taxa in customized profile (from input$in_taxa)
#' @return data of customized profile for downloading
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

download_filtered_customized_ui <- function(id) {
    ns <- NS(id)

    tabPanel(
        "Customized data",
        column(
            12,
            strong(
                em(
                    "NOTE: Depending on your choice in Download filtered data
                    -> Main data, either all or only representative sequences
                    will be downloaded!"
                ),
                style = "color:red"
            ),
            hr()
        ),
        column(
            12,
            dataTableOutput(ns("filtered_custom_data"))
        ),
        column(
            4,
            downloadButton(ns("download_custom_data"),
                           "Download customized data")
        ),
        column(
            4,
            downloadButton(ns("download_custom_fasta"),
                           "Download FASTA sequences"),
            uiOutput(ns("download_custom_fasta.ui"))
        ),
        column(
            4,
            downloadButton(ns("download_custom_long"),
                           "Download data as PhyloProfile input format")
        )
    )
}

download_filtered_customized <- function(input, output, session,
                                         data,
                                         fasta,
                                         in_seq, in_taxa){

    # filtered data for downloading (Customized Profile) -----------------------
    download_custom_data <- reactive({
        data <- as.data.frame(data())

        # get subset of data according to selected genes/taxa
        if (!is.null(in_seq()) | !is.null(in_taxa())) {
            if (in_seq()[1] != "all" & in_taxa()[1] == "all") {
                # select data for selected sequences only
                custom_data <- subset(data, geneID %in% in_seq())
            } else if (in_seq()[1] == "all" & in_taxa()[1] != "all") {
                # select data for selected taxa only
                custom_data <- subset(data, supertaxon %in% in_taxa())
            } else if (in_seq()[1] != "all" & in_taxa()[1] != "all") {
                # select data for selected sequences and taxa
                custom_data <- subset(data, geneID %in% in_seq()
                                      & supertaxon %in% in_taxa())
            } else {
                custom_data <- data
            }
        } else {
            custom_data <- data
        }
        # return data
        custom_data <- as.matrix(custom_data)
        return(custom_data)
    })

    # download data ------------------------------------------------------------
    output$download_custom_data <- downloadHandler(
        filename = function(){
            c("customFilteredData.out")
        },
        content = function(file){
            data_out <- download_custom_data()
            write.table(data_out, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # render download data table -----------------------------------------------
    output$filtered_custom_data <- renderDataTable(rownames = FALSE, {
        data <- download_custom_data()
        data
    })

    # download FASTA -----------------------------------------------------------
    output$download_custom_fasta <- downloadHandler(
        filename = function(){
            c("customFilteredSeq.fa")
        },
        content = function(file){
            fasta_out_df <- fasta()
            write.table(fasta_out_df,
                        file,
                        sep = "\t",
                        col.names = FALSE,
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # download data as long format ---------------------------------------------
    download_custom_data_long <- reactive({
        download_custom_data <- download_custom_data()

        if (ncol(download_custom_data) == 6) {
            download_custom_data_long <- download_custom_data[,c(1,4,2)]
        } else if (ncol(download_custom_data) == 7) {
            download_custom_data_long <- download_custom_data[,c(1,4,2,6)]
        } else if (ncol(download_custom_data) == 8) {
            download_custom_data_long <- download_custom_data[,c(1,4,2,6,7)]
        }

        return(download_custom_data_long)
    })

    output$download_custom_long <- downloadHandler(
        filename = function(){
            c("customFilteredData.phyloprofile")
        },
        content = function(file){
            data_out <- download_custom_data_long()
            write.table(data_out, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(download_custom_data)
}
