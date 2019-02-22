#' Core gene identification
#'
#' @export
#' @param filtered_data full processed main data
#' (from reactive fn "get_data_filtered")
#' @param rank_select selected taxonomy rank (input$rank_select)
#' @param taxa_core selected list of taxa (input$taxa_core)
#' @param percent_core cutoff of percentage taxa present in a supertaxon
#' (input$percent_core)
#' @return list of core genes
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

# source("R/functions.R")

identify_core_gene_ui <- function(id){
    ns <- NS(id)
    dataTableOutput(ns("core_gene.table"))
}

identify_core_gene <- function(input, output, session,
                               filtered_data,
                               rank_select, taxa_core, percent_core,
                               var1_cutoff, var2_cutoff,
                               core_coverage){

    output$core_gene.table <- renderDataTable({
        data <- core_geneDf()
        if (is.null(data)) return()
        else {
            data <- as.data.frame(data)
            data
        }
    })

    core_geneDf <- reactive({
        core_geneDf <- get_core_gene(rank_select(),
                                     taxa_core(),
                                     filtered_data(),
                                     var1_cutoff(), var2_cutoff(),
                                     percent_core(), core_coverage())
        return(core_geneDf)
    })

    if (!is.null(core_geneDf)) return(core_geneDf)
}
