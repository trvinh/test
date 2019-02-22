#' Profile plot
#'
#' @export
#' @param data data for heatmap plot (from reactive fn "dataHeat")
#' @param clusteredDataHeat clustered data (from reactive fn "clusteredDataHeat"
#' @param apply_cluster choose clustered data or not (from input$apply_cluster)
#' @param parameters plot parameters (colors, size, variable names, ...)
#' @param in_seq subset sequences for customized profile (input$in_seq)
#' @param in_taxa subset taxa for customized profile (input$in_taxa)
#' @param rank_select selected taxonomy rank (input$rank_select)
#' @param in_select selected taxon name (input$in_select)
#' @param taxon_highlight highlighted taxon (input$taxon_highlight)
#' @param gene_highlight highlighted gene (input$gene_highlight)
#' @param type_profile either "main_profile" or "customized_profile"
#' @return info for selected point on the profile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

create_profile_plot_ui <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("plot.ui")),
        br(),
        downloadButton(ns("profile_download"),"Download profile",
                       class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        )
    )
}

create_profile_plot <- function(input, output, session,
                                data, clusteredDataHeat,
                                apply_cluster,
                                parameters,
                                in_seq, in_taxa,
                                rank_select, in_select,
                                taxon_highlight, gene_highlight,
                                type_profile,
                                color_by_group) {
    # data for heatmap ---------------------------------------------------------
    dataHeat <- reactive({
        if (is.null(data())) return()

        if (type_profile() == "customized_profile") {
            if (is.null(in_taxa()) | is.null(in_seq())) return()

            data_heat <- data_customized_plot(data(), in_taxa(), in_seq())
            if (apply_cluster() == TRUE) {
                data_heat <- data_customized_plot(clusteredDataHeat(),
                                                  in_taxa(), in_seq())
            }
        } else {
            data_heat <- data_main_plot(data())
            if (apply_cluster() == TRUE) {
                data_heat <- data_main_plot(clusteredDataHeat())
            }
        }
        return(data_heat)
    })

    # render heatmap profile ---------------------------------------------------
    output$plot <- renderPlot({
        if (is.null(data())) return()
        if (type_profile() == "customized_profile") {
            if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
        }

        highlight_profile_plot(
            dataHeat(),
            parameters(),
            taxon_highlight(),
            rank_select(),
            gene_highlight()
        )
    })

    output$plot.ui <- renderUI({
        ns <- session$ns

        if (type_profile() == "customized_profile") {
            if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
            else if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
        }

        withSpinner(
            plotOutput(
                ns("plot"),
                width = parameters()$width,
                height = parameters()$height,
                click = ns("plot_click")
            )
        )
    })

    output$profile_download <- downloadHandler(
        filename = function() {
            c("profile.pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = highlight_profile_plot(dataHeat(),
                                          parameters(),
                                          "none",
                                          rank_select(),
                                          "none"),

                width = parameters()$width * 0.056458333,
                height = parameters()$height * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )
    # get info of clicked point on heatmap plot --------------------------------
    selectedpoint_info <- reactive({

        # get selected supertaxon name
        taxa_list <- phyloprofile::get_name_list()
        # rank_select <- rank_select()
        # rankName <- substr(rank_select, 4, nchar(rank_select))
        rankName <- rank_select()
        in_select <- {
            as.numeric(taxa_list$ncbiID[taxa_list$fullName == in_select()])
        }

        dataHeat <- dataHeat()
        if (is.null(dataHeat)) return()

        if (type_profile() == "customized_profile") {
            # get sub-dataframe of selected taxa and sequences
            dataHeat$supertaxonMod <- substr(
                dataHeat$supertaxon,
                6,
                nchar(as.character(dataHeat$supertaxon))
            )

            if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
            if (in_taxa()[1] == "all" & in_seq()[1] != "all") {
                # select data from dataHeat for selected sequences only
                dataHeat <- subset(dataHeat, geneID %in% in_seq())
            } else if (in_seq()[1] == "all" & in_taxa()[1] != "all") {
                # select data from dataHeat for selected taxa only
                dataHeat <- subset(dataHeat, supertaxonMod %in% in_taxa())
            } else {
                # select data from dataHeat for selected sequences and taxa
                dataHeat <- subset(dataHeat, geneID %in% in_seq()
                                   & supertaxonMod %in% in_taxa())
            }

            # drop all other supertaxon that are not in sub-dataframe
            dataHeat$supertaxon <- factor(dataHeat$supertaxon)
            dataHeat$geneID <- factor(dataHeat$geneID)
        }

        # get values
        if (is.null(input$plot_click$x)) return()
        else {
            # get cooridiate point
            if (parameters()$x_axis == "genes") {
                corX <- round(input$plot_click$y);
                corY <- round(input$plot_click$x)
            } else {
                corX <- round(input$plot_click$x);
                corY <- round(input$plot_click$y)
            }

            # get geneID
            genes <- levels(dataHeat$geneID)
            geneID <- toString(genes[corY])
            # get supertaxon (spec)
            supertaxa <- levels(dataHeat$supertaxon)
            spec <- toString(supertaxa[corX])
            # get var1, percentage of present species and var2 score
            var1 <- NA
            if (!is.na(dataHeat$var1[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var1 <- max(
                    na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec])
                )
            }
            Percent <- NA
            if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                         & dataHeat$supertaxon == spec][1])) {
                Percent <- {
                    max(
                        na.omit(
                            dataHeat$presSpec[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]
                        )
                    )
                }
            }
            var2 <- NA
            if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var2 <- {
                    max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]))
                }
            }

            # get ortholog ID
            orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]
            if (length(orthoID) > 1) {
                orthoID <- paste0(orthoID[1], ",...")
            }

            if (is.na(as.numeric(Percent))) return()
            else {
                info <- c(geneID,
                          as.character(orthoID),
                          as.character(spec),
                          round(as.numeric(var1), 2),
                          round(as.numeric(Percent), 2),
                          round(as.numeric(var2), 2))
                return(info)
            }
        }
    })

    return(selectedpoint_info)
}
