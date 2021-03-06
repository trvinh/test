#' Group Comparison
#'
#' @export
#' @param selected_in_group selected In-group (input$selected_in_group_gc)
#' @param selected_genes_list list of genes to calculate the plots for
#'  (input$list_selected_genes_gc)
#' @param main_rank in input and settings selected taxonomy rank
#'  (input$rank_select)
#' @param selected_variable variable(s) to claculate the plots for
#'  (input$var_name_gc)
#' @param use_common_ancestor boolean if the next common anchestor should be
#'  used (input$use_common_ancestor)
#' @param reference_taxon selected taxon name (input$in_select)
#' @param ncbi_id_list list of ncbi ids (from reactive fn "input_taxonID")
#' @param filtered_data full processed main data
#' (from reactive fn "get_data_filtered")
#' @param right_format_features boolean if the features have the right format
#' (input$right_format_features)
#' @param domain_information dataframe holding the domain input
#' (from reactive fn "get_domain_information)
#' @param plot information if the plots should be generated (input$plot_gc)
#' @param parameter list of parameters needed to generate the plots
#' (from reactive fn "get_parameter_input_gc")
#' @return list of candidate genes
#' @author Carla Mölbert {carla.moelbert@gmx.de}

source("R/functions.R")

group_comparison_ui <- function(id){
    ns <- NS(id)
    fluidPage(
        sidebarPanel(
            withSpinner(uiOutput(ns("candidate_genes"))),
            bsPopover(
                "candidate_genes",
                "",
                "Select gene to show the plots",
                "right"
            ),
            withSpinner(uiOutput(ns("features_of_interest_ui"))),
            bsPopover(
                "features_of_interest_ui",
                "",
                "This function is only use full if the features are
        saved in the right format: featuretype_featurename"
            ),
            sliderInput(
                ns("domains_threshold"),
                "Persentage of proteins in which the domain needs to have an
                instance",
                min = 0 , max = 100, value = 0,
                step = 1, round = FALSE
            ),

            downloadButton(ns("download_plots"), "Download plots"),
            width = 3
        ),
        mainPanel(
            tags$style(
                HTML("#plots_ui { height:650px; overflow-y:scroll}")
            ),
            withSpinner(uiOutput(ns("plots_ui"))),
            width = 9
        )
    )
}

group_comparison <- function(input, output, session,
                             selected_in_group,
                             selected_genes_list,
                             main_rank,
                             selected_variable,
                             use_common_ancestor,
                             reference_taxon,
                             ncbi_id_list,
                             filtered_data,
                             right_format_features,
                             domain_information,
                             plot,
                             parameter,
                             selected_point){
    # Dataframe for the significant Genes ======================================
    #' contains geneID, in_group, out_group, pvalues, features, databases,
    #' rank, var
    candidate_genes <- reactiveValues(plots = NULL)

    # List with all candidate genes ============================================
    output$candidate_genes <- renderUI({
        ns <- session$ns
        plot()

        isolate({
            candidate_genes$plots <- {
                get_significant_genes(selected_in_group(),
                                      selected_genes_list(),
                                      main_rank(),
                                      selected_variable(),
                                      use_common_ancestor(),
                                      reference_taxon(),
                                      parameter(),
                                      ncbi_id_list(),
                                      filtered_data(),
                                      right_format_features(),
                                      domain_information())
            }

            if (is.data.frame(candidate_genes$plots)) {
                significant_genes <- candidate_genes$plots
                x <- as.vector(significant_genes$geneID)
                choices <- c("all", x)

                selectInput(ns("selected_gene"), "Candidate gene(s):",
                            choices,
                            selected = choices[2],
                            multiple = FALSE)

            } else{
                selectInput(ns("selected_gene"), "Candidate gene(s):",
                            NULL,
                            selected = NULL,
                            multiple = FALSE)
            }
        })
    })

    # Output of the plots for the selected gene(s) =============================
    output$plots_ui <- renderUI({
        if (is.character(candidate_genes$plots)) return(candidate_genes$plots)
        get_plots()
    })

    # List with possible features for the selected gene ========================
    output$features_of_interest_ui <- renderUI({
        ns <- session$ns
        input$selected_gene
        isolate({
            gene <- input$selected_gene
            if (!right_format_features()) {
                selectInput(ns("interesting_features"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else if (is.null(gene)) {
                selectInput(ns("interesting_features"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else if (gene == "") {
                selectInput(ns("interesting_features"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else{
                significant_genes <- candidate_genes$plots
                choices <- c("all")
                if (gene == "all") {
                    for (current_gene in significant_genes$geneID) {
                        subset_current_gene <-
                            subset(significant_genes,
                                   significant_genes$geneID == current_gene)
                        choices <- append(
                            choices, unlist(subset_current_gene$databases)
                        )
                    }
                    # show each database only once
                    choices <- choices[!duplicated(choices)]
                } else {

                    subset_gene <- subset(significant_genes,
                                          significant_genes$geneID == gene)

                    choices <- append(choices, unlist(subset_gene$databases))

                }
                selectInput(ns("interesting_features"),
                            "Feature type(s) of interest:",
                            choices,
                            selected = choices[1],
                            multiple = TRUE,
                            selectize = FALSE)
            }
        })
    })

    # download file with the shown plots =======================================
    output$download_plots <- downloadHandler(
        filename = "plotSignificantGenes.zip",
        content = function(file){
            genes <- input$selected_gene
            significant_genes <- candidate_genes$plots
            if ("all" %in% genes) {
                genes <- significant_genes$geneID
            }

            fs <- c()
            #tmpdir <- tempdir()
            setwd(tempdir())

            for (gene in genes) {
                path <- paste(gene, ".pdf", sep = "")
                fs <- c(fs, path)
                pdf(path)
                get_plots_to_download(gene,
                                      parameter(),
                                      input$interesting_features,
                                      significant_genes,input$domains_threshold,
                                      selected_point())
                dev.off()
            }
            zip(zipfile = file, files = fs)
        },
        contentType = "application/zip"
    )

    #' observer for the download functions
    observe({
        if (is.null(selected_in_group())
            | length(selected_genes_list()) == 0) {
            shinyjs::disable("download_plots")
        } else if (plot() == FALSE) {
            shinyjs::disable("download_plots")
        } else if (input$selected_gene == "") {
            shinyjs::disable("download_plots")
        } else {
            shinyjs::enable("download_plots")
        }
    })

    # Deciding which plots will be shown =======================================
    get_plots <- reactive({
        input$interesting_features
        input$domains_threshold
        gene <- as.character(input$selected_gene)
        plot()
        if (is.null(candidate_genes$plots)) return()
        significant_genes <- candidate_genes$plots
        if (gene == "all") {
            plot_output_list <- get_plot_output_list(significant_genes,
                                                     parameter(),
                                                     input$interesting_features,
                                                     input$domains_threshold,
                                                     selected_point())
        }else{
            gene_info <- {
                significant_genes[significant_genes$geneID == gene, ]
            }
            if (nrow(gene_info) == 0) return()
            plot_output_list <- get_plot_output_list(gene_info,
                                                     parameter(),
                                                     input$interesting_features,
                                                     input$domains_threshold,
                                                     selected_point())
        }
        #' List with all plots that will be shown
        return(plot_output_list)
    })

    # List of genes for the customized profile =================================
    gene_list <- reactive({
        if (!is.null(candidate_genes$plots)) {
            significant_genes <- candidate_genes$plots
            return(significant_genes$geneID)
        }
    })
    return(gene_list)
}

# FUNCTIONS ===================================================================

#' Generate the list with all plots -------------------------------------------
#' @export
#' @param genes list of genes
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases to take the features from
#' @return list with all plots
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_plot_output_list <- function(genes, parameters, interesting_features,
                                 domains_threshold, selected_point) {
    if (is.null(genes)) return()
    # Insert plot output objects the list
    plot_output_list <- lapply(1:nrow(genes), function(i) {
        plotname <- paste(genes[i, 1])
        plot_output_object <- renderPlot(get_multiplot(genes[i, ], parameters,
                                                       interesting_features,
                                                       domains_threshold,
                                                       selected_point),
                                         height = 650, width = 700)
    })
    do.call(tagList, plot_output_list) # needed to display properly.
    return(plot_output_list)
}

#' Put the plots for one spicific gene in one multiplot -----------------------
#' @export
#' @param gene_info contains "geneID",  "in_group",  "out_group", "pvalues",
#' "features",  "databases", "var", "rank"
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases to take the features from
#' @return grid arrange with the plots that should be shown for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_multiplot <- function(gene_info, parameters, interesting_features,
                          domains_threshold,
                          selected_point){
    #' Sorting the information to the selected gene ----------------------------
    gene <- as.character(gene_info$geneID)
    in_group <- as.data.frame(gene_info$in_group)
    out_group <- as.data.frame(gene_info$out_group)
    features <- as.data.frame(gene_info$features)

    var <- gene_info$var

    #' Get the barplot ---------------------------------------------------------
    barplot <-  get_barplot_gc(gene,
                               in_group,
                               out_group,
                               features,
                               parameters,
                               interesting_features,
                               domains_threshold)

    if (is.null(barplot)) {
        barplot <- textGrob("The selected domains are not found in the gene")
    }

    #' Get the boxplots  for two variables  ------------------------------------
    if (var == "Both") {
        p_value1 <- round(gene_info$pvalues_1, 2)
        p_value2 <- round(gene_info$pvalues_2, 2)

        #' Check if the p_values should be printed
        if (parameters$show_p_value == TRUE) {
            info_p_value1 <- paste("P-value:", p_value1, sep = " ")
            info_p_value2 <- paste("P-value:", p_value2, sep = " ")
        } else {
            info_p_value1 <- " "
            info_p_value2 <- " "
        }

        #' Get information about the plot colour
        if (parameters$highlight_significant == TRUE) {
            if (is.na(p_value1)) colour1 <- "grey"
            else if (p_value1 < parameters$significance) {
                colour1 <- "indianred2"
            } else colour1 <- "grey"

            if (is.na(p_value2)) colour2 <- "grey"
            else if (p_value2 < parameters$significance) {
                colour2 <- "indianred2"
            } else colour2 <- "grey"
        } else {
            colour1 <- "grey"
            colour2 <- "grey"
        }

        #' Generate the boxplots
        boxplot1 <- get_boxplot_gc(in_group,
                                   out_group,
                                   parameters$var1_id,
                                   gene,
                                   colour1,
                                   info_p_value1, parameters,
                                   selected_point)

        boxplot2 <- get_boxplot_gc(in_group,
                                   out_group,
                                   parameters$var2_id,
                                   gene,
                                   colour2,
                                   info_p_value2, parameters,
                                   selected_point)

        plots <- grid.arrange(textGrob(gene),
                              arrangeGrob(boxplot1, boxplot2, ncol = 2),
                              barplot,
                              heights = c(0.02, 0.45, 0.458), ncol = 1)
    }
    #' get the boxplot if one varibale is selected  ----------------------------
    else {
        p <- round(gene_info$pvalue, 2)

        #' Check if the p_values should be printed
        if (parameters$show_p_value == TRUE) {
            info <- paste("P-value:", p)
        } else {
            info <- " "
        }

        #' Generate the plot
        boxplot <- get_boxplot_gc(in_group,
                                  out_group,
                                  var,
                                  gene,
                                  "grey",
                                  info,
                                  parameters,
                                  selected_point)

        plots <- grid.arrange(textGrob(gene),
                              boxplot,
                              barplot,
                              heights = c(0.02, 0.45, 0.458), ncol = 1)
    }

    #' return the plots --------------------------------------------------------
    return(plots)
}

#' Create a Boxplot -----------------------------------------------------------
#' @export
#' @param in_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param out_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param var variable to consider in the boxplot
#' @param  gene gene for which the plot is generated
#' @param  colour colour of the boxes
#' @param  info info about the p-values
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return boxplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_boxplot_gc <- function(in_group_df,
                           out_group_df,
                           var,
                           gene,
                           colour,
                           info,
                           parameters,
                           selected_point){

    #' pre-processing the data for the boxplot ---------------------------------
    if (var == parameters$var1_id) {
        in_group <- in_group_df$var1
        out_group <- out_group_df$var1
    } else if (var == parameters$var2_id) {
        in_group <- in_group_df$var2
        out_group <- out_group_df$var2
    }

    in_group <- in_group[!is.na(in_group)]
    out_group <- out_group[!is.na(out_group)]

    length_in_group <- length(in_group)
    length_out_group <- length(out_group)

    in_group <- as.data.frame(in_group)
    names(in_group)[1] <- paste("values")
    in_group$group <- "in_group"

    out_group <- as.data.frame(out_group)
    names(out_group)[1] <- paste("values")
    out_group$group <- "Out-Group"

    data_boxplot <- rbind(in_group, out_group)
    data_boxplot <- data_boxplot[complete.cases(data_boxplot), ]

    names <- c(paste("In-Group \n n=", length_in_group, sep = ""),
               paste("Out-Group \n n=", length_out_group, sep = ""))

    #' Generate the boxplot ----------------------------------------------------
    boxplot_gc <- ggplot(data_boxplot, aes(group, values)) +
        geom_violin(position = position_dodge(), scale = "width",
                    fill = colour) +
        labs(x = "", y = var, caption = paste(info), colour = "") +
        scale_x_discrete(labels = names) +
        theme_minimal() +
        stat_summary(aes(colour = selected_point),fun.y = selected_point,
                     geom = "point", size = 3, show.legend = TRUE)

    boxplot_gc <- boxplot_gc +
        theme(axis.text.x = element_text(size = parameters$x_size_gc,
                                         hjust = 1),
              axis.text.y = element_text(size = parameters$y_size_gc),
              axis.title.y = element_text(size = parameters$y_size_gc),
              plot.caption = element_text(size = parameters$p_values_size),
              legend.position = "bottom",
              legend.text = element_text(size = parameters$legend_size_gc ),
              legend.title = element_text(size = parameters$legend_size_gc)) +
        scale_color_manual("", values = c("green"))

    #' return the boxplot ------------------------------------------------------
    return(boxplot_gc)
}


#' Create a Barplot -----------------------------------------------------------
#' @export
#' @param  selected_gene gene for which the plot is generated
#' @param in_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param out_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param features contains "seedID",  "orthoID", "feature", "start",   "end"
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases for which the features should
#' be included
#' @return barplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_barplot_gc <- function(selected_gene,
                           in_group, out_group,
                           features, parameters,
                           interesting_features,
                           domains_threshold){

    #' information about the features ------------------------------------------
    features$feature <- as.character(features$feature)

    #' part in in_group and out-group  -----------------------------------------
    in_group$orthoID <- gsub("\\|", ":", in_group$orthoID)
    in_group_domain_df  <-  {
        subset(features, features$orthoID %in% in_group$orthoID)
    }
    out_group$orthoID <- gsub("\\|", ":", out_group$orthoID)
    out_group_domain_df <- {
        subset(features, features$orthoID %in% out_group$orthoID)
    }

    #' get the dataframes for in- and out-group --------------------------------
    if (nrow(in_group_domain_df) == 0) {
        data_in <- NULL
    } else {
        feature <- unique(in_group_domain_df$feature)
        data_in <- as.data.frame(feature)
        data_in$amount <- 0
        data_in$proteins <- 0
    }

    if (nrow(out_group_domain_df) == 0) {
        data_out <- NULL
    } else {
        feature <- unique(out_group_domain_df$feature)
        data_out <- as.data.frame(feature)
        data_out$amount <- 0
        data_out$proteins <- 0
    }

    seeds <- features$seedID
    seeds <- unique(seeds)

    #' Get the values for the boxplot ------------------ -----------------------
    in_not_empty <- 0
    out_not_empty <- 0

    #' Count for each feature how often it is present in each seed -------------
    for (seed in seeds) {
        #' count the features in the in-group
        if (!is.null(data_in)) {
            in_g <- subset(in_group_domain_df,
                           in_group_domain_df$seedID == seed)
            in_g <- in_g[str_detect(in_g$seedID, in_g$orthoID),]

            if (!empty(in_g)) {
                in_not_empty <- in_not_empty + 1
                in_group_features <-  plyr::count(in_g, "feature")
                for (i in 1:nrow(in_group_features)) {
                    for (j in 1:nrow(data_in)) {
                        if (data_in[j, 1] == in_group_features[i, 1]) {
                            data_in[j, 2] <-
                                data_in[j, 2] + in_group_features[i, 2]
                            data_in[j, 3] <- data_in[j,3] + 1
                        }
                    }
                }
            }
        }

        #' count the featueres in the out-group
        if (!is.null(data_out)) {
            out_g <- subset(out_group_domain_df,
                            out_group_domain_df$seedID == seed)
            out_g <- out_g[str_detect(out_g$seedID, out_g$orthoID),]

            if (!empty(out_g)) {
                out_not_empty <- out_not_empty + 1
                out_group_features <-  plyr::count(out_g, "feature")
                for (i in 1:nrow(out_group_features)) {
                    for (j in 1:nrow(data_out)) {
                        if (data_out[j, 1] == out_group_features[i, 1]) {
                            data_out[j, 2] <-
                                data_out[j, 2] + out_group_features[i, 2]
                            data_out[j, 3] <- data_out[j,3] + 1
                        }
                    }
                }
            }
        }
    }

    if (domains_threshold > 0) {
        threshold_in  <- as.numeric((domains_threshold / 100) * in_not_empty)
        threshold_out <- as.numeric((domains_threshold / 100) * out_not_empty)

        data_in <- data_in[sort(data_in$feature),]
        data_out <- data_out[sort(data_out$feature),]

        features_in <- data_in$feature[!(data_in$feature %in% data_out$feature)]
        features_out <-
            data_out$feature[!(data_out$feature %in% data_in$feature)]

        for (i in features_out) {
            data_in <- rbind(data_in, c(i, 0, 0))
        }

        for (i in features_in) {
            data_out <- rbind(data_out, c(i, 0, 0))
        }

        data_in$proteins <- as.numeric(data_in$proteins)
        data_out$proteins <- as.numeric(data_out$proteins)

        data_in_keep <- data_in[data_in$proteins >= threshold_in,]
        data_out_keep <- data_out[data_out$proteins >= threshold_out, ]

        keep <- rbind(data_in_keep, data_out_keep)
        keep <- keep[!duplicated(keep$feature),]

        data_in <- data_in[data_in$feature %in% keep$feature,]
        data_out <- data_out[data_out$feature %in% keep$feature,]
    }

    #' Calculate the average of appearances for each feature -------------------
    if (!is.null(data_in)) {
        data_in$amount <- as.numeric(data_in$amount)
        data_in$amount <- data_in$amount / in_not_empty
        data_in$type <- "In-Group"
    }

    if (!is.null(data_out)) {
        data_out$amount <- as.numeric(data_out$amount)
        data_out$amount <- data_out$amount / out_not_empty
        data_out$type <- "Out-Group"
    }

    #' Get the data for teh barplot --------------------------------------------
    if (is.null(data_in) & !is.null(data_out)) {
        data_barplot <- data_out
    }  else if (is.null(data_out) & !is.null(data_in)) {
        data_barplot <- data_in
    } else if (!is.null(data_in) & !is.null(data_out)) {
        data_barplot <- rbind(data_in, data_out)
    } else {
        data_barplot <- NULL
    }

    #' only show features that interest the user
    if (!("all" %in% interesting_features)) {
        features_list <- NULL
        for (feature in interesting_features) {
            data_barplot$feature <- as.character(data_barplot$feature)
            subset_features <- subset(data_barplot$feature,
                                      startsWith(data_barplot$feature, feature))
            features_list <- append(features_list, subset_features)
        }

        if (is.null(features_list)) return()
        #' only keep rows in which the feature begins with a element out of the
        #' interesing Features
        features_list <- features_list[!duplicated(features_list)]
        data_barplot <- subset(data_barplot, data_barplot$feature
                               %in% features_list)
    }

    data_barplot <- data_barplot[sort(data_barplot$feature),]

    #' generate the barplot ----------------------------------------------------
    if (!is.null(data_barplot)) {
        barplot_gc <- ggplot(data_barplot,
                             aes(x = feature, y = amount, fill = type ),
                             main  = " ") +
            geom_bar(
                stat = "identity", position = position_dodge(), width = 0.5
            ) +
            scale_fill_grey() +
            labs(x = " ", y = "Average instances per protein", fill = "Group") +
            theme_minimal()

        barplot_gc <- barplot_gc +
            theme(axis.text.x = element_text(size = parameters$x_size_gc,
                                             angle = parameters$angle_gc,
                                             hjust = 1),
                  axis.text.y = element_text(size = parameters$y_size_gc),
                  axis.title.y = element_text(size = parameters$y_size_gc),
                  legend.position = parameters$legend_gc,
                  legend.text = element_text(size = parameters$legend_size_gc ),
                  legend.title = element_text(size = parameters$legend_size_gc))

        #' return the barplot --------------------------------------------------
        return(barplot_gc)
    } else (return(NULL))
}

#' get the plots to download --------------------------------------------------
#' @export
#' @param  gene gene for which the plot is generated
#' @param interesting_features list of databases for which the features should
#' be included
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return arrange grop containing the plots for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_plots_to_download <- function(gene,
                                  parameters,
                                  interesting_features,
                                  significant_genes,
                                  domain_threshold,
                                  selected_point){
    info_gene <- subset(significant_genes,
                        significant_genes$geneID == gene)
    return(get_multiplot(info_gene, parameters, interesting_features,
                         domains_threshold,
                         selected_point))
}

#' Get the p_values to print under the plot -----------------------------------
#' @export
#' @param pvalues list contianing the p-values for a specific gene
#' @return string containing the information about the p-values
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_info_p_values <- function(pvalues) {

    if (is.na(pvalues[1])) info_p_values <- "not enough information"
    else if (length(pvalues) == 1) {
        info_p_values <- paste("Kolmogorov-Smirnov-Test:",
                               as.character(pvalues[1]), sep = " " )
    } else{
        info_p_values1 <- paste("Kolmogorov-Smirnov-Test:",
                                as.character(pvalues[1]), sep = " " )
        info_p_values2 <- paste("Wilcoxon-Mann-Whitney-Test: ",
                                as.character(round(pvalues[2], 6)), sep = " ")
        info_p_values <- paste(info_p_values1, info_p_values2, sep = "\n")
    }

    info_p_values <- paste("p_values:", info_p_values, sep = "\n")
    return(info_p_values)
}
