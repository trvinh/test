#' Profile clustering
#'
#' @param distance_matrix
#' @param cluster_method Method to cluster the distances (input$cluster_method)
#' @param plot_width Width of the generated plot (input$cluster_plot.width)
#' @param plot_height Height of the generated plot (input$cluster_plot.height)

source("R/functions.R")

cluster_profile_ui <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            8,
            downloadButton("download_cluster",
                           "Download plot", class = "butDL"),
            tags$head(tags$style(
                HTML(".butDL{background-color:#476ba3;} .butDL{color: white;}")
            )),
            uiOutput(ns("cluster.ui"))
        ),
        column(
            4,
            downloadButton(ns("download_distance_matrix"),
                           "Download distance matrix"),
            downloadButton(ns("download_cluster_genes"),
                           "Download gene list"),
            tableOutput(ns("brushed_cluster.table"))
        )
    )
}

cluster_profile <- function(input, output, session,
                            distance_matrix,
                            cluster_method,
                            plot_width, plot_height) {
    # Reactive function holding data for clustering ============================
    cluster_data <- reactive({
        if (is.null(distance_matrix)) return()
        df <- cluster_data_dend(distance_matrix(), cluster_method())
        return(df)
    })

    # Dendrogram ===============================================================
    output$dendrogram <- renderPlot({
        if (is.null(cluster_data())) return()
        get_dendrogram(cluster_data())
    })

    output$cluster.ui <- renderUI({
        ns <- session$ns
        withSpinner(plotOutput(
            ns("dendrogram"),
            width = plot_width(),
            height = plot_height(),
            brush = brushOpts(
                id = ns("plot_brush"),
                delay = input$brush_delay,
                delayType = input$brush_policy,
                direction = input$brush_dir,
                resetOnNew = input$brush_reset
            )
        ))
    })

    # download clustered plot ==================================================
    output$download_cluster <- downloadHandler(
        filename = function() {
            "clustered_plot.pdf"
        },
        content = function(file) {
            ggsave(
                file,
                plot = get_dendrogram(cluster_data()),
                dpi = 300,
                device = "pdf",
                limitsize = FALSE
            )
        }
    )

    # Brushed cluster table ====================================================
    #' render brushed_cluster.table based on clicked point on dendrogram plot
    brushed_clusterGene <- reactive({
        dd.col <- cluster_data()
        dt <- dendro_data(dd.col)
        dt$labels$label <- levels(dt$labels$label)

        # get list of selected gene(s)
        if (is.null(input$plot_brush))
            return()
        else{
            top <- as.numeric(round(input$plot_brush$ymin))
            bottom <- as.numeric(round(input$plot_brush$ymax))
            df <- dt$labels[bottom:top,]
        }

        # return list of genes
        df <- df[complete.cases(df), 3]
        return(df)
    })

    output$brushed_cluster.table <- renderTable({
        if (is.null(input$plot_brush$ymin))
            return()

        data <- as.data.frame(brushed_clusterGene())
        data$number <- rownames(data)
        colnames(data) <- c("geneID", "No.")
        data <- data[, c("No.", "geneID")]
        data
    })

    #' download gene list from brushed_cluster.table
    output$download_cluster_genes <- downloadHandler(
        filename = function() {
            c("selectedClusteredGeneList.out")
        },
        content = function(file) {
            data_out <- brushed_clusterGene()
            write.table(
                data_out,
                file,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
        }
    )

    #' download distance matrix
    output$download_distance_matrix <- downloadHandler(
        filename = function() {
            c("distanceMatrixClustering.out")
        },
        content = function(file) {
            data_out <- distance_matrix()
            data_out <- as.matrix(data_out)
            write.table(
                data_out,
                file,
                col.names = TRUE,
                row.names = TRUE,
                quote = FALSE,
                sep = " \t"
            )
        }
    )

    #' Return the brushed genes
    return(brushed_clusterGene)
}
