#' Gene age estimation plot
#'
#' @export
#' @param data for gene age plot (from reactive fn gene_ageDf)
#' @param gene_age_width plot width (from input$gene_age_width)
#' @param gene_age_height plot width (from input$gene_age_height)
#' @param gene_age_text text size (from input$gene_age_text)
#' @return list of genes of a selected age
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

plot_gene_age_ui <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            2,
            downloadButton(ns("gene_age_plot_download"),
                           "Download plot", class = "butDL"),
            tags$head(
                tags$style(HTML(
                    ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
            )
        ),
        column(
            10,
            uiOutput(ns("gene_age.ui")),
            br(),
            em(
                h6("01_Species; 02_Family; 03_Class; 04_Phylum;
            05_Kingdom; 06_Superkingdom; 07_Last universal common ancestor;
           Undef_Genes have been filtered out")
            )
        ),
        hr(),
        column(
            4,
            downloadButton(ns("gene_age_table_download"), "Download gene list")
        ),
        column(
            8,
            tableOutput(ns("gene_age.table"))
        )
    )
}

plot_gene_age <- function(input, output, session,
                          data,
                          gene_age_width, gene_age_height, gene_age_text) {
    # render gene age plot -----------------------------------------------------
    output$gene_agePlot <- renderPlot({
        if (is.null(data())) return()
        create_gene_age_plot(gene_age_plotDf(data()), gene_age_text())
    })

    output$gene_age.ui <- renderUI({
        ns <- session$ns
        withSpinner(
            plotOutput(
                ns("gene_agePlot"),
                width = 600 * gene_age_width(),
                height = 150 * gene_age_height(),
                click = ns("plot_click_gene_age")
            )
        )
    })

    output$gene_age_plot_download <- downloadHandler(
        filename = function() {
            "gene_age_plot.pdf"
        },
        content = function(file) {
            ggsave(
                file,
                plot = create_gene_age_plot(gene_age_plotDf(data()),
                                            gene_age_text()),
                width = 600 * gene_age_width() * 0.056458333,
                height = 150 * gene_age_height() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf"
            )
        }
    )

    # render genAge.table based on clicked point on gene_agePlot ---------------
    selectedgene_age <- reactive({
        if (is.null(data())) return()
        selected_gene <- get_selected_gene_age(data(),
                                               input$plot_click_gene_age$x)
        return(selected_gene)
    })

    output$gene_age.table <- renderTable({
        if (is.null(data())) return()
        if (is.null(input$plot_click_gene_age$x)) return()

        data <- as.data.frame(selectedgene_age())
        data$number <- rownames(data)
        colnames(data) <- c("geneID", "No.")
        data <- data[, c("No.", "geneID")]
        data
    })

    output$gene_age_table_download <- downloadHandler(
        filename = function(){
            c("selectedGeneList.out")
        },
        content = function(file){
            data_out <- selectedgene_age()
            write.table(data_out, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(selectedgene_age)
}

#' get list of gene for a selected age
#' @param gene_ageDf data of estimated gene age (from fn "estimate_gene_age")
#' @param clicked_x x coordinate of selected age

get_selected_gene_age <- function(gene_ageDf, clicked_x){
    data <- gene_ageDf

    # calculate the coordinate range for each age group
    rangeDf <- plyr::count(data, c("age"))

    rangeDf$percentage <- round(rangeDf$freq / sum(rangeDf$freq) * 100)
    rangeDf$rangeStart[1] <- 0
    rangeDf$rangeEnd[1] <- rangeDf$percentage[1]
    if (nrow(rangeDf) > 1) {
        for (i in 2:nrow(rangeDf)) {
            rangeDf$rangeStart[i] <- rangeDf$rangeEnd[i - 1] + 1
            rangeDf$rangeEnd[i] <-
                rangeDf$percentage[i] + rangeDf$rangeEnd[i - 1]
        }
    }

    # get list of gene for selected age
    if (is.null(clicked_x)) return()
    else{
        corX <- 100 - round(clicked_x)
        selectAge <- {
            as.character(rangeDf[rangeDf$rangeStart <= corX
                                 & rangeDf$rangeEnd >= corX, ]$age)
        }
        subData <- subset(data, age == selectAge)
        data <- data[data$age == selectAge, ]
    }

    # return list of genes
    geneList <- levels(as.factor(subData$geneID))
    return(geneList)
}
