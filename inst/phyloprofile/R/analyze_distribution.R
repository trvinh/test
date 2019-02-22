#' Distribution plots
#'
#' @export
#' @param data data for plotting (from reactive fn "presSpecAllDt")
#' @param var_id name of variable (either input$var1_id, input$var2_id or
#' "% present taxa"; from input$selected_dist)
#' @param var_type type of variable (either var1, var2 or presSpec)
#' @param percent percentage cutoff (from input$percent)
#' @param dist_text_size text size of distribution plot
#' @param dist_width width of distribution plot
#' (from input$dist_text_size)
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

analyze_distribution_ui <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            2,
            downloadButton(ns("plot_download_dist"), "Download plot")
        ),
        column(
            10,
            uiOutput(ns("dist_plot.ui"))
        )
    )
}

analyze_distribution <- function(input, output, session,
                                 data,
                                 var_id, var_type,
                                 percent,
                                 dist_text_size, dist_width){

    # render dist_plot.ui ------------------------------------------------------
    output$dist_plot.ui <- renderUI({
        ns <- session$ns
        withSpinner(plotOutput(ns("distribution_plot"),  width = dist_width()))
    })

    output$distribution_plot <- renderPlot(width = dist_width(), height = 356, {
        create_var_dist_plot(data(),
                             var_id(),
                             var_type(),
                             percent(),
                             dist_text_size())
    })

    output$plot_download_dist <- downloadHandler(
        filename = function() {
            paste0("distributionPlot.pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = create_var_dist_plot(
                    data(),
                    var_id(), var_type(),
                    percent(), dist_text_size()
                ),
                width = dist_width() * 0.056458333,
                height = 356 * 0.056458333,
                unit = "cm",
                dpi = 300, device = "pdf", limitsize = FALSE)
        }
    )
}
