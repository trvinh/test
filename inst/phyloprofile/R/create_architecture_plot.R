#' Protein domain architecture plot
#'
#' @export
#' @param point_info() info of clicked point
#' (from reactive fn "point_infoDetail")
#' @param domain_info() domain information
#' (from reactive fn "get_domain_information")
#' @param label_archi_size lable size (from input$label_archi_size)
#' @param title_archi_size title size (from input$title_archi_size)
#' @param archi_height plot height (from input$archi_height)
#' @param archi_width plot width (from input$archi_width)
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

create_architecture_plot_ui <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("archi_plot.ui")),
        downloadButton(ns("archi_download"), "Download plot", class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        ),
        textOutput(ns("selected_domain"))
    )
}

create_architecture_plot <- function(input, output, session,
                                     point_info,
                                     domain_info,
                                     label_archi_size, title_archi_size,
                                     archi_height, archi_width){
    output$archi_plot <- renderPlot({
        if (is.null(nrow(domain_info()))) return()
        g <- create_archi_plot(point_info(),
                               domain_info(),
                               label_archi_size(),
                               title_archi_size())
        if (any(g == "ERR_0")) {
            msg_plot()
        } else {
            grid.draw(g)
        }
    })

    output$archi_plot.ui <- renderUI({
        ns <- session$ns
        if (is.null(nrow(domain_info()))) {
            msg <- paste0(
                "<p><em>Wrong domain file has been uploaded!
        Please check the correct format in
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/",
                "Input-Data#ortholog-annotations-eg-domains\"
        target=\"_blank\" rel=\"noopener\">our PhyloProfile wiki</a>.</em></p>"
            )
            HTML(msg)
        } else {
            withSpinner(
                plotOutput(
                    ns("archi_plot"),
                    height = archi_height(),
                    width = archi_width(),
                    click = ns("archi_click")
                )
            )
        }
    })

    output$archi_download <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            g <- create_archi_plot(point_info(),
                                   domain_info(),
                                   label_archi_size(),
                                   title_archi_size())
            grid.draw(g)
            ggsave(
                file, plot = g,
                width = archi_width() * 0.056458333,
                height = archi_height() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )

    output$selected_domain <- renderText({
        if (is.null(input$archi_click$y)) return()
        convertY(unit(input$archi_click$y, "npc"), "native")
    })
}

#' plot error message
#' @export
#' @param
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

msg_plot <- function() {
    msg <- paste(
        "No information about domain architecture!",
        "Please check:","if you uploaded the correct domain file/folder; or ",
        "if the selected genes (seed & ortholog) do exist in the uploaded file",
        "(please search for the corresponding seedID and hitID)",
        sep = "\n"
    )
    x <- c(1,2,3,4,5)
    y <- c(1,2,3,4,5)
    g <- ggplot(data.frame(x, y), aes(x,y)) +
        geom_point(color = "white") +
        annotate(
            "text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red"
        ) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank()) +
        ylim(0,1)
    return(g)
}
