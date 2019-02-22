#' Create data for main profile plot
#' @export
#' @param data_heat a data frame contains processed profiles
#' @return A data frame contains data for main profile plot.
#' @importFrom stats na.omit
#' @rawNamespace import(data.table, except = c(set, melt))
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{from_input_to_profile}}
#' @examples 
#' data("full_processed_profile", package="phyloprofile")
#' data_main_plot(full_processed_profile)

data_main_plot <- function(data_heat){
    paralogNew <- NULL
    
    # reduce number of inparalogs based on filtered dataHeat
    data_heat_tb <- data.table(na.omit(data_heat))
    data_heat_tb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    data_heat_tb <- data.frame(data_heat_tb[, c("geneID",
                                                "supertaxon",
                                                "paralogNew")])
    
    data_heat <- merge(
        data_heat, data_heat_tb,
        by = c("geneID", "supertaxon"),
        all.x = TRUE
    )
    data_heat$paralog <- data_heat$paralogNew
    data_heat <- data_heat[!duplicated(data_heat), ]
    
    # remove unneeded dots
    data_heat$presSpec[data_heat$presSpec == 0] <- NA
    data_heat$paralog[data_heat$presSpec < 1] <- NA
    data_heat$paralog[data_heat$paralog == 1] <- NA
    
    return(data_heat)
}

#' Create data for customized profile plot
#' @description Create data for customized profile plot based on a selected 
#' list of genes and/or taxa.
#' @export
#' @param data_heat a data frame contains processed profiles
#' @param selected_taxa subset of taxa
#' @param selected_seq subset of sequences
#' @return A data frame contains data for customized profile plot.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{from_input_to_profile}}
#' @examples 
#' data("full_processed_profile", package="phyloprofile")
#' selected_taxa <- c("Mammalia", "Echinoidea", "Gunneridae")
#' selected_seq <- "all"
#' data_customized_plot(full_processed_profile, selected_taxa, selected_seq)

data_customized_plot <- function(data_heat, selected_taxa, selected_seq){
    geneID <- NULL
    supertaxonMod <- NULL
    
    # process data
    data_heat$supertaxonMod <- {
        substr(
            data_heat$supertaxon, 6, nchar(as.character(data_heat$supertaxon))
        )
    }
    
    if (selected_taxa[1] == "all" & selected_seq[1] != "all") {
        # select data from dataHeat for selected sequences only
        data_heat <- subset(data_heat, geneID %in% selected_seq)
    } else if (selected_seq[1] == "all" & selected_taxa[1] != "all") {
        # select data from dataHeat for selected taxa only
        data_heat <- subset(data_heat, supertaxonMod %in% selected_taxa)
    } else {
        # select data from dataHeat for selected sequences and taxa
        data_heat <- subset(data_heat,
                            geneID %in% selected_seq
                            & supertaxonMod %in% selected_taxa)
    }
    
    # remove unneeded dots
    data_heat$presSpec[data_heat$presSpec == 0] <- NA
    data_heat$paralog[data_heat$presSpec < 1] <- NA
    data_heat$paralog[data_heat$paralog == 1] <- NA
    
    return(data_heat)
}


#' Create profile heatmap plot
#' @export
#' @param data data for heatmap plot
#' @param plot_parameter plot parameters (type of x-axis "taxa" or "genes";
#' names of 2 variables; colors for lowest and highest value of variable 1; 
#' colors for lowest and highest value of variable 2; color of co-orthologs; 
#' text sizes for x, y axis and legend; legend position "top", "bottom", 
#' "right", "left" or "none"; zoom ratio of the co-ortholog dots from -1 to 3;
#' angle of x-axis from 0 to 90; 
#' show/hide separate line for reference taxon 1/0; 
#' enable/disable coloring gene categories TRUE/FALSE)
#' @return A profile heatmap plot as ggplot object.
#' @importFrom plyr mapvalues
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 geom_tile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{data_main_plot}}, \code{\link{data_customized_plot}}
#' @examples 
#' data("full_processed_profile", package="phyloprofile")
#' plot_df <- data_main_plot(full_processed_profile)
#' plot_parameter <- list(
#'     "x_axis" = "taxa",
#'     "var1_id" = "FAS",
#'     "var2_id"  = "Traceability",
#'     "low_color_var1" =  "#FF8C00",
#'     "high_color_var1" = "#4682B4",
#'     "low_color_var2" = "#FFFFFF",
#'     "high_color_var2" = "#F0E68C",
#'     "para_color" = "#07D000",
#'     "x_size" = 8,
#'     "y_size" = 8,
#'     "legend_size" = 8,
#'     "main_legend" = "top",
#'     "dot_zoom" = 0,
#'     "x_angle" = 60,
#'     "guideline" = 0,
#'     "color_by_group" = FALSE
#' )
#' 
#' heatmap_plotting(plot_df, plot_parameter)

heatmap_plotting <- function(data, plot_parameter){
    geneID <- NULL
    supertaxon <- NULL
    group <- NULL
    var1 <- NULL
    var2 <- NULL
    presSpec <- NULL
    paralog <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    
    # parameters
    x_axis <- plot_parameter$x_axis
    var1_id <- plot_parameter$var1_id
    var2_id <- plot_parameter$var2_id
    low_color_var1 <- plot_parameter$low_color_var1
    high_color_var1 <- plot_parameter$high_color_var1
    low_color_var2 <- plot_parameter$low_color_var2
    high_color_var2 <- plot_parameter$high_color_var2
    para_color <- plot_parameter$para_color
    x_size <- plot_parameter$x_size
    y_size <- plot_parameter$y_size
    legend_size <- plot_parameter$legend_size
    main_legend <- plot_parameter$main_legend
    dot_zoom <- plot_parameter$dot_zoom
    x_angle <- plot_parameter$x_angle
    guideline <- plot_parameter$guideline
    color_by_group <- plot_parameter$color_by_group
    
    # rescale numbers of paralogs
    data$paralog <- as.numeric(data$paralog)
    if (length(unique(na.omit(data$paralog))) > 0) {
        max_paralog <- max(na.omit(data$paralog))
        data$paralogSize <- (data$paralog / max_paralog) * 3
    }
    
    # remove prefix number of taxa names but keep the order
    data$supertaxon <- {
        mapvalues(
            warn_missing = FALSE,
            data$supertaxon,
            from = as.character(data$supertaxon),
            to = substr(as.character(data$supertaxon),
                        6,
                        nchar(as.character(data$supertaxon)))
        )
    }
    
    # format plot
    if (x_axis == "genes") {
        p <- ggplot(data, aes(x = geneID, y = supertaxon))
    } else{
        p <- ggplot(data, aes(y = geneID, x = supertaxon))
    }
    
    if (color_by_group == TRUE) {
        p <- p + geom_tile(aes(fill = factor(group)), alpha = 0.3)
    } else {
        if (length(unique(na.omit(data$var2))) != 1) {
            p <- p + scale_fill_gradient(
                low = low_color_var2,
                high = high_color_var2,
                na.value = "gray95",
                limits = c(0, 1)
            ) +  #fill color (var2)
                geom_tile(aes(fill = var2))    # filled rect (var2 score)
        }
    }
    
    if (length(unique(na.omit(data$presSpec))) < 3) {
        if (length(unique(na.omit(data$var1))) == 1) {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1),
                                size = data$presSpec * 5 * (1 + dot_zoom),
                                na.rm = TRUE, show.legend = FALSE)
        } else {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1),
                                size = data$presSpec * 5 * (1 + dot_zoom),
                                na.rm = TRUE)
            # color of the corresponding aes (var1)
            p <- p + scale_color_gradient(
                low = low_color_var1,
                high = high_color_var1,
                limits = c(0, 1)
            )
        }
    } else {
        if (length(unique(na.omit(data$var1))) == 1) {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(size = presSpec),
                                color = "#336a98",
                                na.rm = TRUE)
        } else {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1, size = presSpec),
                                na.rm = TRUE)
            # color of the corresponding aes (var1)
            p <- p +
                scale_color_gradient(
                    low = low_color_var1, high = high_color_var1,
                    limits = c(0, 1)
                )
        }
    }
    
    # plot inparalogs (if available)
    if (length(unique(na.omit(data$paralog))) > 0) {
        p <- p + geom_point(data = data,
                            aes(size = paralog),
                            color = para_color,
                            na.rm = TRUE,
                            show.legend = TRUE)
        p <- p + guides(size = guide_legend(title = "# of co-orthologs"))
        
        # to tune the size of circles
        p <- p +
            scale_size_continuous(
                range = c(
                    min(na.omit(data$paralogSize)) * (1 + dot_zoom),
                    max(na.omit(data$paralogSize)) * (1 + dot_zoom)
                )
            )
    } else {
        # remain the scale of point while filtering
        present_vl <- data$presSpec[!is.na(data$presSpec)]
        
        # to tune the size of circles;
        # use "floor(value*10)/10" to round "down" the value with one decimal nr
        p <- p +
            scale_size_continuous(
                range = c(
                    (floor(min(present_vl) * 10) / 10 * 5) * (1 + dot_zoom),
                    (floor(max(present_vl) * 10) / 10 * 5) * (1 + dot_zoom)
                )
            )
    }
    
    if (color_by_group == FALSE) {
        p <- p + guides(fill = guide_colourbar(title = var2_id),
                        color = guide_colourbar(title = var1_id))
    } else {
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = var1_id))
    }
    
    base_size <- 9
    
    # guideline for separating ref species
    if (guideline == 1) {
        if (x_axis == "genes") {
            p <- p + labs(y = "Taxon")
            p <- p + geom_hline(yintercept = 0.5, colour = "dodgerblue4")
            p <- p + geom_hline(yintercept = 1.5, colour = "dodgerblue4")
        } else{
            p <- p + labs(x = "Taxon")
            p <- p + geom_vline(xintercept = 0.5, colour = "dodgerblue4")
            p <- p + geom_vline(xintercept = 1.5, colour = "dodgerblue4")
        }
    }
    
    # format theme
    p <- p + theme_minimal()
    p <- p + theme(
        axis.text.x = element_text(angle = x_angle, hjust = 1, size = x_size),
        axis.text.y = element_text(size = y_size),
        axis.title.x = element_text(size = x_size),
        axis.title.y = element_text(size = y_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size),
        legend.position = main_legend
    )
    
    # return plot
    return(p)
}

#' Highlight gene and/or taxon of interest on the profile plot
#' @export
#' @usage highlight_profile_plot(data, plot_parameter, taxon_highlight, 
#'     rank_name, gene_highlight)
#' @param data data for heatmap plot
#' @param plot_parameter plot parameters (type of x-axis "taxa" or "genes";
#' names of 2 variables; colors for lowest and highest value of variable 1; 
#' colors for lowest and highest value of variable 2; color of co-orthologs; 
#' text sizes for x, y axis and legend; legend position "top", "bottom", 
#' "right", "left" or "none"; zoom ratio of the co-ortholog dots from -1 to 3;
#' angle of x-axis from 0 to 90; 
#' show/hide separate line for reference taxon 1/0; 
#' enable/disable coloring gene categories TRUE/FALSE)
#' @param taxon_highlight taxon of interst
#' @param rank_name working taxonomy rank
#' @param gene_highlight gene of interest
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{data_main_plot}}, \code{\link{data_customized_plot}}
#' @examples 
#' data("full_processed_profile", package="phyloprofile")
#' plot_df <- data_main_plot(full_processed_profile)
#' plot_parameter <- list(
#'     "x_axis" = "taxa",
#'     "var1_id" = "FAS",
#'     "var2_id"  = "Traceability",
#'     "low_color_var1" =  "#FF8C00",
#'     "high_color_var1" = "#4682B4",
#'     "low_color_var2" = "#FFFFFF",
#'     "high_color_var2" = "#F0E68C",
#'     "para_color" = "#07D000",
#'     "x_size" = 8,
#'     "y_size" = 8,
#'     "legend_size" = 8,
#'     "main_legend" = "top",
#'     "dot_zoom" = 0,
#'     "x_angle" = 60,
#'     "guideline" = 0,
#'     "color_by_group" = FALSE
#' )
#' taxon_highlight <- "Mammalia"
#' rank_name <- "class"
#' gene_highlight <- "OG_1019"
#' highlight_profile_plot(
#'     plot_df, plot_parameter, taxon_highlight, rank_name, gene_highlight
#' )

highlight_profile_plot <- function(
    data,
    plot_parameter,
    taxon_highlight,
    rank_name,
    gene_highlight
){
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    
    # get heatmap
    p <- heatmap_plotting(data, plot_parameter)
    
    # highlight taxon
    if (taxon_highlight != "none") {
        # get selected highlight taxon ID
        nameReduced_file <- paste(
            system.file(package="phyloprofile"),
            "phyloprofile/data/taxonNamesReduced.txt",
            sep="/"
        )
        
        if (!file.exists(nameReduced_file)) {
            fileURL <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/taxonNamesReduced.txt"
            )
            res <- tryCatch(
                utils::download.file(
                    fileURL,
                    destfile = nameReduced_file,
                    method="auto"
                ),
                error=function(e) 1
            )
        }
        
        taxa_list <- as.data.frame(read.table(
            nameReduced_file,
            sep = "\t",
            header = TRUE
        ))
        
        taxon_highlight_id <- {
            taxa_list$ncbiID[
                taxa_list$fullName == taxon_highlight 
                & taxa_list$rank == rank_name
            ]
        }
        
        if (length(taxon_highlight_id) == 0L) {
            taxon_highlight_id <- {
                taxa_list$ncbiID[taxa_list$fullName == taxon_highlight]
            }
        }
        
        # get taxonID together with it sorted index
        highlight_taxon <- {
            toString(
                data[data$supertaxonID == taxon_highlight_id, 2][1]
            )
        }
        
        # get index
        selected_index <- as.numeric(
            as.character(substr(highlight_taxon, 2, 4))
        )
        
        # draw a rect to highlight this taxon's column
        if (plot_parameter$x_axis == "taxa") {
            rect <- data.frame(
                xmin = selected_index - 0.5,
                xmax = selected_index + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
        } else {
            rect <- data.frame(
                ymin = selected_index - 0.5,
                ymax = selected_index + 0.5,
                xmin = -Inf,
                xmax = Inf
            )
        }
        
        p <- p + geom_rect(
            data = rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "yellow",
            alpha = 0.3,
            inherit.aes = FALSE
        )
    }
    
    # highlight gene
    if (gene_highlight != "none") {
        # get selected highlight gene ID
        gene_highlight <- gene_highlight
        
        # get index
        all_genes <- levels(data$geneID)
        selected_index <- match(gene_highlight, all_genes)
        
        # draw a rect to highlight this taxon's column
        if (plot_parameter$x_axis == "taxa") {
            rect <- data.frame(
                ymin = selected_index - 0.5,
                ymax = selected_index + 0.5,
                xmin = -Inf,
                xmax = Inf
            )
        } else {
            rect <- data.frame(
                xmin = selected_index - 0.5,
                xmax = selected_index + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
        }
        
        p <- p + geom_rect(
            data = rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "yellow",
            alpha = 0.3,
            inherit.aes = FALSE
        )
    }
    
    return(p)
}
