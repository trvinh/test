#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @param info a list contains seed and orthologs IDs
#' @param domain_df dataframe contains domain info
#' @param label_archi_size lable size (in px)
#' @param title_archi_size title size (in px)
#' @importFrom dplyr filter
#' @importFrom gridExtra arrangeGrob
#' @return A plot as arrangeGrob object. Use grid::grid.draw(plot) to render.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{domain_plotting}}, \code{\link{sort_domains}},
#' \code{\link{parse_domain_input}}, \code{\link{get_qual_col_for_vector}}
#' @examples
#' seedID <- "OG_1009"
#' orthoID <- "A.thaliana@3702@241736"
#' info <- c(seedID, orthoID)
#' domain_file <- system.file(
#'     "extdata", "domain_files/OG_1009.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' domain_df <- parse_domain_input(seedID, domain_file, "file")
#' plot <- create_archi_plot(info, domain_df, 9, 9)
#' grid::grid.draw(plot)

create_archi_plot <- function(
    info,
    domain_df,
    label_archi_size, title_archi_size
){
    orthoID <- NULL

    # info
    group <- as.character(info[1])
    ortho <- as.character(info[2])

    # get sub dataframe based on selected group_id and orthoID
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")

    subdomain_df <- domain_df[grep(grepID, domain_df$seedID), ]
    subdomain_df$feature <- as.character(subdomain_df$feature)

    if (nrow(subdomain_df) < 1) {
        return(paste0("ERR_0"))
    } else {

        # ortho domains df
        ortho_df <- dplyr::filter(subdomain_df, orthoID == ortho)

        # seed domains df
        seed_df <- dplyr::filter(subdomain_df, orthoID != ortho)

        if (nrow(seed_df) == 0) seed_df <- ortho_df

        seed <- as.character(seed_df$orthoID[1])

        # return ERR_0 if seed_df and ortho_df are empty
        if (nrow(seed_df) == 0) return(paste0("ERR_0"))

        # change order of one dataframe's features
        # based on order of other df's features
        if (length(ortho_df$feature) < length(seed_df$feature)) {
            ordered_ortho_df <- ortho_df[order(ortho_df$feature), ]
            ordered_seed_df <- sort_domains(ordered_ortho_df, seed_df)
        } else {
            ordered_seed_df <- seed_df[order(seed_df$feature), ]
            ordered_ortho_df <- sort_domains(ordered_seed_df, ortho_df)
        }

        # join weight values and feature names
        if ("weight" %in% colnames(ordered_ortho_df)) {
            ordered_ortho_df$yLabel <- paste0(
                ordered_ortho_df$feature,
                " (",
                round(ordered_ortho_df$weight, 2),
                ")"
            )
        } else {
            ordered_ortho_df$yLabel <- ordered_ortho_df$feature
        }
        if ("weight" %in% colnames(ordered_seed_df)) {
            ordered_seed_df$yLabel <- paste0(
                ordered_seed_df$feature,
                " (",
                round(ordered_seed_df$weight, 2),
                ")"
            )
        } else {
            ordered_seed_df$yLabel <- ordered_seed_df$feature
        }

        # create color scheme for all features
        # the same features in seed & ortholog will have the same colors
        feature_seed <- levels(as.factor(ordered_seed_df$feature))
        feature_ortho <- levels(as.factor(ordered_ortho_df$feature))
        all_features <- c(feature_seed, feature_ortho)
        all_colors <- get_qual_col_for_vector(all_features)

        color_scheme <- structure(
            all_colors,
            .Names = all_features
        )

        # plotting
        sep <- "|"

        if ("length" %in% colnames(subdomain_df)) {
            plot_ortho <- domain_plotting(
                ordered_ortho_df,
                ortho,
                sep,
                label_archi_size,
                title_archi_size,
                min(subdomain_df$start),
                max(c(subdomain_df$end, subdomain_df$length)),
                color_scheme
            )
            plot_seed <- domain_plotting(
                ordered_seed_df,
                seed,
                sep,
                label_archi_size,
                title_archi_size,
                min(subdomain_df$start),
                max(c(subdomain_df$end, subdomain_df$length)),
                color_scheme
            )

        } else{
            plot_ortho <- domain_plotting(
                ordered_ortho_df,
                ortho,
                sep,
                label_archi_size,
                title_archi_size,
                min(subdomain_df$start),
                max(subdomain_df$end),
                color_scheme
            )
            plot_seed <- domain_plotting(
                ordered_seed_df,
                seed,
                sep,
                label_archi_size,
                title_archi_size,
                min(subdomain_df$start),
                max(subdomain_df$end),
                color_scheme
            )
        }

        # grid.arrange(plot_seed,plot_ortho,ncol=1)
        if (ortho == seed) {
            g <- gridExtra::arrangeGrob(plot_seed, ncol = 1)
        } else {
            seed_height <- length(levels(as.factor(ordered_seed_df$feature)))
            ortho_height <- length(levels(as.factor(ordered_ortho_df$feature)))

            g <- gridExtra::arrangeGrob(plot_seed, plot_ortho, ncol = 1,
                                        heights = c(seed_height, ortho_height))
        }
        return(g)
    }
}

#' Create architecure plot for a single protein
#' @usage domain_plotting(df, geneID, sep, label_size, title_size, min_start, 
#'     max_end, color_scheme)
#' @param df domain dataframe for ploting
#' @param geneID ID of seed or orthologous protein
#' @param sep separate indicator for title
#' @param label_size lable size
#' @param title_size title size
#' @param min_start the smallest start position of all domains
#' @param max_end the highest stop position of all domains
#' @param color_scheme color scheme for all domain types
#' @return A ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{get_qual_col_for_vector}},
#' \code{\link{parse_domain_input}}
#' @examples
#' \dontrun{
#' # get domain data
#' domain_file <- system.file(
#'     "extdata", "domain_files/OG_1009.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' domain_df <- parse_domain_input(seedID, domain_file, "file")
#' df <- domain_df[domain_df$orthoID == "A.thaliana@3702@241736",]
#' # create color scheme for all domain types
#' all_features <- levels(as.factor(df$feature))
#' all_colors <- get_qual_col_for_vector(all_features)
#' color_scheme <- structure(
#'     all_colors,
#'     .Names = all_features
#' )
#' # other parameters
#' geneID <- "A.thaliana@3702@241736"
#' sep <- "|"
#' label_size <- 9
#' title_size <- 9
#' min_start <- min(df$start)
#' max_end <- max(df$end)
#' # do plotting
#' domain_plotting(
#'     df,
#'     geneID,
#'     sep,
#'     label_size, title_size,
#'     min_start, max_end,
#'     color_scheme
#' )
#' }


domain_plotting <- function(df,
                            geneID,
                            sep,
                            label_size, title_size,
                            min_start, max_end,
                            color_scheme){
    feature <- NULL
    end <- NULL
    start <- NULL

    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        geom_segment(
            data = df,
            aes(y = feature, yend = feature, x = min_start, xend = max_end),
            color = "white",
            size = 0
        ) +
        scale_color_manual(values = color_scheme)

    # draw lines for representing sequence length
    if ("length" %in% colnames(df)) {
        gg <- gg + geom_segment(
            data = df,
            aes(x = 0, xend = length, y = feature, yend = feature),
            size = 1,
            color = "#b2b2b2"
        )
    }

    # draw line and points
    gg <- gg + geom_segment(
        data = df,
        aes(x = start, xend = end, y = feature, yend = feature),
        size = 1.5
    )
    gg <- gg + geom_point(
        data = df,
        aes(y = feature, x = start),
        color = "#b2b2b2",
        size = 3,
        shape = 3
    )
    gg <- gg + geom_point(
        data = df,
        aes(y = feature, x = end),
        color = "#edae52",
        size = 3,
        shape = 5
    )

    # draw dashed line for domain path
    gg <- gg + geom_segment(data = df[df$path == "Y", ],
                            aes(x = start, xend = end,
                                y = feature, yend = feature),
                            size = 3,
                            linetype = "dashed")

    # # add text above
    # gg <- gg + geom_text(data = df,
    #                      aes(x = (start + end) / 2,
    #                          y = feature, label = round(weight,2)),
    #                        color = "#9fb059",
    #                        size = descSize,
    #                        vjust = -0.75,
    #                        fontface = "bold",
    #                        family = "serif")

    # theme format
    title_mod <- gsub(":", sep, geneID)
    gg <- gg + scale_y_discrete(expand = c(0.075, 0),
                                breaks = df$feature,
                                labels = df$yLabel)
    gg <- gg + labs(title = paste0(title_mod), y = "Feature")
    gg <- gg + theme_minimal()
    gg <- gg + theme(panel.border = element_blank())
    gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(
        plot.title = element_text(face = "bold", size = title_size)
    )
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_text(size = label_size),
        axis.title.y = element_text(size = label_size),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
    # return plot
    return(gg)
}

#' Sort one domain dataframe based on the other domain dataframe
#' @description Sort domain dataframe of one protein (either seed or ortholog)
#' based on the dataframe of the its paired protein, in order to bring the
#' common domain feature in the same order which make it easy for comparing.
#' @param seed_df data of seed protein
#' @param ortho_df data of ortholog protein
#' @return sorted domain list (dataframe)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' # get domain data
#' domain_file <- system.file(
#'     "extdata", "domain_files/OG_1009.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' domain_df <- parse_domain_input(seedID, domain_file, "file")
#' # get seed_df and ortho_df
#' sub_df <- domain_df[domain_df$seedID == "OG_1009#A.thaliana@3702@241736",]
#' ortho_df <- dplyr::filter(sub_df, orthoID == "A.thaliana@3702@241736")
#' seed_df <- dplyr::filter(sub_df, orthoID != "A.thaliana@3702@241736")
#' # sort
#' sort_domains(seed_df, ortho_df)
#' }

sort_domains <- function(seed_df, ortho_df){
    orderNo <- NULL
    # get list of features in seed_df
    feature_list <- as.data.frame(levels(as.factor(seed_df$feature)))
    colnames(feature_list) <- c("feature")
    # and add order number to each feature
    feature_list$orderNo <- seq(length(feature_list$feature))

    # merge those info to ortho_df
    ordered_ortho_df <- merge(ortho_df, feature_list, all.x = TRUE)

    # sort ortho_df
    index <- with(ordered_ortho_df, order(orderNo))
    ordered_ortho_df <- ordered_ortho_df[index, ]

    #turn feature column into a character vector
    ordered_ortho_df$feature <- as.character(ordered_ortho_df$feature)
    #then turn it back into an ordered factor (to keep this order when plotting)
    ordered_ortho_df$feature <- factor(
        ordered_ortho_df$feature,
        levels = unique(ordered_ortho_df$feature)
    )
    #return sorted df
    return(ordered_ortho_df)
}
