# CREATE MANUAL COLOR SCHEME FOR A LIST

#' Create qualitative colours
#' @param n number of colors
#' @param light light colors TRUE or FALSE
#' @return list of n different colors
#' @importFrom grDevices colorRampPalette
#' @source Modified based on https://gist.github.com/peterk87/6011397
#' @examples
#' \dontrun{
#' qualitative_colours(5)
#' }

qualitative_colours <- function(n, light=FALSE) {
    # Get a specified number of qualitative colours if possible.
    # This function will default to a continuous color scheme if there are more
    # than 21 colours needed.

    # rainbow12equal <- c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D",
    #                     "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF",
    #                     "#BF4DBF", "#BF4D86")
    rich12equal <- c(
        "#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466",
        "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300"
    )

    # Qualitative colour schemes by Paul Tol
    ifelse(
        n >= 19 & n <= 21,
        # return 21 qualitative color scheme
        return(
            grDevices::colorRampPalette(
                c(
                    "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                    "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
                    "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77",
                    "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
                    "#DD7788"
                )
            )(n)
        ),
        ifelse(
            n >= 16 & n <= 18,
            # up to 18 qualitative color scheme
            return(
                grDevices::colorRampPalette(
                    c(
                        "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                        "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711",
                        "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
                        "#771122", "#AA4455", "#DD7788"
                    )
                )(n)
            ),
            ifelse(
                n == 15,
                # 15 qualitative color scheme
                return(
                    grDevices::colorRampPalette(
                        c(
                            "#771155", "#AA4488", "#CC99BB", "#114477",
                            "#4477AA", "#77AADD", "#117777", "#44AAAA",
                            "#77CCCC", "#777711", "#AAAA44", "#DDDD77",
                            "#774411", "#AA7744", "#DDAA77", "#771122",
                            "#AA4455", "#DD7788"
                        )
                    )(n)
                ),
                ifelse(
                    n > 12 & n <= 14,
                    # 14 qualitative color scheme
                    return(
                        grDevices::colorRampPalette(
                            c(
                                "#882E72", "#B178A6", "#D6C1DE", "#1965B0",
                                "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                                "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D",
                                "#E8601C", "#DC050C"
                            )
                        )(n)
                    ),
                    ifelse(
                        n > 9 & n <= 12,
                        ifelse(
                            light,
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Set3'
                            )),
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Paired'
                            ))
                        ),
                        ifelse(
                            n <= 9,
                            ifelse(
                                light,
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Pastel1'
                                )),
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Set1'
                                ))
                            ),
                            # else(n > 21,
                            # If there are more than 21 qualitative colours,
                            # default to a continuous colour scheme,
                            # rich12equal in this case
                            return(grDevices::colorRampPalette(rich12equal)(n))
                        )
                    )
                )
            )
        )
    )
}

#' Get color for a list of items
#' @return list of colors for each element (same elements will have the same
#' color)
#' @param x input list
#' @export
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{qualitative_colours}}
#' @examples
#' items <- c("a", "b", "c")
#' get_qual_col_for_vector(items)

get_qual_col_for_vector <- function(x) {
    types <- unique(x)
    types <- types[!is.na(types)]
    type_colors <- qualitative_colours(length(types))

    colors_types <- as.vector(x)

    count_types <- 1
    count_colors <- 1

    while (count_types <= length(types)) {
        if (count_colors > length(type_colors)) {
            count_colors <- 1
        }

        colors_types[colors_types == types[count_types]] <-
            type_colors[count_colors]

        count_colors <- count_colors + 1
        count_types <- count_types + 1
    }

    return(unlist(colors_types))
}
