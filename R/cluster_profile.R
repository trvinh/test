#' Get data for calculating the distance matrix
#' @export
#' @usage get_data_clustering(data, profile_type, var1_aggregate_by, 
#'     var2_aggregate_by)
#' @param data processed profile data
#' @param profile_type type of data used for calculating the distance matrix. 
#' Either "binary" (consider only the presence/absence status of orthlogs), or
#' "var1"/"var2" for taking values of the additional variables into account.
#' @param var1_aggregate_by aggregate method for VAR1 (min, max, mean or median)
#' @param var2_aggregate_by aggregate method for VAR2 (min, max, mean or median)
#' @return A wide dataframe contains values for calculating distance matrix.
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @note Documented by Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{from_input_to_profile}}
#' @examples 
#' data("full_processed_profile", package="phyloprofile")
#' data <- full_processed_profile
#' profile_type <- "binary"
#' var1_aggregate_by <- "max"
#' var2_aggregate_by <- "mean"
#' get_data_clustering(data, profile_type, var1_aggregate_by, var2_aggregate_by)

get_data_clustering <- function(
    data, profile_type, var1_aggregate_by, var2_aggregate_by
) {
    supertaxon <- NULL
    presSpec <- NULL
    
    # remove lines where there is no found ortholog
    sub_data_heat <- subset(data, data$presSpec > 0)
    
    # transform data into wide matrix
    if (profile_type == "binary") {
        sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
        sub_data_heat$presSpec[sub_data_heat$presSpec > 0] <- 1
        sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
        wide_data <- spread(sub_data_heat, supertaxon, presSpec)
    } else {
        var <- profile_type
        
        sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", var)]
        sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
        
        # aggreagte the values by the selected method
        if (var == "var1") aggregate_by <- var1_aggregate_by
        else aggregate_by <- var2_aggregate_by
        
        sub_data_heat <- aggregate(
            sub_data_heat[, var],
            list(sub_data_heat$geneID, sub_data_heat$supertaxon),
            FUN = aggregate_by
        )
        
        colnames(sub_data_heat) <- c("geneID", "supertaxon", var)
        wide_data <- spread(sub_data_heat, supertaxon, var)
    }
    
    # set name for wide matrix as gene IDs
    dat <- wide_data[, 2:ncol(wide_data)]
    rownames(dat) <- wide_data[, 1]
    dat[is.na(dat)] <- 0
    
    return(dat)
}

#' Calculate the distance matrix
#' @export
#' @param profiles profile data for distance calculating
#' @param method distance calculation method ("euclidean", "maximum", 
#' "manhattan", "canberra", "binary", "distance_correlation", 
#' "mutual_information" or "pearson" for binary data; "distance_correlation" or
#' "mutual_information" for non-binary data)
#' @return A distance matrix for input phylogenetic profiles.
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @note Documented by Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{get_data_clustering}}
#' @examples 
#' data("full_processed_profile_large", package="phyloprofile")
#' data <- full_processed_profile_large
#' profile_type <- "binary"
#' profiles <- get_data_clustering(
#'     data, profile_type, var1_aggregate_by, var2_aggregate_by)
#' method <- "mutual_information"
#' get_distance_matrix(profiles, method)

get_distance_matrix <- function(profiles, method) {
    
    profiles <-  profiles[, colSums(profiles != 0) > 0]
    profiles <-  profiles[rowSums(profiles != 0) > 0, ]
    
    dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
    if (method %in% dist_methods) {
        distance_matrix <- dist(profiles, method = method)
    } else if (method == "distance_correlation") {
        matrix <- data.frame()
        for (i in seq_len(nrow(profiles))) { # rows
            for (j in seq_len(nrow(profiles))) { # columns
                if (i == j) {
                    matrix[i,i] = 1 # if this cell NA as.dist not work probably
                    break
                }
                dist <- dcor(unlist(profiles[i,]), unlist(profiles[j,]))
                # Swich the value so that the profiles with a high correlation
                # are clustered together
                matrix[i,j] <- 1 - dist
            }
        }
        
        profile_names <- rownames(profiles)
        colnames(matrix) <- profile_names[seq_len(length(profile_names)) - 1]
        rownames(matrix) <- profile_names
        distance_matrix <- as.dist(matrix)
    } else if (method == "mutual_information") {
        distance_matrix <- mutualInfo(as.matrix(profiles))
        distance_matrix <- max(distance_matrix, na.rm = TRUE) - distance_matrix
    } else if (method == "pearson") {
        distance_matrix <-  cor.dist(as.matrix(profiles))
    }
    
    return(distance_matrix)
}

#' Create a dendrogram tree from the distance matrix
#' @export
#' @param distance_matrix calculated distance matrix
#' @param cluster_method clustering method ("single", "complete", 
#' "average" for UPGMA, "mcquitty" for WPGMA, "median" for WPGMC, 
#' or "centroid" for UPGMC)
#' @return A dendrogram tree object
#' @importFrom stats as.dendrogram
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{get_data_clustering}}, 
#' \code{\link{get_distance_matrix}}, \code{\link{hclust}}
#' @examples
#' data("full_processed_profile_large", package="phyloprofile")
#' data <- full_processed_profile_large
#' profile_type <- "binary"
#' profiles <- get_data_clustering(
#'     data, profile_type, var1_aggregate_by, var2_aggregate_by)
#' dist_method <- "mutual_information"
#' distance_matrix <- get_distance_matrix(profiles, dist_method)
#' cluster_method <- "complete"
#' cluster_data_dend(distance_matrix, cluster_method)

cluster_data_dend <- function(distance_matrix, cluster_method) {
    if (is.null(distance_matrix)) return()
    dd.col <- as.dendrogram(hclust(distance_matrix, method = cluster_method))
    return(dd.col)
}

#' Plot dendrogram tree
#' @export
#' @param dd dendrogram object
#' @return A dendrogram plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{cluster_data_dend}}
#' @examples 
#' data("full_processed_profile_large", package="phyloprofile")
#' data <- full_processed_profile_large
#' profile_type <- "binary"
#' profiles <- get_data_clustering(
#'     data, profile_type, var1_aggregate_by, var2_aggregate_by)
#' dist_method <- "mutual_information"
#' distance_matrix <- get_distance_matrix(profiles, dist_method)
#' cluster_method <- "complete"
#' dd <- cluster_data_dend(distance_matrix, cluster_method)
#' get_dendrogram(dd)

get_dendrogram <- function(dd) {
    if (is.null(dd)) return()
    py <- dendextend::as.ggdend(dd)
    p <- ggplot(py, horiz = TRUE, theme = theme_minimal()) +
        theme(axis.title = element_blank(), axis.text.y = element_blank())
    return(p)
}