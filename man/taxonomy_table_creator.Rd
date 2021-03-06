% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_taxonomy_matrix.R
\name{taxonomy_table_creator}
\alias{taxonomy_table_creator}
\title{Align ncbi taxonomy IDs of all taxa into a sorted rank list}
\usage{
taxonomy_table_creator(idListFile, rankListFile)
}
\arguments{
\item{idListFile}{a text file whose each row is a rank+ID list of a taxon}

\item{rankListFile}{a text file whose each row is a rank list of a taxon}
}
\value{
An aligned taxonomy dataframe
}
\description{
Align ncbi taxonomy IDs of all taxa into a sorted rank list
}
\examples{
idListFile <- system.file(
    "extdata", "data/idList.txt", package = "phyloprofile", mustWork = TRUE
)
rankListFile <- system.file(
    "extdata", "data/rankList.txt", package = "phyloprofile", mustWork = TRUE
)
taxonomy_table_creator(idListFile, rankListFile)
}
\seealso{
\code{\link{rankIndexing}}
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
