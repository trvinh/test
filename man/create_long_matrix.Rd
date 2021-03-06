% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_main_input.R
\name{create_long_matrix}
\alias{create_long_matrix}
\title{Create a long matrix format for all kinds of input file}
\usage{
create_long_matrix(input_file)
}
\arguments{
\item{input_file}{input file in orthoXML, multiple FASTA, wide or long
matrix format.}
}
\value{
A data frame that contains input data in long-format
}
\description{
Create a long matrix format for all kinds of input file
}
\examples{
input_file <- system.file(
    "extdata", "test.main.wide", package = "phyloprofile", mustWork = TRUE
)
create_long_matrix(input_file)
}
\seealso{
\code{\link{xml_parser}}, \code{\link{fasta_parser}},
\code{\link{wide_to_long}}
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
