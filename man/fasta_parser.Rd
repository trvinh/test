% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_main_input.R
\name{fasta_parser}
\alias{fasta_parser}
\title{Parse multi-fasta input file}
\usage{
fasta_parser(input_file)
}
\arguments{
\item{input_file}{input multiple fasta file}
}
\value{
A data frame containing input data in long-format
}
\description{
Parse multi-fasta input file
}
\examples{
\dontrun{
input_file <- system.file(
    "extdata", "test.main.fasta", package = "phyloprofile", mustWork = TRUE
)
fasta_parser(input_file)
}
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
