% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_oma_browser.R
\name{create_profile_from_oma}
\alias{create_profile_from_oma}
\title{Create phylogenetic profile from a raw OMA dataframe}
\usage{
create_profile_from_oma(final_oma_df)
}
\arguments{
\item{final_oma_df}{raw OMA data for a list of input IDs}
}
\value{
phylogenetic profiles in long format
}
\description{
Create phylogenetic profile from a raw OMA dataframe
}
\examples{
oma_data <- get_data_for_one_oma("HUMAN29397", "OG")
create_profile_from_oma(oma_data)
}
\seealso{
\code{\link{get_data_for_one_oma}}
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
