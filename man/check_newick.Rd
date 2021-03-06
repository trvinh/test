% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_phylotree.R
\name{check_newick}
\alias{check_newick}
\title{Check the validity of input newick tree}
\usage{
check_newick(tree, input_taxonID)
}
\arguments{
\item{tree}{input newick tree}

\item{input_taxonID}{list of all input taxon ID}
}
\value{
checking result (1 = missing parenthesis; 2 = missing comma;
3 = tree has singleton; or list of missing taxa in profile)
}
\description{
Check the validity of input newick tree
}
\examples{
data("pp_tree", package="phyloprofile")
check_newick(pp_tree, c("ncbi3702", "ncbi3711", "ncbi7029"))
}
\seealso{
\code{\link{get_input_taxa_id}} for getting input taxon IDs,
\code{\link{pp_tree}} for an example of input tree
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
