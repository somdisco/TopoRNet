% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trnobj_doc.R
\name{prune_CONN_edge_weight}
\alias{prune_CONN_edge_weight}
\title{Prune CONN by a minimum edge weight}
\usage{
TRNobj$prune_CONN_edge_weight(min_weight)
}
\arguments{
\item{min_weight}{minimum weight allowed in the pruned graph. 
Any edges in \code{CONN_EL} whose \code{CONN} value is < \code{min_weight} will be pruned.}
}
\value{
None
}
\description{
Prune the CONN graph by removing edges whose weights are \strong{strictly less than} a given minimum edge weight.
}
