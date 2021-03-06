% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trnobj_doc.R
\name{load_list}
\alias{load_list}
\title{Populate a TRN object from a list}
\usage{
TRNobj$load_list(TRNList)
}
\arguments{
\item{TRNList}{a TRN object converted to a list, e.g., with \code{as_list}}
}
\value{
None
}
\description{
This method populates all fields of a TRN object from the fields of an R list object. 
The list must have field names which exactly match the TRN field names.
}
