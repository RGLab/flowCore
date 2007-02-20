\name{filterDetails-methods}
\docType{methods}
\alias{filterDetails-methods}
\alias{filterDetails}
\alias{filterDetails,filterResult,missing-method}
\alias{filterDetails,filterResult,ANY-method}
\alias{filterDetails<-,}
\alias{filterDetails<-, filterResult,character,setOperationFilter-method}
\alias{filterDetails<-, filterResult, character,filter-method}
\alias{filterDetails<-, filterResult, character,ANY-method}

\title{Obtain details about a filter operation}
\description{
A filtering operation captures details about its metadata and stores it in a \code{filterDetails} slot that is accessed using the
\code{filterDetails} method. Each set of metadata is indexed by the \code{filterId} of the filter allowing for all the metadata in
a complex filtering operation to be recovered after the final filtering. 
}
\section{Methods}{
\describe{

\item{result = "filterResult", filterId = "missing"}{When no particular \code{filterId} is specified all the details are returned}

\item{result = "filterResult", filterId = "ANY"}{You can also obtain a particular subset of details}
}}
\keyword{methods}

