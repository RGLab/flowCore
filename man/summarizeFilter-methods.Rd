\name{summarizeFilter-methods}

\docType{methods}

\alias{summarizeFilter-methods}
\alias{summarizeFilter}
\alias{summarizeFilter,filterResult,filter-method}
\alias{summarizeFilter,filterResult,filterReference-method}
\alias{summarizeFilter,filterResult,parameterFilter-method}
\alias{summarizeFilter,filterResult,subsetFilter-method}
\alias{summarizeFilter,logicalFilterResult,norm2Filter-method}
\alias{summarizeFilter,logicalFilterResult,parameterFilter-method}
\alias{summarizeFilter,multipleFilterResult,curv1Filter-method}
\alias{summarizeFilter,multipleFilterResult,curv2Filter-method}
\alias{summarizeFilter,multipleFilterResult,parameterFilter-method}


\title{ Methods for function summarizeFilter }

\description{
  
  Internal methods to populate the \code{filterDetails} slot of a
  \code{\link{filterResult}} object.
  
}

\section{Methods}{
  \describe{
  
    \item{result = "filterResult", filter = "filter"}{
      \code{summarizeFilter} methods are called during the process of
      filtering. Their output is a list, and it can be arbitraty data
      that should be stored along with the results of a filtering
      operation. }

\item{result = "filterResult", filter = "filterReference"}{ see above }

\item{result = "filterResult", filter = "parameterFilter"}{ see above }

\item{result = "filterResult", filter = "subsetFilter"}{ see above }

\item{result = "logicalFilterResult", filter = "norm2Filter"}{ see above }

\item{result = "logicalFilterResult", filter = "parameterFilter"}{ see above }

\item{result = "multipleFilterResult", filter = "curv1Filter"}{ see above }

\item{result = "multipleFilterResult", filter = "curv2Filter"}{ see above }

\item{result = "multipleFilterResult", filter = "parameterFilter"}{ see above }
}}

\keyword{methods}

