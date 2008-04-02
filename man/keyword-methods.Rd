\name{keyword-methods}
\docType{methods}
\alias{keyword}
\alias{keyword-methods}
\alias{keyword,flowFrame,missing-method}
\alias{keyword,flowFrame,function-method}
\alias{keyword,flowFrame,character-method}
\alias{keyword,flowFrame,list-method}
\alias{keyword<-}
\alias{keyword<-,flowFrame,list-method}
\alias{keyword<-,flowFrame,ANY-method}
\alias{keyword<-,flowFrame,character-method}

\title{Methods to retrieve keywords of a flowFrame}

\description{

  Accessor and replacement methods for items in the description slot
  (usually read in from a FCS file header). It lists the \code{keywords}
  and its values for a flowFrame specified by a character
  vector. Additional methods for \code{function} and \code{lists} exists
  for more programmatic access to the keywords.

}

\usage{
keyword(object, keyword)
}
\arguments{
  
  \item{object}{Object of class \code{\link{flowFrame}}.}
  
  \item{keyword}{character vector or list of potential keywords or
    function. If missing all keywords are returned.}
  
}

\section{Methods}{
  \describe{

    \item{object = "flowFrame", keyword = "character"}{Return values for
      all keywords from the \code{description} slot in \code{object}
      that match the character vector \code{keyword}.}

    \item{object = "flowFrame", keyword = "function"}{Apply the function in
      \code{keyword} on the \code{description} slot of \code{object}.}

    \item{object = "flowFrame", keyword = "list"}{Combine characters and
      functions in a list to select keyword values.}

    \item{object = "flowFrame", keyword = "missing"}{This is essentially an
      alias for \code{\link{description}} and returns all keyword-value
      pairs.}
    
  }
}

\author{N. LeMeur, F. Hahne, B. Ellis}

\seealso{\code{\link{description}}}

\examples{

 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 keyword(samp)

 keyword(samp,"FCSversion")

}

\keyword{methods}

