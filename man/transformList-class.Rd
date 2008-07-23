\name{transformList-class}
\docType{class}
\alias{transformList-class}
\alias{transformList}
\alias{colnames,transformList-method}
\alias{c,transformList-method}

\title{Class "transformList"}
\description{Class "transformList"}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("transformList",
    ...)}, by calling the \code{\link{transform}} method with key-value
    pair arguments of the form \code{key} equals character and
    \code{value} equals function, or by using the constructor
    \code{transformList}. See below for details

  }

\usage{

transformList(from, tfun, to=from)

}
  
\arguments{

  \item{from,to}{ Characters giving the names of the measurement parameter
    on which to transform on and into which the result is supposed to be
    stored. If both are equal, the existing parameters will be
    overwritten. }
    
  \item{tfun}{ A list if functions or a character vector of the names of
    the functions used to transform the data. R's recycling rules apply,
    so a single function can be given to be used on all parameters. }

}


\section{Slots}{
  \describe{
    \item{\code{transforms}:}{Object of class \code{"list"}, where each
      list item is of class \code{\link{transformMap}}. }
  }
}

\section{Methods}{

  \item{colnames}{\code{signature(x = "transformList")}: This returns
    the names of the parameters that are to be transformed. }

  \item{c}{\code{signature(x = "transformList")}: Concatenate
    \code{transformList}s or regular lists and \code{transformLists}. }

  \item{\%on\%}{\code{signature(e1 = "transformList", e2 =
    "flowFrame")}: Perform a transformation using the
    \code{transformList} on a \code{\{link{flowFrame}} or
    \code{\{link{flowSet}}. }
  
}

\author{B. Ellis, F. Hahne}

\seealso{
  \code{\link{transform}}, \code{\link{transformMap}}
}
\examples{

tl <- transformList(c("FSC-H", "SSC-H"), list(log, asinh))
colnames(tl)
c(tl, transformList("FL1-H", "linearTransform"))

}
\keyword{classes}
