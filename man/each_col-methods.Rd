\name{each_col}
\docType{methods}
\alias{each_row}
\alias{each_row-methods}
\alias{each_row,flowFrame-method}
\alias{each_col}
\alias{each_col-methods}
\alias{each_col,flowFrame-method}

\title{Method to apply functions over flowFrame margins}

\description{
  
  Returns a vector or array of values obtained by applying a function to
  the margins of a flowFrame. This is equivalent of running
  \code{\link{apply}} on the output of \code{exprs(flowFrame)}.
  
}
\usage{
each_col(x, FUN, ...)
each_row(x, FUN, ...)
}

\arguments{
  \item{x}{Object of class \code{\linkS4class{flowFrame}}.}
  \item{FUN}{the function to be applied. In the case of
    functions like '+', '\%*\%', etc., the function name must be
    backquoted or quoted.}
  \item{...}{optional arguments to 'FUN'.}
}

\author{B. Ellis, N. LeMeur, F. Hahne}

\seealso{
  \code{\link{apply}}
}

\examples{

samp <- read.FCS(system.file("extdata", "0877408774.B08", package="flowCore"),
transformation="linearize")
each_col(samp, summary)

}

\keyword{methods}
