\name{each_col}
\alias{each_col}
\alias{each_row}

\title{Function to apply methods over flowFrame margins}
\description{
  Returns a vector or array  of values obtained by applying a
  function to margins of flowFrame.
}
\usage{
each_col(x, FUN, ...)
each_row(x, FUN, ...)
}

\arguments{
  \item{x}{Object of class \code{\link[flowCore]{flowFrame}}.}
  \item{FUN}{the function to be applied: see Details. In the case of
    functions like '+', '%*%', etc., the function name must be
    backquoted or quoted.}
  \item{...}{optional arguments to 'FUN'.}
}

\author{B. Ellis, N. LeMeur}

\seealso{
  \code{apply}
}

\examples{
samp <- read.FCS(system.file("extdata",
        "0877408774.B08", package="flowCore"),
transformation="linearize")
each_col(samp, summary)
}
\keyword{methods}
