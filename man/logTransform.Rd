\name{logTransform}
\alias{logTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logTransform(transformationId, logbase = 10, r, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transformationId}{ ~~Describe \code{transformationId} here~~ }
  \item{logbase}{ ~~Describe \code{logbase} here~~ }
  \item{r}{ ~~Describe \code{r} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(transformationId,logbase=10,r,d){
     if(!is.double(r) || r <= 0)
       stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
       stop("d must be numeric")
      if(!is.double(r) || r <=0)
       stop("r must be numeric and positive")
    if(!is.double(logbase) || logbase <= 1)
       stop("logabse must be a pnumeric greater than 1")
    new("transform",.Data=function(x){
        x <- log(x,logbase)*(r/d)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
