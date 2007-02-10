\name{logicleTransform}
\alias{logicleTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logicleTransform(w = 0, r = 262144, d = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{ ~~Describe \code{w} here~~ }
  \item{r}{ ~~Describe \code{r} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function(w=0,r=262144,d=5,...) {
  if(w>d) stop("Negative range decades must be smaller than total number of decades")
  w = w*log(10)
  d = d*log(10)
  p = if(w==0) 1 else uniroot(function(p) -w+2*p*log(p)/(p+1),c(.Machine$double.eps,2*(w+d)))$root
  ##new("biexponentialTransformation",a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
 biexponentialTransform(a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
