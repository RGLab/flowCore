\name{biexponentialTransform}
\alias{biexponentialTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
biexponentialTransform(transformationId, a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transformationId}{ ~~Describe \code{transformationId} here~~ }
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{c}{ ~~Describe \code{c} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{f}{ ~~Describe \code{f} here~~ }
  \item{w}{ ~~Describe \code{w} here~~ }
  \item{tol}{ ~~Describe \code{tol} here~~ }
  \item{maxit}{ ~~Describe \code{maxit} here~~ }
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
function(transformationId,a=.5,b=1,c=.5,d=1,f=0,w=0,
           tol=.Machine$double.eps^0.25,maxit=as.integer(5000)){
    new("transform",.Data=function(x){
        x <- .Call(biexponential_transform,x,a,b,c,d,f,w,tol,maxit)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
