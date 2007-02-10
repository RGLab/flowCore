\name{quadraticTransform}
\alias{quadraticTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
quadraticTransform(transformationId, a, b, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transformationId}{ ~~Describe \code{transformationId} here~~ }
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{c}{ ~~Describe \code{c} here~~ }
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
function(transformationId,a,b,c){
  if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
  if(!is.double(c))
       stop("c must be numeric")
    new("transform",.Data=function(x){
        x <- a*x^2 + b*x + c
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
