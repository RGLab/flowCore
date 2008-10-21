\name{compensatedParameter-class}
\docType{class}
\alias{compensatedParameter-class}
\alias{compensatedParameter}
\alias{eval,compensatedParameter,missing,missing-method}

\title{Class "compensatedParameter" }

\description{

  Emission spectral overlap can be corrected by subtracting the the
  amount of spectral overlap from the total detected signals. This
  compensation process can be described by using spillover matrices.

  compensatedParameter objects allow for compensation of specific
  parameters the user is interested in by creating compensatedParameter
  objects and evaluating them.This allows for use of
  compensatedParameter in gate definitions.

}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form
  \code{compensatedParameter(parameters,spillRefId,transformationId,searchEnv)}.
	
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"character"} -flow parameters to
 be  compensated }
    \item{\code{spillRefId}:}{Object of class \code{"character"} -name of the
 compensation object (The compensation object contains the spillover Matrix) }
    \item{\code{searchEnv}:}{Object of class \code{"environment"} -environment in
 which the compensation object is defined }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique Id to
  reference the compensatedParameter object }
  }
}
\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}
\section{Methods}{
No methods defined with class "compensatedParameter" in the signature.
}
\references{ ~~ }
\author{Gopalakrishnan N,F.Hahne}
\note{  The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed
  parameters are returned as a matrix with a single column. (See example below)}

\seealso{
compensation

}
\examples{
samp   <- read.flowSet(path=system.file("extdata", "compdata", "data", package="flowCore"))
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
comp.mat

## create a compensation object 
comp <- compensation(comp.mat,compensationId="comp1")
## create a compensated parameter object 
cPar1<-compensatedParameter(c("FL1-H","FL3-H"),"comp",searchEnv=.GlobalEnv)
compOut<-eval(cPar1)(exprs(samp[[1]]))

}
\keyword{classes}
