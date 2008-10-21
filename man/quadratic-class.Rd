\name{quadratic-class}
\docType{class}
\alias{quadratic-class}
\alias{quadratic}
\alias{eval,quadratic,missing,missing-method}

\title{Class "quadratic" }

\description{
Quadratic transform class defines a transformation 
defined by the function \deqn{f(parameter,a)=a*parameter^2}	
}

\section{Objects from the Class}{
Objects can be created by calls to the constructor
\code{quadratic(parameters,a,transformationId)}
}

\section{Slots}{
\describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"}-non zero mutiplicative 
    constant}
    \item{\code{parameters}:}{Object of class \code{"transformation"}-flow 
    parameter to be transformed}
    \item{\code{transformationId}:}{Object of class \code{"character"}-unique 
    ID to reference the transformation}
}
}

\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", 
            distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform",
            distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform",
            distance 4.
}

\section{Methods}{
No methods defined with class "quadratic" in the signature.
}

\references{Gating-ML Candidate Recommendation for Gating Description in Flow
Cytometry V 1.5 
}

\author{Gopalakrishnan N, F.Hahne}

\note{ The quadratic transformation object can be evaluated using the eval method
 by passing the data frame as an argument.The transformed
parameters are returned as a column vector. (See example below)  }

\seealso{
dg1polynomial,ratio,squareroot
}

\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  quad1<-quadratic(parameters="FSC-H",a=2,transformationId="quad1")
  transOut<-eval(quad1)(exprs(dat))
}
\keyword{classes}
