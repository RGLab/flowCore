\name{squareroot-class}
\docType{class}
\alias{squareroot-class}
\alias{squareroot}
\alias{eval,squareroot,missing,missing-method}

\title{Class "squareroot" }
\description{

Square root transform class defines a transformation defined 
by the function \deqn{f(parameter,a)= \sqrt{ |{\frac{parameter}{a}|}}}	

} 
\section{Objects from the Class}{
Objects can be created by calls to the constructor
\code{squareroot(parameters,a,transformationId)}
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} -non zero multiplicative 
constant }
    \item{\code{parameters}:}{Object of class \code{"transformation"}
-flow parameter to be transformed  }
    \item{\code{transformationId}:}{Object of class \code{"character"} 
 -unique ID to reference the transformation }
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}
\section{Methods}{
No methods defined with class "squareroot" in the signature.
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{The squareroot transformation object can be evaluated using the eval method
 by passing the data frame as an argument.The transformed
parameters are returned as a column vector. (See example below)  }



\seealso{

}
\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  sqrt1<-squareroot(parameters="FSC-H",a=2,transformationId="sqrt1")
  transOut<-eval(sqrt1)(exprs(dat))
}
\keyword{classes}
