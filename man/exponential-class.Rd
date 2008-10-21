\name{exponential-class}
\docType{class}
\alias{exponential-class}
\alias{exponential}
\alias{eval,exponential,missing,missing-method}

\title{Class "exponential" }
\description{ Exponential transform class defines a trasformation given by the
function
          \deqn{f(parameter,a,b)=e^{parameter/b}*\frac{1}{a}}

}
\section{Objects from the Class}{
Objects can be created by calls to the constructor\code{exponential(parameters,a,b)}.

}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"}- non zero constant }
    \item{\code{b}:}{Object of class \code{"numeric"}- non zero constant }
    \item{\code{parameters}:}{Object of class \code{"transformation"}- flow 
          parameter to be transformed }
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
No methods defined with class "exponential" in the signature.
}

\references{Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5  }

\author{ Gopalakrishnan N, F.Hahne }

\note{ The exponential transformation object can be evaluated using the eval method
 by passing the data frame as an argument.The transformed
parameters are returned as a matrix with a single column}


\seealso{
logarithm
}

\examples{
 dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  exp1<-exponential(parameters="FSC-H",a=1,b=37,transformationId="exp1")
  transOut<-eval(exp1)(exprs(dat))
}
\keyword{classes}
