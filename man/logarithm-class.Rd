\name{logarithm-class}
\docType{class}
\alias{logarithm-class}
\alias{logarithm}
\alias{eval,logarithm,missing,missing-method}

\title{Class "logarithm"  }
\description{	Logartithmic transformation of an argument is a transforation defined by the function 

 \deqn{f(parameter,a,b)=  ln(a*prarameter)*b \;\;\;\;\;\;\;\;\;a*parameter>0}
                            \deqn{0 \;\;\;\;\;\;\;a*parameter<=0}
}
\section{Objects from the Class}{
Objects can be created by calls to the constructor
\code{logarithm(parameters,a,b,transformationId)}
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} -non zero multiplicative constant }
    \item{\code{b}:}{Object of class \code{"numeric"} -non zero multiplicative constant  }
    \item{\code{parameters}:}{Object of class \code{"transformation"}-flow parameters to be transformed}
    \item{\code{transformationId}:}{Object of class \code{"character"}-unique ID to reference the transformation   }
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}

\section{Methods}{
No methods defined with class "logarithm" in the signature.
}

\references{Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5   }

\author{ 
Gopalakrishnan N, F.Hahne}

\note{ The logarithm transformation object can be evaluated using the eval method
 by passing the data frame as an argument.The transformed
parameters are returned as a matrix with a single column. (See example below) }

\seealso{
exponential,quadratic
}
\examples{
 dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  lg1<-logarithm(parameters="FSC-H",a=2,b=1,transformationId="lg1")
  transOut<-eval(lg1)(exprs(dat))
}
\keyword{classes}
