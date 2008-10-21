\name{asinht-class}
\docType{class}
\alias{asinht-class}
\alias{asinht}
\alias{eval,asinht,missing,missing-method}

\title{Class "asinht"  }
\description{ Inverse hyperbolic sin transfromation is defined by the function
\deqn{f(parameter,a,b)=sinh^{-1}(a*parameter)*b}
}
\section{Objects from the Class}{
Objects can be created by calls  to the constructor
\code{asinht(parameter,a,b,transformationId)}

}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} -non zero constant }
    \item{\code{b}:}{Object of class \code{"numeric"} -non zero constant }
    \item{\code{parameters}:}{Object of class \code{"transformation"} -flow parameter
to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to reference the transformation}
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}
\section{Methods}{
No methods defined with class "asinht" in the signature.
}
\references{Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5  }

\author{ Gopalakrishnan N, F.Hahne }

\note{ The inverse hyperbolic sin transformation object can be evaluated using the eval method
 by passing the data frame as an argument.The transformed parameters are returned as a matrix with a single column. (See example below)  }

\seealso{
sinht
}

\examples{
 dat <- read.FCS(system.file("extdata","0877408774.B08",  package="flowCore"))
  asinh1<-asinht(parameters="FSC-H",a=2,b=1,transformationId="asinH1")
  transOut<-eval(asinh1)(exprs(dat))
}
\keyword{classes}
