\name{sinht-class}
\docType{class}
\alias{sinht-class}
\alias{sinht}
\alias{eval,sinht,missing,missing-method}

\title{Class "sinht" }

\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sinht", parameters, ...)}.
	 ~~ describe objects here ~~ 
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"}- non zero constant }
    \item{\code{b}:}{Object of class \code{"numeric"}- non zero constant}
    \item{\code{parameters}:}{Object of class \code{"transformation"} -flow parameters 
          to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to 
           reference the transformation }
  }
}

\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}

\section{Methods}{
No methods defined with class "sinht" in the signature.
}

\references{ Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5  }

\author{Gopalakrishnan N, F.Hahne  }

\note{ The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed parameters are returned as a matrix with a single column.(See example below)}

\seealso{
asinht	
}
\examples{
 dat <- read.FCS(system.file("extdata","0877408774.B08",  package="flowCore"))
 sinh1<-sinht(parameters="FSC-H",a=1,b=2000,transformationId="sinH1")
 transOut<-eval(sinh1)(exprs(dat))
}
\keyword{classes}
