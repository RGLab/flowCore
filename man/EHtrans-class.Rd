\name{EHtrans-class}
\docType{class}
\alias{EHtrans-class}
\alias{EHtrans}
\alias{eval,EHtrans,missing-method}

\title{Class "EHtrans" }
\description{	
 EH transformation of a parameter is defined by the function
      \deqn{EH(parameter,a,b)= 10^{(\frac{parameter}{a})} + \frac{b*parameter}{a}-1   ~~~~parameter>=0}
                            
      \deqn{-10^{(\frac{-parameter}{a})} + \frac{b*parameter}{a}+1  ~~~~parameter<0}

}
\section{Objects from the Class}{
Objects can be created by calls to the constructor
    \code{EHtrans(parameters,a,b,transformationId)}
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} -numeric constant greater than zero }
    \item{\code{b}:}{Object of class \code{"numeric"} -numeric constant greater than zero }
    \item{\code{parameters}:}{Object of class \code{"transformation"} - flow parameter to be 
    transformed}
    \item{\code{transformationId}:}{Object of class \code{"character"}- unique ID to reference the transformation }
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}
\section{Methods}{
No methods defined with class "EHtrans" in the signature.
}
\references{Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5 }
\author{ Gopalakrishnan N, F.Hahne}
\note{The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed
  parameters are returned as a matrix with a single column. (See example below) }

 

\seealso{
hyperlog
}
\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
                  package="flowCore"))
  eh1<-EHtrans("FSC-H",a=1250,b=4,transformationId="eh1")
  transOut<-eval(eh1)(exprs(dat))
}
\keyword{classes}
