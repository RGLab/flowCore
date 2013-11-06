\name{ratiotGml2-class}
\docType{class}
\alias{ratiotGml2-class}
\alias{ratiotGml2}
\alias{eval,ratiotGml2,missing-method}
\alias{initialize,ratiotGml2-method}
\alias{parameters,ratiotGml2-method}

\title{Class "ratiotGml2"  }
\description{ 
Ratio transformation as parameterized in Gating-ML 2.0 is defined by the following function:
\deqn{f(p1, p2, A, B, C) = A * (p1 - B) / (p2 - C)}
}

\section{Objects from the Class}{
Objects can be created by calls to the constructor

\code{ratiotGml2(p1, p2, A, B, C, transformationId)}
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{function} ~~ }
    \item{\code{numerator}:}{Object of class \code{"transformation"} 
      -flow parameter to be used as numerator in the transformation function }
    \item{\code{denominator}:}{Object of class \code{"transformation"} 
      -flow parameter to be used as denominator in the transformation function }    
    \item{\code{pA}:}{Object of class \code{numeric} constant A }
    \item{\code{pB}:}{Object of class \code{numeric} constant B }
    \item{\code{pC}:}{Object of class \code{numeric} constant C }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to reference the transformation}
  }
}

\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}

\section{Methods}{
No methods defined with class "ratiotGml2" in the signature.
}

\references{
Gating-ML 2.0: International Society for Advancement of Cytometry (ISAC) standard for 
representing gating descriptions in flow cytometry. \url{http://flowcyt.sourceforge.net/gating/20130122.pdf}
}

\author{Spidlen, J.}

\note{
The ratiotGml2 transformation object can be evaluated using the eval method by passing the data frame as an argument.
The transformed parameters are returned as matrix with one column. (See example below)
}

\seealso{
\code{\link{ratio}}, \code{\link{transform-class}}, \code{\link{transform}}
}

\examples{
myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
    package="flowCore"))
myRatioT <- ratiotGml2("FSC-H", "SSC-H", pA = 2, pB = 3, 
    pC = -10, transformationId = "myRatioT")
transOut <- eval(myRatioT)(exprs(myDataIn))
}

\keyword{classes}
