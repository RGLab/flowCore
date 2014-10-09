\name{logtGml2-class}
\docType{class}
\alias{logtGml2-class}
\alias{logtGml2}
\alias{eval,logtGml2,missing-method}

\title{Class logtGml2}
\description{Log transformation as parameterized in Gating-ML 2.0 is defined by the following function:
\deqn{bound(f, boundMin, boundMax) = max(min(f,boundMax),boundMin))}
where
\deqn{f(parameter, T, M) = (1/M) * log10(x/T) + 1}

This transformation provides a logarithmic display that maps scale values from the \eqn{(0, T]} interval to the 
\eqn{(-Inf, 1]} interval such that the data value T is mapped to 1 and M decades of data are mapped into the interval. 
Also, the limit for x going to 0 is -Inf.

In addition, if a boundary is defined by the boundMin and/or boundMax parameters, then the result of this
transformation is restricted to the [boundMin,boundMax] interval. Specifically, should the result of the f
function be less than boundMin, then let the result of this transformation be boundMin. Analogically, should the 
result of the f function be more than boundMax, then let the result of this transformation be boundMax.
The boundMin parameter shall not be greater than the boundMax parameter.
} 

\section{Objects from the Class}{
Objects can be created by calls to the constructor

\code{logtGml2(parameter, T, M, transformationId, boundMin, boundMax)}
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{function} ~~ }
    \item{\code{T}:}{Object of class \code{numeric} -positive constant (top of scale value) }
    \item{\code{M}:}{Object of class \code{numeric} -positive constant (number of decades) }
    \item{\code{parameters}:}{Object of class \code{"transformation"} -flow parameter to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to reference the transformation}
    \item{\code{boundMin}:}{Object of class \code{numeric} -lower bound of the transformation, default -Inf}
    \item{\code{boundMax}:}{Object of class \code{numeric} -upper bound of the transformation, default Inf}
  }
}

\section{Extends}{
Class \code{\linkS4class{singleParameterTransform}}, directly.
Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
}

\section{Methods}{
No methods defined with class logtGml2 in the signature.
}

\references{
Gating-ML 2.0: International Society for Advancement of Cytometry (ISAC) standard for 
representing gating descriptions in flow cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
}

\author{Spidlen, J.}

\note{
The log transformation object can be evaluated using the eval method by passing the data frame as an argument.
The transformed parameters are returned as a matrix with a single column. (See example below)  
}

\seealso{
\code{\link{logTransform}}, \code{\link{transform-class}}, \code{\link{transform}}
}

\examples{
myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
    package="flowCore"))
myLogTr1 <- logtGml2(parameters = "FSC-H", T = 1023, M = 4.5, 
    transformationId="myLogTr1")
transOut <- eval(myLogTr1)(exprs(myDataIn))
}

\keyword{classes}
