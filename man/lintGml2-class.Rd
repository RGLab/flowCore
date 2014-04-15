\name{lintGml2-class}
\docType{class}
\alias{lintGml2-class}
\alias{lintGml2}
\alias{eval,lintGml2,missing-method}

\title{Class lintGml2}
\description{Linear transformation as parameterized in Gating-ML 2.0 is defined by the following function:
\deqn{f(parameter, T, A) = (parameter + A) / (T + A)}

This transformation provides a linear display that maps scale values from the \eqn{[-A, T]} interval to the \eqn{[0, 1]} interval. 
However, it is defined for all \eqn{x in R} including outside of the \eqn{[-A, T]} interval.
} 

\section{Objects from the Class}{
Objects can be created by calls to the constructor

\code{lintGml2(parameter, T, A, transformationId)}
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{function} ~~ }
    \item{\code{T}:}{Object of class \code{numeric} -positive constant (top of scale value) }
    \item{\code{A}:}{Object of class \code{numeric} -non-negative constant that is less than or equal
                     to T; it is determining the bottom end of the transformation }
    \item{\code{parameters}:}{Object of class \code{"transformation"} -flow parameter to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to reference the transformation}
  }
}

\section{Extends}{
Class \code{\linkS4class{singleParameterTransform}}, directly.
Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
}

\section{Methods}{
No methods defined with class lintGml2 in the signature.
}

\references{
Gating-ML 2.0: International Society for Advancement of Cytometry (ISAC) standard for 
representing gating descriptions in flow cytometry. \url{http://flowcyt.sourceforge.net/gating/20130122.pdf}
}

\author{Spidlen, J.}

\note{
The linear transformation object can be evaluated using the eval method by passing the data frame as an argument.
The transformed parameters are returned as a matrix with a single column. (See example below)  
}

\seealso{
\code{\link{linearTransform}}, \code{\link{transform-class}}, \code{\link{transform}}
}

\examples{
myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
    package="flowCore"))
myLinTr1 <- lintGml2(parameters = "FSC-H", T = 1000, A = 0, 
    transformationId="myLinTr1")
transOut <- eval(myLinTr1)(exprs(myDataIn))
}

\keyword{classes}
