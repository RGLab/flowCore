\name{hyperlogtGml2-class}
\docType{class}
\alias{hyperlogtGml2-class}
\alias{hyperlogtGml2}
\alias{eval,hyperlogtGml2,missing-method}

\title{Class hyperlogtGml2}
\description{Hyperlog transformation parameterized according to Gating-ML 2.0 is defined by the following function:
\deqn{hyperlog(x, T, W, M, A) = root(EH(y, T, W, M, A) - x)}
and \eqn{EH} is defined as:
\deqn{EH(y, T, W, M, A) = ae^{by} + cy - f}
where
\itemize{ 
  \item x is the value that is being transformed (an FCS dimension value). Typically, x is less than or 
    equal to T, although the transformation function is also defined for x greater than T.
  \item y is the result of the transformation.
  \item T is greater than zero and represents the top of scale value.
  \item M is greater than zero and represents the number of decades that the true logarithmic scale 
    approached at the high end of the Hyperlog scale would cover in the plot range.
  \item W is positive and not greater than half of M and represents the number of such decades 
    in the approximately linear region.
  \item A is the number of additional decades of negative data values to be included. A shall be 
    greater than or equal to \eqn{-W}, and less than or equal to \eqn{M - 2W} 
  \item root is a standard root finding algorithm (e.g., Newton's method) that finds y such as 
    \eqn{B(y, T, W, M, A) = x}.
}
and \eqn{a}, \eqn{b}, \eqn{c} and \eqn{f} are defined by means of \eqn{T}, \eqn{W}, \eqn{M}, 
\eqn{A}, \eqn{w}, \eqn{x0}, \eqn{x1}, \eqn{x2}, \eqn{e0}, \eqn{ca} and \eqn{fa} as:
\deqn{w = W/(M+A)}
\deqn{x2 = A/(M+A)}
\deqn{x1 = x2 + w}
\deqn{x0 = x2 + 2*w}
\deqn{b = (M + A)*ln(10)}
\deqn{e0 = e^{b*x0}}
\deqn{ca = e0/w}
\deqn{fa = e^{b*x1} + ca*x1}
\deqn{a = T / (e^b + ca - fa)}
\deqn{c = ca * a}
\deqn{f = fa * a}

} 

\section{Objects from the Class}{
Objects can be created by calls to the constructor

\code{hyperlogtGml2(parameter, T, M, W, A, transformationId)}
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{function} ~~ }
    \item{\code{T}:}{Object of class \code{numeric} -positive constant (top of scale value) }
    \item{\code{M}:}{Object of class \code{numeric} -positive constant (desired number of decades) }
    \item{\code{W}:}{Object of class \code{numeric} -positive constant that is not greater than half of M
      (the number of such decades in the approximately linear region) }
    \item{\code{A}:}{Object of class \code{numeric} -a constant that is greater than or equal to -W, and also
      less than or equal to M-2W. (A represents the number of additional decades of negative data values to 
      be included.) }
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
No methods defined with class hyperlogtGml2 in the signature.
}

\references{
Gating-ML 2.0: International Society for Advancement of Cytometry (ISAC) standard for 
representing gating descriptions in flow cytometry. \url{http://flowcyt.sourceforge.net/gating/20130122.pdf}
}

\author{Spidlen, J., Moore, W.}

\note{
That \code{hyperlogtGml2} transformation brings "reasonable" data values to the scale of \eqn{[0,1]}. 
The transformation is somewhat similar to \code{\link{logicletGml2}}. (See Gating-ML 2.0 for detailed comparison) 

The hyperlog transformation object can be evaluated using the eval method by passing the data frame as an argument.
The transformed parameters are returned as a matrix with a single column. (See example below)  
}

\seealso{
\code{\link{hyperlog}}, \code{\link{logicleTransform}}, \code{\link{transform-class}}, \code{\link{transform}}
}

\examples{
myDataIn  <- read.FCS(system.file("extdata", "0877408774.B08", 
    package="flowCore"))
myHyperLg <- hyperlogtGml2(parameters = "FSC-H", T = 1023, M = 4.5, 
    W = 0.5, A = 0, transformationId="myHyperLg")
transOut  <- eval(myHyperLg)(exprs(myDataIn))
}

\keyword{classes}
