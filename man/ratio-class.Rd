\name{ratio-class}
\docType{class}
\alias{ratio-class}
\alias{ratio}
\alias{eval,ratio,missing-method}
\alias{initialize,ratio-method}

\title{Class "ratio"  }
\description{ 
ratio transform calculates the ratio of two parameters defined by the function
\deqn{f(parameter_1,parameter_2)=\frac{parameter_1}{parameter_2}}
}

\section{Objects from the Class}{
Objects can be created by calls to the constructor
\code{ratio(parameter1,parameter2,transformationId)
}.
}

\section{Slots}{
\describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{numerator}:}{Object of class \code{"transformation"} -flow parameter
     to be transformed }
    \item{\code{denominator}:}{Object of class \code{"transformation"} -flow parameter
         to be transformed}
    \item{\code{transformationId}:}{Object of class \code{"character"} unique ID to
    reference the transformation }
}
}

\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}

\section{Methods}{
No methods defined with class "ratio" in the signature.
}

\references{ 
Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5 
}

\author{
Gopalakrishnan N, F.Hahne
}

\note{
The ratio transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed parameters are returned as matrix with one column. (See example below)
}

\seealso{
 dg1polynomial,quadratic,squareroot
}

\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  rat1<-ratio("FSC-H","SSC-H",transformationId="rat1")
  transOut<-eval(rat1)(exprs(dat))
}

\keyword{classes}
