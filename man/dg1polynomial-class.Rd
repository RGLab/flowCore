\name{dg1polynomial-class}
\docType{class}
\alias{dg1polynomial-class}
\alias{dg1polynomial}
\alias{eval,dg1polynomial,missing,missing-method}
\alias{initialize,dg1polynomial-method}
\alias{parameters<-,dg1polynomial,character-method}
\alias{parameters<-,dg1polynomial,parameters-method}


\title{Class "dg1polynomial"}

\description{ 
dg1polynomial allows for scaling,linear combination and translation within a single transformation defined by the function

\deqn{
f(parameter_1,...,parameter_n,a_1,...,a_n,b) = b + \Sigma_{i=1}^n 
       a_i*parameter_i
}
}


\section{Objects from the Class}{
Objects can be created by using the constructor \code{dg1polynomial(parameter,a,b,transformationId)}.
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"parameters"}
          -the flow parameters that are to be transformed}
    \item{\code{a}:}{Object of class \code{"numeric"} - coefficients of length equal to the number of flow parameters }
    \item{\code{b}:}{Object of class \code{"numeric"} - coefficient of length 1 that performs the translation }
    \item{\code{transformationId}:}{Object of class \code{"character"} unique ID to reference the transformation }
    }
}

\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}

\section{Methods}{
No methods defined with class "dg1polynomial" in the signature.
}

\references{
Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry V 1.5 
}

\author{ 
Gopalakrishnan N, F.Hahne
}

\note{ 
The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed parameters are returned as a matrix with a single column.(See example below)
}

\seealso{
ratio,quadratic,squareroot
}
\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  dg1<-dg1polynomial(c("FSC-H","SSC-H"),a=c(1,2),b=1,transformationId="dg1")
  transOut<-eval(dg1)(exprs(dat))
}

\keyword{classes}
