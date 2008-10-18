\name{hyperlog-class}
\docType{class}
\alias{hyperlog-class}
\alias{hyperlog}
\alias{eval,hyperlog,missing,missing-method}

\title{Class "hyperlog"  }
\description{ Hyperlog transformation of a parameter is defined by 
the function 
\deqn{f(parameter,a,b)=root{EH(y,a,b)-parameter}}

where EH is a function defined by
\deqn{EH(y,a,b)= 10^{(\frac{y}{a})} + \frac{b*y}{a}-1      \;\;\;\;\;\;\;\;\;y>=0}
                            
\deqn{-10^{(\frac{-y}{a})} + \frac{b*y}{a}+1 \;\;\;\;\;\;\;y<0}

}
\section{Objects from the Class}{
Objects can be created by calls to the constuctor
\code{hyperlog(parameter,a,b,transformationId)}

}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} - numeric constant
treater than zero }
    \item{\code{b}:}{Object of class \code{"numeric"} numeric constant greater than zero }
    \item{\code{parameters}:}{Object of class \code{"transformation"} -flow parameter to be 
transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"} - unique ID to 
reference the transformation}
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}
\section{Methods}{
No methods defined with class "hyperlog" in the signature.
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed
parameters are returned as a matrix with a single column. (See example below) }

\seealso{

}
\examples{

  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  hlog1<-hyperlog("FSC-H",a=1,b=1,transformationId="hlog1")
  transOut<-eval(hlog1)(exprs(dat))
}
\keyword{classes}
