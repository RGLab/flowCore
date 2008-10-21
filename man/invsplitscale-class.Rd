\name{invsplitscale-class}
\docType{class}
\alias{invsplitscale-class}
\alias{invsplitscale}
\alias{eval,invsplitscale,missing,missing-method}


\title{Class "invsplitscale"  }
\description{ The inverse split scale transformation is defined by the function 
\deqn{f(parameter,r,maxValue,transitionChannel)= \frac{(parameter-b)}{a}
      \;\;\;\;\;\;\;\;\;parameter<=t*a + b}
     
\deqn{ \frac{10^{parameter*\frac{d}{r}}}{c}   \;\;\;\;\;\;\;parameter > t*a+b }  

 where, 
 \deqn{b=\frac{transitionChannel}{2}}
\deqn{d=\frac{2*log_{10}(e)*r}{transitionChannel} + log_{10}(maxValue) }
\deqn{t=10^{log_{10}t}}
\deqn{a= \frac{transitionChannel}{2*t}}
\deqn{log_{10}ct=\frac{(a*t+b)*d}{r}}
\deqn{c=10^{log_{10}ct}}
}

\section{Objects from the Class}{
Objects can be created by calls to the constructor \newline
 \code{invsplitscale(parameters,r,maxValue,transitionChannel,transformationId)}
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{r}:}{Object of class \code{"numeric"} -a positive value indicating
     the range of the logarithmical part of the display  }
    \item{\code{maxValue}:}{Object of class \code{"numeric"}  -a positive value 
           indicating the maximum value the transformation is applied to}
    \item{\code{transitionChannel}:}{Object of class \code{"numeric"} -non negative 
          value that indicates where to split the linear vs. logarithmical transformation}
    \item{\code{parameters}:}{Object of class \code{"transformation"}  - flow parameter
        to be transformed}
    \item{\code{transformationId}:}{Object of class \code{"character"} -unique ID to
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
No methods defined with class "invsplitscale" in the signature.
}

\references{Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry  }

\author{Gopalakrishnan N,F.Hahne}

\note{ The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed
  parameters are returned as a matrix with a single column. (See example below) }

\seealso{splitscale

}

\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
  sp1<-invsplitscale("FSC-H",r=512,maxValue=2000,transitionChannel=512)
  transOut<-eval(sp1)(exprs(dat))
}
\keyword{classes}
