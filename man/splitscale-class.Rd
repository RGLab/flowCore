\name{splitscale-class}
\docType{class}
\alias{splitscale-class}
\alias{splitscale}
\alias{eval,splitscale,missing,missing-method}

\title{Class "splitscale"  }
\description{ The split scale transformation class defines a transformation that
 has a logarithmic scale at high values and a linear scale at low values. The transition points are chosen so that the slope of the  
trasformation is continuos at the transition points.

The split scale transformation is defined by the function

\deqn{f(parameter,r,maxValue,transitionChannel)= a*parameter+ b      \;\;\;\;\;\;\;\;\;parameter<=t}
     
\deqn{ log_{10}(c*parameter)*\frac{r}{d}  \;\;\;\;\;\;\;parameter > t }                       
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
 \code{splitscale(parameters,r,maxValue,transitionChannel,transformationId)}

}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{r}:}{Object of class \code{"numeric"}-a positive value indicating the range of the logarithmical 
                    part of the display }
    \item{\code{maxValue}:}{Object of class \code{"numeric"} -a positive value indicating the maximum value the transformation
      is applied to}
    \item{\code{transitionChannel}:}{Object of class \code{"numeric"} -non negative value that indicates where to 
          split the linear vs. logarithmical transformation}
    \item{\code{parameters}:}{Object of class \code{"transformation"} - flow parameter to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"}-unique ID to reference the transformation}
  }
}
\section{Extends}{
Class \code{"\linkS4class{singleParameterTransform}"}, directly.
Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
}
\section{Methods}{
No methods defined with class "splitscale" in the signature.
}
\references{ Gating-ML Candidate Recommendation for Gating Description in Flow Cytometry

}
\author{ Gopalakrishnan N, F.Hahne}
\note{  The transformation object can be evaluated using the eval method by passing the data frame as an argument.The transformed
  parameters are returned as a matrix with a single column. (See example below)}

\seealso{invsplitscale
	
}
\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
  sp1<-splitscale("FSC-H",r=768,maxValue=10000,transitionChannel=256)
  transOut<-eval(sp1)(exprs(dat))
}
\keyword{classes}
