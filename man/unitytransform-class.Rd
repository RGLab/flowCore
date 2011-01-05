\name{unitytransform-class}
\docType{class}
\alias{unitytransform-class}
\alias{unitytransform}
\alias{show,unitytransform-method}
\alias{eval,unitytransform,missing,missing-method}

\title{Class "unitytransform" }
\description{
Unity transform class transforms parameters names provided as characters into unity transform objects which can be evaluated to retrieve the corresponding columns from the
data frame	
}

\section{Objects from the Class}{
Objects can be created by calls to the constructor \code{unitytransform(parameters,transformationId)}.
	
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"character"} -flow parameters to be transformed }
    \item{\code{transformationId}:}{Object of class \code{"character"}-unique ID to
         reference the transformation  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}
\section{Methods}{
No methods defined with class "unitytransform" in the signature.
}

\author{ Gopalakrishnan N, F.Hahne }

\seealso{
dg1polynomial, ratio
}
\examples{
  dat <- read.FCS(system.file("extdata","0877408774.B08",
  package="flowCore"))
  un1<-unitytransform(c("FSC-H","SSC-H"),transformationId="un1")
  transOut<-eval(un1)(exprs(dat))
}
\keyword{classes}
