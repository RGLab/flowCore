\name{transformMap-class}
\docType{class}
\alias{transformMap-class}

\title{A class for mapping transforms between parameters}
\description{
This class provides a mapping between parameters and transformed parameters via a 
function.
}
\section{Objects from the Class}{
Objects of this type are not usually created by the user, except perhaps in special
circumstances. They are generally automatically created by the inline 
\code{\link[flowCore]{transform}} process during the creation of a \code{\link[flowCore]{transformFilter}}
}
\section{Slots}{
	 \describe{
    \item{\code{output}:}{Name of the transformed parameter}
    \item{\code{input}:}{Name of the parameter to transform}
    \item{\code{f}:}{Function used to accomplish the transform}
  }
}
\author{ B. Ellis }
\seealso{\code{\link[flowCore]{transform}}}
\examples{
new("transformMap",input="FSC-H",output="FSC-H",f=log)
}
\keyword{classes}
