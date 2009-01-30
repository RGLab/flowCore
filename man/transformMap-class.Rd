\name{transformMap-class}
\docType{class}
\alias{transformMap-class}
\alias{transformMap}
\alias{show,transformMap-method}

\title{A class for mapping transforms between parameters}

\description{
  
  This class provides a mapping between parameters and transformed
  parameters via a function.

}

\section{Objects from the Class}{
  
  Objects of this type are not usually created by the user, except
  perhaps in special circumstances. They are generally automatically
  created by the inline \code{\link[flowCore]{transform}} process during
  the creation of a \code{\link{transformFilter}}, or by a call to the
  \code{\link{transformList}} constructor.

}

\section{Slots}{
  \describe{
    \item{\code{output}:}{Name of the transformed parameter}
    \item{\code{input}:}{Name of the parameter to transform}
    \item{\code{f}:}{Function used to accomplish the transform}
  }
}

\section{Methods}{
\describe{
  \item{show}{\code{signature(object = "transformList")}: Print details
    about the object. }
}
}

\author{ B. Ellis, F. Hahne }
    
\seealso{

  \code{\link{transform}}, \code{\link{transformList}}

}

\examples{

new("transformMap", input="FSC-H", output="FSC-H", f=log)

}

\keyword{classes}
