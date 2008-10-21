\name{normalization-class}
\docType{class}
\alias{normalization-class}
\alias{normalization}
\alias{normalize}
\alias{add,workFlow,normalization-method}
\alias{identifier<-,normalization,character-method}
\alias{identifier,normalization-method}
\alias{normalize,flowSet,normalization-method}
\alias{parameters,normalization-method}

\title{Class "normalization"
}

\description{

 Class and methods to normalize a a \code{flowSet} using a potentially
 complex normalization function.

}

\section{Objects from the Class}{
  
  Objects should be created using the contructor
  \code{normalization()}. See the \code{Usage} and \code{Arguments}
  sections for details.

}
\section{Slots}{
  \describe{
    
    \item{\code{parameters}:}{Object of class \code{"character"}. The
      flow parameters that are supposed to be normalized by the
      normalization function. }
    
    \item{\code{normalizationId}:}{Object of class \code{"character"}. An
    identifier for the object.
    }

    \item{\code{normFunction}:}{Object of class \code{"function"} The
      normalization function. It has to take two mandatory arguments:
      \code{x}, the \code{flowSet}, and \code{parameters}, a character
      of parameter names that are to be normalized by the
      function. Additional arguments have to be passed in vio
      \code{arguments}.}
    
    \item{\code{arguments}:}{Object of class \code{"list"} A names list
      of additional arguments. Can be \code{NULL}. }

  }
}





\usage{

normalization(parameters, normalizationId="defaultNormalization",
normFunction, arguments=list())

normalize(data, x)

}

\arguments{

  \item{parameters}{ Character vector of parameter names. }

  
  \item{normalizationId}{ The identifier for the normalization object. }

  \item{x}{ An object of class \code{\linkS4class{flowSet}}. }

  \item{normFunction} {The normalization function}
  
  \item{arguments} {The list of additional arguments to
    \code{normFunction}}
    
  \item{data} {The \code{flowSet} to normalize. }
  
}

\value{

  A \code{normalization} object for the constructor.

  A \code{\linkS4class{flowSet}} for the \code{normalize} methods.
  
}


\section{Methods}{
  
  \describe{
    
    \item{add}{\code{signature(wf = "workFlow", action =
	"normalization")}: The constructor for the workFlow. }
    
    \item{identifier<-}{\code{signature(object = "normalization", value
	= "character")}: Set method for the identifier slot. }
    
    \item{identifier}{\code{signature(object = "normalization")}: Get
      method for the identifier slot. }
    
    \item{normalize}{\code{signature(data = "flowSet", x =
	"normalization")}: Apply a normalization to a \code{\linkS4class{flowSet}}. }
    
    \item{parameters}{\code{signature(object = "normalization")}: The
	 more generic constructor. }
	 }
}


\details{
  Data normalization of a \code{flowSet} is a rather fuzzy concept and
  the class mainly existst for method dispatch in the workflow
  tools. The idea is to have a rather general function that takes a
  \code{flowSet} and a list of parameter names as input and applies any
  kind of normalization to the respective data columns. The output of
  the function has to be a \code{flowSet} again. Although we don't
  formally check for it, the dimensions of the input and of the output
  set should remain the same. Additional arguments may be passed to the
  normalization function via the \code{arguments} list. Internally we
  evaluate the function using \code{\link{do.call}} and one should check
  its documentation for details.

  Currently, the most prominent example for a normalization function is
  warping, as provided by the \code{flowStats} package.


}

\author{F. Hahne}


\keyword{methods}  
\keyword{classes}
