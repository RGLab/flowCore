\name{compensation-class}

\docType{class}

\alias{compensation-class}
\alias{compensation}
\alias{compensate,flowFrame,compensation-method}
\alias{identifier,compensation-method}
\alias{identifier<-,compensation,character-method}
\alias{show,compensation-method}

\title{Class "compensation"}

\description{Class and methods to describe a compensation operation. }


\section{Objects from the Class}{
  
  Objects can be created by calls of the form \code{new("compensation",
    ...)} or via the contructor \code{compensation}.
  
}
\section{Slots}{
  
  \describe{
    
    \item{\code{spillover}:}{Object of class \code{"matrix"}. The
      spillover matrix. }
    
    \item{\code{invert}:}{Object of class \code{"logical"}. Switch to
      indicate whether the matrix in \code{spillover} is a spillover
      matrix or a compensation matrix (essentially the inverse)}
    
    \item{\code{compensationId}:}{Object of class \code{"character"}. A
      identifier for the object. }
    
  }
}



\usage{

compensation(spillover, invert = TRUE, compensationId = "default")

}

\arguments{

  \item{spillover}{ The spillover or compensation matrix. }

  \item{invert}{ A logical to swich between the two matrix
    representation. }

  \item{compensationId}{ The name of the object. }
  
}

\value{

  A \code{compensation} object.
  
}


\section{Methods}{
  
  \describe{
    
    \item{compensate}{\code{signature(x = "flowFrame", spillover =
	"compensation")}: Apply the compensation on a flowFrame. }
    
    \item{identifier, identifier<-}{\code{signature(object =
      "compensation")}: Accessor and replacement methods for the
      \code{compensationId} slot. }
    
    \item{show}{\code{signature(object = "compensation")}: Print details
      about the object. }
    
  }
}

\details{

  Essentially, the \code{compensation} class is just a wrapper around a matrix.
  
}


\author{ Florian Hahne }



\seealso{
	
  \code{\linkS4class{transform}},
  \code{\linkS4class{flowFrame}},
  
}

\examples{
showClass("view")
}

\keyword{classes}

