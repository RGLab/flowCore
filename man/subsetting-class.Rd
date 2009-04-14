\name{subsetting-class}
\docType{class}
\alias{subsetting-class}
\alias{subsetting}
\alias{add,workFlow,subsetting-method}
\alias{identifier<-,subsetting,character-method}
\alias{identifier,subsetting-method}
\alias{show,subsetting-method}

\title{Class "subsetting"
}

\description{

 Class and methods to subset a a \code{flowSet}. This is only needed for
 method dispatch in the workFlow framework.

}

\section{Objects from the Class}{
  
  Objects should be created using the contructor
  \code{subsetting()}. See the \code{Usage} and \code{Arguments}
  sections for details.

}
\section{Slots}{
  \describe{
    
    \item{\code{subsettingId}:}{Object of class \code{"character"}. An
    identifier for the object.}

    \item{\code{indices}:}{Object of class
      \code{"numericOrCharacter"}. Indices into a \code{flowSet}.}

  }
}





\usage{

subsetting(indices, subsettingId="defaultSubsetting")


}

\arguments{

  \item{indices}{ Character or numeric vector of sample names. }

  
  \item{subsettingId}{ The identifier for the subsetting object. }
  
}

\value{

  A \code{subsetting} object.

}


\section{Methods}{
  
  \describe{
    
    \item{add}{\code{signature(wf = "workFlow", action =
	"subsetting")}: The constructor for the workFlow. }
    
    \item{identifier<-}{\code{signature(object = "subsetting", value
	= "character")}: Set method for the identifier slot. }
    
    \item{identifier}{\code{signature(object = "subsetting")}: Get
      method for the identifier slot. }
    
    \item{show}{\code{signature(object = "subsetting")}: Show details
      about the object. }
    
  }
    
}


\details{
  The class mainly existst for method dispatch in the workflow
  tools.
}

\author{F. Hahne}


\keyword{methods}  
\keyword{classes}
