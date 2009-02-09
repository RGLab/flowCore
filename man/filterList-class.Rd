\name{filterResultList-class}

\docType{class}

\alias{filterList-class}
\alias{filterList}
\alias{show,filterList-method}
\alias{identifier,filterList-method}
\alias{identifier<-,filterList,character-method}

\title{Class "filterList" }

\description{Container for a list of \code{\link{filter}} objects. The
  class mainly exists for method dispatch. }

\section{Objects from the Class}{
  
  Objects are created from regular lists using the constructor \code{filterList}.

}
\section{Slots}{
  
  \describe{
    
    \item{\code{.Data}:}{Object of class \code{"list"}. The class
      directly extends \code{list}, and this slot holds the list data.}
    
    \item{\code{filterId}:}{Object of class \code{"character"}. The
      identifier for the object. }
    
  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
}


\usage{

filterList(x, filterId=identifier(x[[1]]))

}

\arguments{

  \item{x}{ A list of \code{\link{filter}} objects. }

  \item{filterId}{ The global identifier of the filter list. As default,
    we take the filterId of the first \code{filter} object in \code{x}. }
  
}

\value{

  A \code{filterList} object for the constructor

}



\section{Methods}{
  \describe{
    
    \item{show}{\code{signature(object = "filterList")}: Print
      details about the object. }
    
    \item{identifier, identifier<-}{\code{signature(object =
	"filteList")}: Accessor and replacement method for the object's
      filterId slot. }
    
  }
}


\author{Florian Hahne}

\seealso{
  \code{\link[flowCore]{filter}},
}

\examples{

f1 <- rectangleGate(FSC=c(100,200), filterId="testFilter")
f2 <- rectangleGate(FSC=c(200,400))
fl <- filterList(list(a=f1, b=f2))
fl
identifier(fl)

  }


\keyword{classes}
