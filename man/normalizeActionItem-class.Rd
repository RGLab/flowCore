\name{normalizeActionItem-class}

\docType{class}

\alias{normalizeActionItem-class}
\alias{normalizeActionItem}
\alias{print,normalizeActionItem-method}
\alias{Rm,normalizeActionItem,workFlow,character-method}
\alias{show,normalizeActionItem-method}

\title{Class "normalizeActionItem"}

\description{
  
  Class and method to capture normalization operations in a flow
  cytometry workflow.

}

\section{Objects from the Class}{
  
  Objects should be created using the \code{add} method, which creates a
  \code{normalizeActionItem} from a \code{\link{normalization}} object
  and directly assigns it to a \code{\link{workFlow}}. Alternatively,
  one can use the \code{normalizeActionItem} constructor function for
  more programmatic access.
  
}

\section{Slots}{
  
  \describe{

    \item{\code{normalization}:}{Object of class
      \code{"fcNormalizationReference"}. A reference to the
      \code{\link{normalization}} object that is used for the
      compensation operation. }
    
    \item{\code{ID}:}{Object of class \code{"character"}. A unique
      identifier for the \code{actionItem}. }
    
    \item{\code{name}:}{Object of class \code{"character"}. A more
      human-readable name}
    
    \item{\code{parentView}:}{Object of class
      \code{"fcViewReference"}. A reference to the parent
      \code{\link{view}} the \code{normalizeActionItem} is applied on. }

    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment in the \code{\link{workFlow}}. }
   
  }
}

\section{Extends}{
Class \code{"\linkS4class{actionItem}"}, directly.
}


\usage{

normalizeActionItem(ID = paste("normActionRef", guid(), sep = "_"),
name = paste("action", identifier(get(normalization)), sep = "_"),
parentView, normalization, workflow)

}

\arguments{

  \item{workflow}{ An object of class \code{\link{workFlow}} for which a
    view is to be created. }

  \item{ID}{ A unique identifier of the view, most likely created by
    using the internal \code{guid} function. }

  \item{name}{ A more human-readable name of the view. }
  
  \item{parentView, normalization}{ References to the parent
    \code{\link{view}} and \code{\link{normalization}} objects,
    respectively. }
  
}

\value{

  A reference to the \code{normalizeActionItem} that is created inside the
  \code{\link{workFlow}} environment as a side effect of calling the
  \code{add} method.

  A \code{normalizeActionItem} object for the constructor.
  
}


\section{Methods}{
  
  \describe{
    
    \item{print}{\code{signature(x = "normalizeActionItem")}: Print
      details about the object. }
    
    \item{Rm}{\code{signature(symbol = "normalizeActionItem", envir =
	"workFlow", subSymbol = "character")}: Remove a
	\code{normalizeActionItem} from a \code{\link{workFlow}}. This
	method is recursive and will also remove all dependent
	\code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
    
    \item{show}{\code{signature(object = "normalizeActionItem")}: Print
      details about the object. }
    
  }
}

\details{

  \code{normalizeActionItems} provide a means to bind normalization
  operations like warping in a workflow. Each
  \code{normalizeActionItem} represents a single
  \code{\link{normalization}}.
  
}

\author{ Florian Hahne }

\seealso{
	
  \code{\linkS4class{workFlow}},
  \code{\linkS4class{actionItem}},
  \code{\linkS4class{gateActionItem}},
  \code{\linkS4class{transformActionItem}},
  \code{\linkS4class{compensateActionItem}},
  \code{\linkS4class{view}}
  
}

\examples{
showClass("view")
}

\keyword{classes}
