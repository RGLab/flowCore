\name{subsettingActionItem-class}

\docType{class}

\alias{subsettingActionItem-class}
\alias{subsettingActionItem}
\alias{print,subsettingActionItem-method}
\alias{Rm,subsettingActionItem,workFlow,character-method}
\alias{show,subsettingActionItem-method}

\title{Class "subsettingActionItem"}

\description{
  
  Class and method to capture subsetting operations in a flow
  cytometry workflow.

}

\section{Objects from the Class}{
  
  Objects should be created using the \code{add} method, which creates a
  \code{subsettingActionItem} from a \code{\link{normalization}} object
  and directly assigns it to a \code{\link{workFlow}}. Alternatively,
  one can use the \code{subsettingActionItem} constructor function for
  more programmatic access.
  
}

\section{Slots}{
  
  \describe{

    \item{\code{subsetting}:}{Object of class
      \code{"fcSubsettingReference"}. A reference to the
      \code{\link{subsetting}} object that is used for the
      operation. }
    
    \item{\code{ID}:}{Object of class \code{"character"}. A unique
      identifier for the \code{actionItem}. }
    
    \item{\code{name}:}{Object of class \code{"character"}. A more
      human-readable name}
    
    \item{\code{parentView}:}{Object of class
      \code{"fcViewReference"}. A reference to the parent
      \code{\link{view}} the \code{subsettingActionItem} is applied on. }

    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment in the \code{\link{workFlow}}. }
   
  }
}

\section{Extends}{
Class \code{"\linkS4class{actionItem}"}, directly.
}


\usage{

subsettingActionItem(ID = paste("subActionRef", guid(), sep = "_"),
name = paste("action", identifier(get(subsetting)), sep = "_"),
parentView, subsetting, workflow)

}

\arguments{

  \item{workflow}{ An object of class \code{\link{workFlow}} for which a
    view is to be created. }

  \item{ID}{ A unique identifier of the view, most likely created by
    using the internal \code{guid} function. }

  \item{name}{ A more human-readable name of the view. }
  
  \item{parentView, subsetting}{ References to the parent
    \code{\link{view}} and \code{\link{subsetting}} objects,
    respectively. }
  
}

\value{

  A reference to the \code{subsettingActionItem} that is created inside the
  \code{\link{workFlow}} environment as a side effect of calling the
  \code{add} method.

  A \code{subsettingActionItem} object for the constructor.
  
}


\section{Methods}{
  
  \describe{
    
    \item{print}{\code{signature(x = "subsettingActionItem")}: Print
      details about the object. }
    
    \item{Rm}{\code{signature(symbol = "subsettingActionItem", envir =
	"workFlow", subSymbol = "character")}: Remove a
	\code{subsettingActionItem} from a \code{\link{workFlow}}. This
	method is recursive and will also remove all dependent
	\code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
    
    \item{show}{\code{signature(object = "subsettingActionItem")}: Print
      details about the object. }
    
  }
}

\details{

  \code{subsettingActionItems} provide a means to bind subsetting
  operations in a workflow. Each
  \code{subsettingActionItem} represents a single
  \code{\link{subsetting}}.
  
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
