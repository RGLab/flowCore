\name{compensateActionItem-class}

\docType{class}

\alias{compensateActionItem-class}
\alias{compensateActionItem}
\alias{print,compensateActionItem-method}
\alias{Rm,compensateActionItem,workFlow,character-method}
\alias{show,compensateActionItem-method}

\title{Class "compensateActionItem"}

\description{
  
  Class and method to capture compensation operations in a flow
  cytometry workflow.

}

\section{Objects from the Class}{
  
  Objects should be created using the \code{add} method, which creates a
  \code{compensateActionItem} from a \code{\link{compensation}} object
  and directly assigns it to a \code{\link{workFlow}}. Alternatively,
  one can use the \code{compensateActionItem} constructor function for
  more programmatic access.
  
}

\section{Slots}{
  
  \describe{

    \item{\code{compensate}:}{Object of class
      \code{"fcCompensateReference"}. A reference to the
      \code{\link{compensation}} object that is used for the
      compensation operation. }
    
    \item{\code{ID}:}{Object of class \code{"character"}. A unique
      identifier for the \code{actionItem}. }
    
    \item{\code{name}:}{Object of class \code{"character"}. A more
      human-readable name}
    
    \item{\code{parentView}:}{Object of class
      \code{"fcViewReference"}. A reference to the parent
      \code{\link{view}} the \code{compensateActionItem} is applied on. }

    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment in the \code{\link{workFlow}}. }
   
  }
}

\section{Extends}{
Class \code{"\linkS4class{actionItem}"}, directly.
}


\usage{

compensateActionItem(ID = paste("compActionRef", guid(), sep = "_"),
name = paste("action", identifier(get(compensate)), sep = "_"),
parentView, compensate, workflow)

}

\arguments{

  \item{workflow}{ An object of class \code{\link{workFlow}} for which a
    view is to be created. }

  \item{ID}{ A unique identifier of the view, most likely created by
    using the internal \code{guid} function. }

  \item{name}{ A more human-readable name of the view. }
  
  \item{parentView, compensate}{ References to the parent
    \code{\link{view}} and \code{\link{compensation}} objects,
    respectively. }
  
}

\value{

  A reference to the \code{compensateActionItem} that is created inside the
  \code{\link{workFlow}} environment as a side effect of calling the
  \code{add} method.

  A \code{compensateActionItem} object for the constructor.
  
}


\section{Methods}{
  
  \describe{
    
    \item{print}{\code{signature(x = "compensateActionItem")}: Print
      details about the object. }
    
    \item{Rm}{\code{signature(symbol = "compensateActionItem", envir =
	"workFlow", subSymbol = "character")}: Remove a
	\code{compensateActionItem} from a \code{\link{workFlow}}. This
	method is recursive and will also remove all dependent
	\code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
    
    \item{show}{\code{signature(object = "compensateActionItem")}: Print
      details about the object. }
    
  }
}

\details{

  \code{compensateActionItems} provide a means to bind compensation
  operations in a workflow. Each \code{compensateActionItem} represents a
  single \code{\link{compensation}}.
  
}

\author{ Florian Hahne }

\seealso{
	
  \code{\linkS4class{workFlow}},
  \code{\linkS4class{actionItem}},
  \code{\linkS4class{gateActionItem}},
  \code{\linkS4class{transformActionItem}},
  \code{\linkS4class{view}}
  
}

\examples{
showClass("view")
}

\keyword{classes}
