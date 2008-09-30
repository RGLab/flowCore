\name{actionItem-class}

\docType{class}

\alias{actionItem-class}
\alias{identifier,actionItem-method}
\alias{names,actionItem-method}
\alias{alias,actionItem-method}
\alias{parent,actionItem-method}
\alias{Rm,actionItem,workFlow,character-method}

\title{Class "actionItem"}


\description{ Class and method to capture standard operations in a flow
cytometry workflow. }

\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}

\section{Slots}{
  
  \describe{

    \item{\code{ID}:}{Object of class \code{"character"}. A unique
      identifier for the \code{actionItem}. }
    
    \item{\code{name}:}{Object of class \code{"character"}. A more
      human-readable name}
    
    \item{\code{parentView}:}{Object of class \code{"fcViewReference"}. A
      reference to the parent \code{\link{view}} the \code{actionItem}
      is applied on. }

    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment in the \code{\link{workFlow}}. }

  }
}


\section{Methods}{
  
  \describe{
    
    \item{identifier}{\code{signature(object = "actionItem")}: Accessor
      for the \code{ID} slot. }
    
    \item{names}{\code{signature(x = "actionItem")}: Accessor
      for the \code{name} slot.}
    
    \item{parent}{\code{signature(object = "actionItem")}: Accessor for
      the \code{parentView} slot. Note that the reference is resolved,
      i.e., the \code{\link{view}} object is returned. }
    
    \item{alias}{\code{signature(object = "actionItem")}: Get the alias table
      from a \code{actionItem}. }
    
    \item{Rm}{\code{signature(symbol = "actionItem", envir = "workFlow",
	subSymbol = "character")}:  Remove a \code{actionItem} from a
      \code{\link{workFlow}}. This method is recursive and will also
      remove all dependent \code{\link[view-class]{views}} and
      \code{actionItems}. }

  }
}

\details{

  \code{actionItems} provide a means to bind standard operations on flow
  cytometry data in a workflow. Usually, the user doesn't have to create
  these objects, instead they will be automatically created when
  allpying one of the standard operations (gating, transformation,
  compensation) to a \code{\link{workFlow}} object. Each
  \code{actionItem} creates one or several new
  \code{\link[view-class]{views}}, which again can be the basis to apply
  further operations. One can conceptualize \code{actionItems} being the
  edges in the workflow tree connecting \code{\link[view-class]{views}},
  which are the nodes of the tree. There are more specific subclasses
  for the three possible types of operation:
  \code{\link{gateActionItem}} for gating oprations,
  \code{\link{transformActionItem}} for transformations, and
  \code{\link{compensateActionItem}} for compensation operations. See
  their documentation for details.
  
}



\author{ Florian Hahne }



\seealso{
	
  \code{\linkS4class{workFlow}},
  \code{\linkS4class{gateActionItem}},
  \code{\linkS4class{transformActionItem}},
  \code{\linkS4class{compensateActionItem}},
  \code{\linkS4class{view}}
  
}

\examples{
showClass("view")
}

\keyword{classes}


