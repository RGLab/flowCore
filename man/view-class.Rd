\name{view-class}

\docType{class}

\alias{view-class}
\alias{view}
\alias{action,view-method}
\alias{action}
\alias{alias}
\alias{alias,view-method}
\alias{Data,view-method}
\alias{Data}
\alias{names,view-method}
\alias{parent}
\alias{parent,view-method}
\alias{parent,NULL-method}
\alias{print,view-method}
\alias{Rm,view,workFlow,character-method}
\alias{show,view-method}
\alias{identifier,view-method}
\alias{identifier,NULL-method}
\alias{xyplot,formula,view-method}
\alias{xyplot,view,missing-method}

\title{Class "view" }

\description{ Class and method to capture the results of standard
  operations (called "views" here) in a flow cytometry workflow. }

\section{Objects from the Class}{

  Objects should be created using the constructor \code{view}, which
  also assigns the view to a \code{\link{workFlow}} object.
  
}
       
\section{Slots}{
  
  \describe{
    
    \item{\code{ID}:}{Object of class \code{"character"}. A unique
      identifier for the view. }
    
    \item{\code{name}:}{Object of class \code{"character"}. A more
      human-readable name}
    
    \item{\code{action}:}{Object of class \code{"fcActionReference"}. A
      reference to the
      \code{\link[flowCore:actionItem-class]{actionItem}} that generated
      the view. }

    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment in the \code{\link{workFlow}}. }

     \item{\code{alias}:}{Object of class \code{"fcAliasReference"}. A
      reference to the alias table.}
    
    \item{\code{data}:}{Object of class \code{"fcDataReference"} A
      reference to the data that is associated to the view. See
      \code{\link{gateView}} for details on copying and subsetting of the raw
      data in the context of gating.}
  }
}


\usage{

view(workflow, ID=paste("viewRef", guid(), sep="_"),
                 name="default", data, action)

parent(object)

Data(object)

action(object)

alias(object, \dots)

}

\arguments{

  \item{workflow}{ An object of class \code{\link{workFlow}} for which a
    view is to be created. }

  \item{object}{ An object of class \code{view} or one of its
    subclasses. }
    
  \item{ID}{ A unique identifier of the view, most likely created by
    using the internal \code{guid} function. }

  \item{name}{ A more human-readable name of the view. }
  
  \item{data, action}{ References to the data and
    \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively. }

  \item{\dots}{ Further arguments that get passed to the generic. }

}

\value{

  A reference to the view that is created inside the
  \code{\link{workFlow}} environment as a side effect of calling the
  constructor.

  The parent view (i.e., the view based on which the current view was
  created) for the parent method.
  
}

\section{Methods}{
  \describe{
    \item{action}{\code{signature(object = "view")}: Accessor for the
      \code{action} slot. Note that this returns the actual
      \code{\link[flowCore:actionItem-class]{actionItem}} object, i.e.,
  the reference gets resolved. } 

    \item{Data}{\code{signature(object = "view")}:  Accessor for the
      \code{data} slot. Note that this returns the actual data  object,
      i.e., the reference gets resolved. }
    
    \item{names}{\code{signature(x = "view")}: Accessor to the
      \code{name} slot. }

    \item{alias}{\code{signature(object = "view")}: Get the alias table
      from a \code{view}. }
    
    \item{parent}{\code{signature(object = "view")}: The parent view,
      i.e., the view based on which the current view was created. }
    
    \item{print}{\code{signature(x = "view")}: Print details about the
      object. }
    
    \item{Rm}{\code{signature(symbol = "view", envir = "workFlow",
	subSymbol = "character")}: Remove a \code{view} from a
      \code{\link{workFlow}}. This method is recursive and will also
      remove all dependent \code{views} and
      \code{\link[flowCore:actionItem-class]{actionItems}}. }
    
    \item{show}{\code{signature(object = "view")}: Print details about the
      object. }
      
    \item{xyplot}{\code{signature(x = "formula", data = "view")}: Plot
      the data underlying the view. }
    
    \item{xyplot}{\code{signature(x = "view", data = "missing")}: Plot
      the data underlying the view. }
  }
}

\details{

  \code{Views} provide a means to bind the results of standard
  operations on flow cytometry data in a workflow. Each view can be
  considered the outcome of one operation. There are more specific
  subclasses for the three possible types of operation:
  \code{\link{gateView}} for gating operations,
  \code{\link{transformView}} for transformations, and
  \code{\link{compensateView}} for compensation operations. See their
  documentation for details.
  
}


\author{ Florian Hahne }



\seealso{
	
  \code{\linkS4class{workFlow}},
  \code{\linkS4class{gateView}},
  \code{\linkS4class{transformView}},
  \code{\linkS4class{compensateView}},
  \code{\linkS4class{actionItem}}
  
}

\examples{
showClass("view")
}

\keyword{classes}
