\name{workFlow-class}

\docType{class}

\alias{workFlow-class}
\alias{workFlow}
\alias{add}
\alias{add,workFlow,concreteFilter-method}
\alias{add,workFlow,filterList-method}
\alias{add,workFlow,transformList-method}
\alias{add,workFlow,compensation-method}
\alias{assign,ANY,ANY,missing,workFlow,missing,missing-method}
\alias{assign,missing,ANY,workFlow,missing,missing,missing-method}
\alias{assign,missing,ANY,missing,workFlow,missing,missing-method}
\alias{assign,character,ANY,workFlow,missing,missing,missing-method}
\alias{assign,fcReference,ANY,workFlow,missing,missing,missing-method}
\alias{[,workFlow,ANY-method}
\alias{[[,workFlow,ANY-method}
\alias{$,workFlow,character-method}
\alias{get,character,workFlow,missing,missing,missing-method}
\alias{get,character,missing,workFlow,missing,missing-method}
\alias{ls}
\alias{ls,workFlow,missing,missing,missing,missing-method}
\alias{ls,workFlow,missing,missing,missing,character-method}
\alias{mget}
\alias{mget,character,workFlow,missing,missing,missing-method}
\alias{names,workFlow-method}
\alias{plot,workFlow,missing-method}
\alias{Rm,character,workFlow,character-method}
\alias{show,workFlow-method}
\alias{views}
\alias{views,workFlow-method}
\alias{actions}
\alias{actions,workFlow-method}
\alias{alias,workFlow-method}
\alias{alias,environment-method}
\alias{summary,workFlow-method}
\alias{journal}
\alias{journal,workFlow-method}
\alias{journal,environment-method}
\alias{tree}
\alias{tree,workFlow-method}
\alias{tree,environment-method}
\alias{undo}

\title{Class "workFlow"}

\description{Class and methods to organize standard flow cytometry data
  analysis operation in a concise workflow.}


\section{Objects from the Class}{

  Objects should be created using the constructor \code{workFlow}, which
  takes a \code{\link{flowFrame}} or \code{\link{flowSet}} as only
  mandatory input and creates a basic view for that. 
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{name}:}{Object of class \code{"character"}. The name of
      the workFlow object. }
    
    \item{\code{tree}:}{Object of class \code{"fcTreeReference"}. A
      reference to the \code{\link[graph:graphNEL-class]{graphNEL}}
      objects representing the view structure of the workflow. }
    
    \item{\code{env}:}{Object of class \code{"environment"}. The
      evaluation environment for the workflow in which all objects will be
      stored. }

  }
}

\usage{

workFlow(data, name = "default", env = new.env(parent = emptyenv()))

undo(wf, n=1)

}

\arguments{

  \item{data}{ An object of class \code{\link{flowFrame}} or
    \code{\link{flowSet}} for which a basic \code{\link{view}} is
    created. }

  \item{name}{ A more human-readable name of the view. }
  
  \item{env}{ Object of class \code{environment}. The
    evaluation environment used for the \code{\link{workFlow}}. }

  \item{wf}{ Object of class \code{workFlow}. }

  
  \item{n}{ The number of operations to undo. }
  
}

\value{

  A \code{workFlow} object for the constructor

  Both \code{applyParentFilter} and \code{undo} are called for their
  side-effects.
  
}


\section{Methods}{
  
  \describe{
    
    \item{add}{\code{signature(wf = "workFlow", action =
	"concreteFilter")}: Create a new \code{\link{gateActionItem}} and
      \code{\link{gateView}} from a \code{\link{filter}} and assign those to
      the workflow. }

    \item{add}{\code{signature(wf = "workFlow", action =
	"filterList")}: Create a new \code{\link{gateActionItem}} and
      \code{\link{gateView}} from a \code{\link{filterList}}
	and assign those to the workflow. }

    \item{add}{\code{signature(wf = "workFlow", action =
	"transformList")}: Create a new \code{\link{transformActionItem}} and
      \code{\link{transformView}} from a \code{\link{transform}}
	and assign those to the workflow. }
    
    \item{add}{\code{signature(wf = "workFlow", action =
	"compensation")}: Create a new \code{\link{compensateActionItem}} and
      \code{\link{compensateView}} from a \code{\link{compensation}}
      and assign those to the workflow.}
    
    \item{assign}{\code{signature(x = "ANY", value = "ANY", pos =
	"missing", envir = "workFlow", inherits = "missing", immediate =
	"missing")}: Assign an object to the environment in the
      \code{workFlow} object and return a \code{\link{fcReference}} to
      it. The symbol for the object is created as a unique identifier. }

    \item{assign}{\code{signature(x = "missing", value = "ANY", pos =
	"workFlow", envir = "missing", inherits = "missing", immediate =
	"missing")}: see above }

    \item{assign}{\code{signature(x = "missing", value = "ANY", pos =
	"missing", envir = "workFlow", inherits = "missing", immediate =
	"missing")}: same as above, but provide custom symbol for the
      assignment.}
    
    \item{assign}{\code{signature(x = "character", value = "ANY", pos =
    "workFlow", envir = "missing", inherits = "missing", immediate =
    "missing")}: see above}

    \item{assign}{\code{signature(x = "fcReference", value = "ANY", pos
	= "workFlow", envir = "missing", inherits = "missing", immediate
	= "missing")}: same as above, but assign object using an
      existing \code{\link{fcReference}}. Note that assigning
      \code{NULL} esentially removes the original object. }

    \item{[}{\code{signature(x = "workFlow", i = "ANY")}: Cast a useful
      error message. }

    \item{[[}{\code{signature(x = "workFlow", i = "ANY")}: Treat the
      \code{workFlow} object as a regular environment. Essentially, this
      is equivalent to \code{get(x, i)}.  }

    \item{\$}{\code{signature(x = "workFlow", name = "character")}:
      Allow for list-like access. Note that completetion is only
      available for \code{\link[view-class]{views}} since all other
      objects in the environment are considered to be internal. }

    \item{get}{\code{signature(x = "character", pos = "workFlow", envir
	= "missing", mode = "missing", inherits = "missing")}: Get an
      object identified by symbol \code{x} from the environment in the
      \code{workFlow}. }
    
    \item{get}{\code{signature(x = "character", pos = "missing", envir =
	"workFlow", mode = "missing", inherits = "missing")}: see above }
    
    \item{ls}{\code{signature(name = "workFlow", pos = "missing", envir
	= "missing", all.names = "missing", pattern = "missing")}: List
      the content of the environment in the \code{workFlow}. }
    
    \item{ls}{\code{signature(name = "workFlow", pos = "missing", envir
	= "missing", all.names = "missing", pattern = "character")}: see
      above}
    
    \item{mget}{\code{signature(x = "character", envir = "workFlow",
	mode = "missing", ifnotfound = "missing", inherits =
	"missing")}: Get multiple objects identified by the symbols in
      \code{x} from the environment in the \code{workFlow}. }

    \item{names}{\code{signature(x = "workFlow")}: List the identifiers
      for all \code{\link[view-class]{views}} and
      \code{\link[actionItem-class]{actionItems}} in the
      \code{workFlow}. }

    \item{plot}{\code{signature(x = "workFlow", y = "missing")}: Plot
      the structure of the \code{workFlow} tree. }
    
    \item{Rm}{\code{signature(symbol = "character", envir = "workFlow",
	subSymbol = "character")}: Remove the object identified by the
      symbol \code{symbol} from the \code{workFlow}. }

    \item{undo}{\code{signature(wf = "workFlow", n = "numeric")}: Undo
      the last \code{n} operations on the  \code{workFlow}. }
    
    \item{show}{\code{signature(object = "workFlow")}: Print details
      about the object. }
    
    \item{summary}{\code{signature(object = "workFlow")}: Summarize a
      view in the \code{workFlow}. }

    \item{actions}{\code{signature(x = "workFlow")}: List the names of
      the \code{\link[actionItem-class]{actionItems}} in the
      \code{workFlow}. }

    \item{views}{\code{signature(x = "workFlow")}: List the names of
      only the \code{\link[view-class]{views}} in the
      \code{workFlow}. }

    \item{alias}{\code{signature(object = "workFlow")}: Return the alias
      table for the \code{workFlow}.}

    \item{alias}{\code{signature(object = "environment")}: Return the
      alias table from a generic environment. The method tries to find
      'fcAliasRef' among the object symbols in the environment. }
    
    \item{journal}{\code{signature(object = "workFlow")}: Return the
      journal for the \code{workFlow}.}

    \item{journal}{\code{signature(object = "environment")}: Return the
      journal from a generic environment. The method tries to find
      'fcJournalRef' among the object symbols in the environment. }

    \item{tree}{\code{signature(object = "workFlow")}: Return the
      tree of the \code{workFlow}.}
     
     \item{journal}{\code{signature(object = "environment")}: Return the
       tree from a generic environment. The method tries to find
       'fcTreeRef' among the object symbols in the environment. }
      
  }
}

\details{

  \code{workFlow} objects organize standard flow data analysis
  operations like gating, compensation and transformation in one single
  object. The user can interact with a \code{workFlow} object
  (e.g. adding operations, removing them, summarizing the results)
  without having to keep track of intermediate objects and names.

  The integral part of a \code{workFlow} is an evaluation environment
  which holds all objects that are created during the analysis. The
  structure of the whole workflow is a tree, where nodes represent
  \code{link{view}}s (or results of an operation) and edges represent
  \code{\link[actionItem-class]{actionItems}} (or the operations
  themselve).
  
}


\author{ Florian Hahne }



\seealso{
	
  \code{\linkS4class{view}},
  \code{\linkS4class{actionItem}}
  
}

\examples{
showClass("view")
}

\keyword{classes}
