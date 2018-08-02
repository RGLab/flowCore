\name{fcReference-class}

\docType{class}

\alias{fcReference-class}
\alias{fcReference}
\alias{fcStructureReference-class}
\alias{fcTreeReference-class}
\alias{fcTreeReference}
\alias{fcAliasReference-class}
\alias{fcAliasReference}
\alias{fcDataReference-class}
\alias{fcDataReference}
\alias{fcActionReference-class}
\alias{fcActionReference}
\alias{fcViewReference-class}
\alias{fcViewReference}
\alias{fcFilterResultReference-class}
\alias{fcFilterResultReference}
\alias{fcFilterReference-class}
\alias{fcFilterReference}
\alias{fcCompensateReference-class}
\alias{fcCompensateReference}
\alias{fcTransformReference-class}
\alias{fcTransformReference}
\alias{fcNormalizationReference-class}
\alias{fcNormalizationReference}
\alias{fcNullReference-class}
\alias{fcSubsettingReference}
\alias{fcSubsettingReference-class}
\alias{assign}
\alias{fcNullReference}
\alias{get,fcReference,missing,missing-method}
\alias{get,fcNullReference,missing,missing-method}
\alias{identifier,fcReference-method}
\alias{isNull}
\alias{isNull,fcReference-method}
\alias{Rm}
\alias{Rm,fcReference,missing,character-method}
\alias{Rm,fcReference,workFlow,character-method}
\alias{Rm,fcNullReference,missing,character-method}
\alias{show,fcReference-method}
\alias{show,fcNullReference-method}




\title{Class "fcReference" and its subclasses}

\description{Classes and methods to provide reference-based semantics for flow
  cytometry workflows. }


\section{Extends}{

  \code{fcStructureReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.
  
  \code{fcTreeReference}:
  
  Class \code{"\linkS4class{fcStructureReference}"}, directly.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcStructureReference", distance 2.

  \code{fcAliasReference}:
  
  Class \code{"\linkS4class{fcStructureReference}"}, directly.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcStructureReference", distance 2.

  \code{fcDataReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

  \code{fcActionReference}:
  
  Class \code{"\linkS4class{fcStructureReference}"}, directly.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcStructureReference", distance 2.

  \code{fcViewReference}:
  
  Class \code{"\linkS4class{fcStructureReference}"}, directly.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcStructureReference", distance 2.

  \code{fcFilterResultReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

  \code{fcFilterReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

  \code{fcCompensateReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

  \code{fcTransformReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

   \code{fcNormalizationReference}:
  
  Class \code{"\linkS4class{fcReference}"}, directly.

  \code{fcNullReference}:
  
  Class \code{"\linkS4class{fcDataReference}"}, directly.
  Class \code{"\linkS4class{fcActionReference}"}, directly.
  Class \code{"\linkS4class{fcViewReference}"}, directly.
  Class \code{"\linkS4class{fcFilterResultReference}"}, directly.
  Class \code{"\linkS4class{fcFilterReference}"}, directly.
  Class \code{"\linkS4class{fcCompensateReference}"}, directly.
  Class \code{"\linkS4class{fcTransformReference}"}, directly.
  Class \code{"\linkS4class{fcNormalizationReference}"}, directly.
  Class \code{"\linkS4class{fcTreeReference}"}, directly.
  Class \code{"\linkS4class{fcAliasReference}"}, directly.
  Class \code{"\linkS4class{fcReference}"}, by class "fcDataReference",
  distance2.
  Class \code{"\linkS4class{fcStructureReference}"}, by class
  "fcActionReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcActionReference", distance 3.
  Class \code{"\linkS4class{fcStructureReference}"}, by class
  "fcViewReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class "fcViewReference",
  distance3.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcFilterResultReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcFilterReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcCompensateReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class
  "fcTransformReference", distance 2.
  Class \code{"\linkS4class{fcStructureReference}"}, by class
  "fcTreeReference", distance 2.
  Class \code{"\linkS4class{fcReference}"}, by class "fcTreeReference",
  distance 3.

}

\section{Objects from the Class}{
  
  Objects should be created via the \code{assign} constructor. Whenever
  an object is assigned to a \code{\link{workFlow}} using the
  \code{assign} method, an appropriate instance of class
  \code{fcReference} or one of its subclasses is returned. In addition,
  there are the usual constructor functions of same names as the classes
  that can be used for object instantiation without assignment. Note
  that this might lead to unresolvable references unless the object
  referenced to is available in the environment.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{ID}:}{Object of class \code{"character"} The name of the
      object in \code{env} referenced to. }
    
    \item{\code{env}:}{Object of class \code{"environment"} An
      environment that contains the referenced objects. Usually, this
      will be the environment that's part of a \code{\link{workFlow}}
      object.}
    
  }
}

\usage{

fcReference(ID=paste("genericRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcTreeReference(ID=paste("treeRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcAliasReference(ID=paste("aliasRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcDataReference(ID=paste("dataRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcActionReference(ID=paste("actionRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcViewReference(ID=paste("viewRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcFilterResultReference(ID=paste("fresRef", guid(), sep="_"),
    env=new.env(parent=emptyenv()))

fcFilterReference(ID=paste("filterRef", guid(), sep="_"), env=new.env(parent=emptyenv()))

fcCompensateReference(ID=paste("compRef",
                                           guid(), sep="_"),
                                  env=new.env(parent=emptyenv()))

fcNormalizationReference(ID=paste("normRef",
                                           guid(), sep="_"),
                                  env=new.env(parent=emptyenv()))

fcSubsettingReference(ID=paste("subRef",
                                           guid(), sep="_"),
                                  env=new.env(parent=emptyenv()))

fcTransformReference(ID=paste("transRef",
                                          guid(), sep="_"),
                                 env=new.env(parent=emptyenv()))

fcNullReference(...)


assign(x, value, pos = -1, envir = as.environment(pos), inherits = FALSE, 
    immediate = TRUE)

\S4method{get}{fcReference,missing,missing}(x, pos = -1, envir = as.environment(pos), mode = "any", 
    inherits = TRUE)

isNull(f)

Rm(symbol, envir, subSymbol, ...)

}

\arguments{

  \item{x, f, symbol}{ An object of class or inheriting from
    class \code{fcReference}.}

  \item{ID}{ The reference identifier.}
    
  \item{value}{ An arbitrary R object which is supposed to be assigned to
    the environment in the \code{\link{workFlow}} object and to which a
    reference is returned. }

  \item{env}{ An environment, usually within a \code{\link{workFlow}}
    object. }
  
  \item{pos, envir}{ Objects of class \code{\link{workFlow}}. }

  \item{inherits, immediate, mode, subSymbol, \dots}{Further arguments
    from the generics that are not used in this context. }

}


\value{

  An object of class \code{fcReference} or one of its subclasses for the
  \code{assign} constructor.

  The object referenced to for the \code{get} method.

  A character string of the object symbol for the \code{identifier}
  method.

  A logical scalar for the \code{isNull} method.
  
}

\section{Methods}{
  
  \describe{

    \item{get}{\code{signature(x = "fcReference", pos = "missing", envir
	= "missing", mode = "missing", inherits = "missing")}: Resolve
      the reference, i.e., get the object from the environment. }

    \item{get}{\code{signature(x = "fcNullReference", pos = "missing",
	envir = "missing", mode = "missing", inherits = "missing")}: Resolve
      the reference. This always returns \code{NULL}. }
    
      \item{identifier}{\code{signature(object = "fcReference")}: Return
      a character string of the object name. }
      
    \item{isNull}{\code{signature(f = "fcReference")}: Check whether a
      \code{fcReference} is a \code{NULL} reference. Note that this is
      different from a unresolvable reference. }
    
    \item{Rm}{\code{signature(symbol = "fcReference", envir = "missing",
	subSymbol = "character")}: Remove the object referenced to by a
      \code{fcReference} from its environment. The argument
      \code{subSymbol} will be automatically set by the generic and
      should never be provided by the user. }
    
    \item{Rm}{\code{signature(symbol = "fcReference", envir =
	"workFlow", subSymbol = "character")}: Remove the object referenced to by a
      \code{fcReference} from a \code{\link{workFlow}}. The argument
      \code{subSymbol} will be automatically set by the generic and
      should never be provided by the user. }

    \item{Rm}{\code{signature(symbol = "fcNullReference", envir =
	"missing", subSymbol = "character")}: Essentially, this doesn't do
      anything since there is no object referenced to. }
    
    \item{show}{\code{signature(object = "fcReference")}: Print details
      about the object. }
    
    \item{show}{\code{signature(object = "fcNullReference")}: Print details
      about the object. }
    
  }
}

\details{

  These classes provide references to objects within an R environment
  and allow for method dispatch based on the nature of the referenced
  object. The parent \code{fcReference} class is used for references to
  all R objects, unless there exists a more specific
  subclass. \code{fcTreeReference}, \code{fcViewReference}, and
  \code{fcActionReference} are used to reference to
  \code{\link[graph:graphNEL-class]{graphNEL}}, \code{\link{view}}, and
  \code{\link[flowCore:actionItem-class]{actionItem}} objects,
  respectively. \code{fcDataReference} should be used for
  \code{\link{flowFrame}} or \code{\link{flowSet}} objects, whereas
  \code{fcFilterResultReference}, \code{fcFilterReference},
  \code{fcTransformReference}, \code{fcCompensateReference}, and
  \code{fcNormalizationReference} link to \code{\link{filterResult}},
  \code{\link{filter}}, \code{\linkS4class{transform}} and
  \code{\linkS4class{compensation}} objects. \code{fsStructureReference}
  only exists to jointly dispatch on certain subgroups of references.
  
}



\author{ Florian Hahne }

\seealso{
  
  \code{\link{workFlow}}
  
}

\examples{
showClass("fcReference")
}

\keyword{classes}
