\name{filterResultList-class}

\docType{class}

\alias{filterResultList-class}
\alias{filterResultList}
\alias{[,filterResultList,ANY-method}
\alias{[[,filterResultList,ANY-method}
\alias{names,filterResultList-method}
\alias{parameters,filterResultList-method}
\alias{show,filterResultList-method}
\alias{split,flowSet,filterResultList-method}
\alias{summary,filterResultList-method}

\title{Class "filterResultList" }

\description{Container to store the result of applying a \code{filter} on a
  \code{flowSet} object}

\section{Objects from the Class}{
  
  Objects are created by applying a \code{\link{filter}} on a
  \code{\link{flowSet}}. The user doesn't have to deal with manual
  object instantiation.

}
\section{Slots}{
  
  \describe{
    
    \item{\code{.Data}:}{Object of class \code{"list"}. The class
      directly extends \code{list}, and this slot holds the list data.}
    
    \item{\code{frameId}:}{Object of class \code{"character"} The IDs of
      the \code{\link[flowCore:flowFrame-class]{flowFrames}} in the filtered
      \code{\link{flowSet}}.}
  
    \item{\code{filterDetails}:}{Object of class \code{"list"}. Since
      \code{filterResultList} inherits from \code{\link{filterResult}},
      this slot has to be set. It contains only the input filter. }
    
    \item{\code{filterId}:}{Object of class \code{"character"}. The
      identifier for the object. }
    
  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{filterResult}"}, directly.
Class \code{"\linkS4class{concreteFilter}"}, by class "filterResult", distance 2.
Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 3.
}


\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "filterResultList", i = "ANY")}: Subset
      to \code{filterResultList}. }
    \item{[[}{\code{signature(x = "filterResultList", i = "ANY")}: Subset
      to individual \code{\link{filterResult}}. }
    
    \item{names}{\code{signature(x = "filterResultList")}: Accessor to
      the frameId slot. }
    
    \item{parameters}{\code{signature(object = "filterResultList")}:
      Return parameters on which data has been filtered. }
    
    \item{show}{\code{signature(object = "filterResultList")}: Print
      details about the object. }
    
    \item{split}{\code{signature(x = "flowSet", f =
	"filterResultList")}: Split a \code{\link{flowSet}} based on the
      results in the \code{filterResultlIst}. See \code{\link{split}}
      for details. }
    
    \item{summary}{\code{signature(object = "filterResultList")}:
      Summarize the filtering operation. This creates a
      \code{\link[flowCore:filterSummaryList-class]{filterSummaryList}}
      object. } 
    
  }
}


\author{Florian Hahne}

\seealso{
  \code{\linkS4class{filter}},
  \code{\linkS4class{filterResult}},
  \code{\linkS4class{logicalFilterResult}},
  \code{\linkS4class{multipleFilterResult}},
  \code{\linkS4class{randomFilterResult}}
}


\examples{

## Loading example data and creating a curv1Filter
data(GvHD)
dat <- GvHD[1:3]
c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)

## applying the filter
fres <- filter(dat, c1f)
fres

## subsetting the list
fres[[1]]
fres[1:2]

## details about the object
parameters(fres)
names(fres)
summary(fres)

## splitting based on the filterResults
split(dat, fres)
 
}

\keyword{classes}
