\name{filterSummary-class}
\docType{class}
\alias{filterSummary-class}
\alias{filterSummary}
\alias{summary}
\alias{[[,filterSummary,numeric-method}
\alias{[[,filterSummary,character-method}
\alias{$,filterSummary-method}
\alias{coerce,filterSummary,data.frame-method}
\alias{length,filterSummary-method}
\alias{names,filterSummary-method}
\alias{print,filterSummary-method}
\alias{show,filterSummary-method}
\alias{toTable}
\alias{toTable,filterSummary-method}

\title{Class "filterSummary"}

\description{ Class and methods to handle the summary information of a
  gating operation. }

\section{Objects from the Class}{
  
  Objects are created by calling \code{summary} on a
  \code{link{filterResult}} object. The user doesn't have to deal with
  manual object instantiation.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{name}:}{Object of class \code{"character"} The name(s) of
      the populations created in the filtering operation. For a
      \code{\link{logicalFilterResult}} this is just a single value; the
      name of the \code{link{filter}}. }
    
    \item{\code{true}:}{Object of class \code{"numeric"}. The number of
      events within the population(s). }
    
    \item{\code{count}:}{Object of class \code{"numeric"}. The total
      number of events in the gated \code{\link{flowFrame}}. }
    
    \item{\code{p}:}{Object of class \code{"numeric"} The percentage of
      cells in the population(s). }
  }
}

\usage{

summary(object, \dots)

}

\arguments{

  \item{object}{ An object inheriting from class
    \code{\link{filterResult}} which is to be summarized. }

  \item{\dots}{ Further arguments that are passed to the generic. }

}

\value{

  An object of class \code{filterSummary} for the \code{summary}
  constructor, a named list for the subsetting operators. The \code{$}
  operator returns a named vector of the respective value, where each
  named element corresponds to one sub-population.

}



\section{Methods}{
  \describe{
    
    \item{[[}{\code{signature(x = "filterSummary", i = "numeric")}:
      Subset the \code{filterSummary} to a single population. This only
      makes sense for
      \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResults}}.
      The output is a list of summary statistics. }
    
    \item{[[}{\code{signature(x = "filterSummary", i = "character")}:
      see above }
    
    \item{\$}{\code{signature(x = "filterSummary", name = "ANY")}: A
      list-like accessor to the slots and more. Valid values are
      \code{n} and \code{count} (those are identical), \code{true} and
      \code{in} (identical), \code{false} and \code{out} (identical),
      \code{name}, \code{p} and \code{q} (\code{1-p}).  }
    
    \item{coerce}{\code{signature(from = "filterSummary", to =
	"data.frame")}: Coerce object to \code{data.frame}. }
    
    \item{length}{\code{signature(x = "filterSummary")}: The number of
      populations in the \code{fitlerSummary}. }
    
    \item{names}{\code{signature(x = "filterSummary")}: The names of the
      populations in the \code{filterSummary}. }
    
    \item{print}{\code{signature(x = "filterSummary")}: Print details
      about the object. }
    
    \item{show}{\code{signature(object = "filterSummary")}: Print
      details about the object.}
    
    \item{toTable}{\code{signature(x = "filterSummary")}: Coerce object
      to \code{data.frame}. }
  }
}


\details{

  Calling \code{summary} on a \code{\link{filterResult}} object prints
  summary information on the screen, but also creates objects of class
  \code{filterSummary} for computational access.
  
}

\author{ Florian Hahne, Byron Ellis }

\seealso{
	
  \code{\linkS4class{filterResult}},
  \code{\linkS4class{logicalFilterResult}},
  \code{\linkS4class{multipleFilterResult}},
  \code{\linkS4class{flowFrame}}
   \code{\linkS4class{filterSummaryList}}
  
}

\examples{

## Loading example data, creating and applying a curv1Filter
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))
c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)
fres <- filter(dat, c1f)

## creating and showing the summary
summary(fres)
s <- summary(fres)

## subsetting
s[[1]]
s[["peak 2"]]

##accessing details
s$true
s$n
toTable(s)

}

\keyword{classes}



