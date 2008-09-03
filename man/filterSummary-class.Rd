\name{filterSummary-class}
\docType{class}
\alias{filterSummary-class}
\alias{[[,filterSummary,character-method}
\alias{[[,filterSummary,numeric-method}
\alias{$,filterSummary,ANY-method}
\alias{show,filterSummary-method}
\alias{toTable}
\alias{toTable,filterSummary-method}
\alias{coerce,filterSummary,data.frame-method}
\alias{summary}
\alias{summary,subsetFilter-method}
\alias{names,filterSummary-method}
\alias{length,filterSummary-method}


\title{Class "filterSummary" }

\description{

  Class and methods to summarize a \code{\link{filterResult}} after
  applying a filter

}


\section{Objects from the Class}{
Objects are created internally by calling the \code{\link{summary}}
method on objects inheriting from class \code{\link{filterResult}}.
The user doesn't have to deal with manual object instantiation.
}

\section{Slots}{
  \describe{
    
    \item{\code{name}:}{Object of class \code{"character"} A vector of
      population names. In the case of a
      \code{\link{logicalFilterResult}}, this will only be one.}
       
    \item{\code{true}:}{Object of class \code{"numeric"} The number of
      events in the filter for each population.}
    
    \item{\code{count}:}{Object of class \code{"numeric"} The number of
      total events to which the \code{\link{filter}} was applied.}
    
    \item{\code{p}:}{Object of class \code{"numeric"} The percentage of
      events in the filter (as a ratio in \code{[0,1]}. }
    
  }
}


\usage{
summary(object, \dots)
}


\arguments{

  \item{object}{ Object inheriting from class
  \code{\link{filterResult}}. }
    
  \item{\dots}{ Further arguments, currently ignored. }
  
}


\value{

  An object of class \code{filterSummary} for the \code{summary}
  constructor, a named list for the subsetting operators. The \code{$}
  operator returns a named vector of the respective value, where each
  named element corresponds to one sub-population.
  
}


\section{Methods}{
  \describe{
    
    \item{[[}{\code{signature(x = "filterSummary", i = "character")}:
      Subsetting of  \code{filterSummary} objects for the case of
      multiple populations.}
   
     
    \item{\$}{\code{signature(x = "filterSummary", name = "ANY")}: A
      list-like accessor to the slots and more. Valid values are
      \code{n} and \code{count} (those are identical), \code{true} and
      \code{in} (identical), \code{false} and \code{out} (identical),
      \code{name}, \code{p} and \code{q} (\code{1-p}).}
  
    \item{show}{\code{signature(object = "filterSummary")}: A
      human-readable output}
    
    \item{toTable}{\code{signature(x = "filterSummary")}: Output as
      \code{data.frame} for more programmatic access.}
  }
}


\author{Florian Hahne, Byron Ellis}


\seealso{
  
  \code{\linkS4class{filterResult}},
  \code{\linkS4class{multipleFilterResult}}
  
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
