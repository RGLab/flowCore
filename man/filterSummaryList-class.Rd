\name{filterSummaryList-class}

\docType{class}

\alias{filterSummaryList-class}
\alias{filterSummaryList}
\alias{toTable,filterSummaryList-method}

\title{Class "filterSummaryList"}

\description{
  
  Class and methods to handle summary statistics for from filtering
  operations on whole \code{\link[flowSet-class]{flowSets}}.
  
}


\section{Objects from the Class}{
  
  Objects are created by calling \code{summary} on a
  \code{link{filterResultList}} object. The user doesn't have to deal
  with manual object instantiation.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{.Data}:}{Object of class \code{"list"}. The class
      directly extends \code{list}, and this slot holds the list data.}
       
  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
}

\section{Usage}{

summary(object, \dots)

}

\arguments{

  \item{object}{ An object of class
    \code{\link[filterResultList-class]{filterResultList}} which is to
    be summarized. }

  \item{\dots}{ Further arguments that are passed to the generic. }

}

\value{

  An object of class \code{filterSummaryList}.

}


\section{Methods}{
  \describe{
    
    \item{toTable}{\code{signature(x = "filterSummaryList")}: Coerce
      object to \code{data.frame}. Additional factors are added to
      indicate list items in the original object. }

  }
}

   

\details{

  Calling \code{summary} on a \code{\link{filterResult}} object prints
  summary information on the screen, but also creates objects of class
  \code{filterSummary} for computational access.
  
}

\author{ Florian Hahne }

\seealso{
	
  \code{\linkS4class{filterResult}},
  \code{\linkS4class{filterResultList}},
  \code{\linkS4class{logicalFilterResult}},
  \code{\linkS4class{multipleFilterResult}},
  \code{\linkS4class{flowFrame}}
  \code{\linkS4class{filterSummary}}
  
}

\examples{

## Loading example data, creating and applying a curv1Filter
data(GvHD)
dat <- GvHD[1:3]
c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)
fres <- filter(dat, c1f)

## creating and showing the summary
summary(fres)
s <- summary(fres)

## subsetting
s[[1]]

##accessing details
toTable(s)

}

\keyword{classes}



