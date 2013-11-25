\name{flowSet-class}
\docType{class}
\alias{flowSet-class}
\alias{flowSet}
\alias{[,flowSet-method}
\alias{[,flowSet,ANY-method}
\alias{$,flowSet-method}
\alias{[[,flowSet-method}
\alias{[[,flowSet,ANY-method}
\alias{[[<-,flowSet-method}
\alias{[[<-,flowSet,ANY,ANY,flowFrame-method}
\alias{[[<-,flowFrame-method}
\alias{fsApply,flowSet-method}
\alias{show,flowSet-method}
\alias{length,flowSet-method}
\alias{colnames,flowSet-method}
\alias{colnames<-,flowSet-method}
\alias{identifier,flowSet-method}
\alias{identifier<-,flowSet,ANY-method}
\alias{compensate,flowSet,ANY-method}
\alias{sampleNames,flowSet-method}
\alias{sampleNames<-,flowSet,ANY-method}
\alias{phenoData,flowSet-method}
\alias{phenoData<-,flowSet,ANY-method}
\alias{phenoData<-,flowSet,phenoData-method}
\alias{pData,flowSet-method}
\alias{pData<-,flowSet,data.frame-method}
\alias{plot,flowSet,ANY-method}
\alias{plot,flowSet-method}
\alias{varLabels,flowSet-method}
\alias{varLabels<-,flowSet-method}
\alias{varLabels<-,flowSet,ANY-method}
\alias{varMetadata,flowSet-method}
\alias{varMetadata<-,flowSet,ANY-method}
\alias{compensate,flowSet,matrix-method}
\alias{transform,flowSet-method}
\alias{rbind2,flowSet,missing}
\alias{rbind2,flowSet,flowSet-method}
\alias{rbind2,flowSet,flowSet,missing-method}
\alias{rbind2,flowSet,flowFrame-method}
\alias{rbind2,flowFrame,flowSet-method}
\alias{rbind2,flowSet,missing-method} 
\alias{summary,flowSet-method}



\title{'flowSet': a class for storing flow cytometry raw data from
quantitative cell-based assays}

\description{This class is a container for a set of
  \code{\linkS4class{flowFrame}} objects}


\section{Creating Objects}{

  Objects can be created using\cr
  \code{  new('flowSet',}\cr
  \code{    frames = ....,   # environment with flowFrames}\cr
  \code{    phenoData = .... # object of class AnnotatedDataFrame}\cr
  \code{    colnames = ....  # object of class character}\cr
  \code{  )}\cr

  or via the constructor \code{flowSet}, which takes arbitrary numbers
  of flowFrames, either as a list or directly as arguments, along with
  an optional
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} for
  the \code{phenoData} slot and a \code{character} scalar for the
  \code{name} by which the object can be referenced.

  \code{flowSet(..., phenoData)}

  Alternatively, \code{flowSets} can be coerced from \code{list} and
  \code{environment} objects.

  \code{as(list("A"=frameA,"B"=frameB),"flowSet")}

  The safest and easiest way to create \code{flowSet}s directly from
  \acronym{FCS} files is via the \code{\link{read.flowSet}} function,
  and there are alternative ways to specify the files to read. See the
  separate documentation for details.

}

\section{Slots}{
  
  \describe{
    
    \item{\code{frames}:}{An \code{\link[base:environment]{environment}}
      containing one or more \code{\linkS4class{flowFrame}} objects.}

    \item{\code{phenoData}:}{A
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
      containing the phenotyoic data for the whole data set. Each row
      corresponds to one of the \code{\linkS4class{flowFrame}}s in the
      \code{frames} slot.  The \code{sampleNames} of \code{phenoData}
      (see below) must match the names of the
      \code{\linkS4class{flowFrame}} in the \code{frames} environment.
      }

    \item{\code{colnames}:}{A \code{character} object with the (common)
      column names of all the data matrices in the
      \code{\linkS4class{flowFrame}}s.}
    
  }
}


\section{Methods}{
  
  \describe{

    \item{[, [[}{Subsetting. \code{x[i]} where \code{i} is a scalar,
      returns a \code{flowSet} object, and \code{x[[i]]} a
      \code{\linkS4class{flowFrame}} object. In this respect the
      semantics are similar to the behavior of the subsetting operators
      for lists. \code{x[i, j]} returns a \code{flowSet} for which the
      parameters of each \code{\linkS4class{flowFrame}} have been subset
      according to \code{j}, \code{x[[i,j]]} returns the subset of a
      single \code{\linkS4class{flowFrame}} for all parameters in
      \code{j}. Similar to data frames, valid values for \code{i} and
      \code{j} are logicals, integers and characters.
      
      \emph{Usage:}

      \code{   flowSet[i]}
      
      \code{   flowSet[i,j]}
      
      \code{   flowSet[[i]]}
      
    }
    
    \item{\$}{Subsetting by frame name. This will return a single
      \code{\linkS4class{flowFrame}} object. Note that names may have to
      be quoted if they are no valid R symbols
      (e.g. \code{flowSet$"sample 1"}}
    
    \item{colnames, colnames<-}{Extract or replace the \code{colnames}
      slot.
      
      \emph{Usage:}
      
      \code{   colnames(flowSet)}
      
      \code{   colnames(flowSet) <- value}
      
    }

    \item{identifier, identifier<-}{Extract or replace the \code{name}
      item from the environment.
      
      \emph{Usage:}
      
      \code{   identifier(flowSet)}
      
      \code{   identifier(flowSet) <- value}
      
    }
    
    
    \item{phenoData, phenoData<-}{Extract or replace the
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
      from the \code{phenoData} slot.
      
      \emph{Usage:}
      
      \code{   phenoData(flowSet)}
      
      \code{   phenoData(flowSet) <- value}
      
    }
    
    \item{pData, pData<-}{Extract or replace the data frame (or columns
      thereof) containing actual phenotypic information from the
      \code{phenoData} slot.
      
      \emph{Usage:}
      
      \code{   pData(flowSet)}

      \code{   pData(flowSet)$someColumn <- value}
      
    }
    
    \item{varLabels, varLabels<-}{ Extract and set varLabels in the
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
      of the \code{phenoData} slot.
      
      \emph{Usage:}
      
      \code{   varLabels(flowSet)}
      
      \code{   varLabels(flowSet) <- value}
      
    }
    
    \item{sampleNames}{Extract and replace sample names from the
      \code{phenoData} object. Sample names correspond to frame
      identifiers, and replacing them will also replace the \code{GUID}
      slot for each frame. Note that \code{sampleName} need to be
      unique.
      
      \emph{Usage:}
      
      \code{   sampleNames(flowSet)}
      
      \code{   sampleNames(flowSet) <- value}
      
    }
    
    \item{keyword}{Extract or replace keywords specified in a character
      vector or a list from the \code{description} slot of each
      frame. See \code{\link{keyword}} for details.
      
      \emph{Usage:}
      
      \code{   keyword(flowSet, list(keywords))}
      
      \code{   keyword(flowSet, keywords)}

      \code{   keyword(flowSet) <- list(foo="bar") }
            
    }
    
    \item{length}{number of \code{\linkS4class{flowFrame}} objects in
      the set.
      
      \emph{Usage:}
      
      \code{   length(flowSet)}
      
    }
    
    \item{show}{display object summary.}

     \item{summary}{Return descriptive statistical summary (min, max,
      mean and quantile) for each channel of each
      \code{\linkS4class{flowFrame}}

      \emph{Usage:}
      
      \code{   summary(flowSet)}
      
    }
    
    
    \item{fsApply}{Apply a function on all frames in a \code{flowSet}
      object. Similar to \code{\link{sapply}}, but with additional
      parameters. See separate documentation for details.
      
      \emph{Usage:}
      
      \code{   fsApply(flowSet, function, ...)}
      
      \code{   fsApply(flowSet, function, use.exprs=TRUE, ...)}
      
    }
    
    \item{compensate}{Apply a compensation matrix on all frames in a
      \code{flowSet} object. See separate documentation for details.
      
      \emph{Usage:}
      
      \code{   compensate(flowSet, matrix)}
      
    }
    
    \item{transform}{Apply a transformation function on all frames of a
      \code{flowSet} object. See separate documentation for details.
      
      \emph{Usage:}
      
      \code{   transform(flowSet, ...)}
      
    }
    
    \item{filter}{Apply a filter object on a \code{flowSet}
      object. There are methods for \code{\linkS4class{filter}}s,
      \code{\link{filterSet}}s and lists of filters. The latter has to
      be a named list, where names of the list items are matching
      sampleNames of the \code{flowSet}. See \code{\linkS4class{filter}}
      for details.
      
      \emph{Usage:}
      
      \code{   filter(flowSet, filter)}
      
      \code{   filter(flowSet, list(filters))}
      
    }
    
    \item{split}{Split all \code{flowSet} objects according to a
      \code{\link{filter}}, \code{\link{filterResult}} or a list of such
      objects, where the length of the list has to be the same as the
      length of the \code{flowSet}. This returns a list of
      \code{\linkS4class{flowFrame}}s or an object of class
      \code{flowSet} if the \code{flowSet} argument is set to
      \code{TRUE}. Alternatively, a \code{flowSet} can be split into
      separate subsets according to a factor (or any vector that can be
      coerced into factors), similar to the behaviour of
      \code{\link[base]{split}} for lists. This will return a list of
      \code{flowSet}s. See \code{\link{split}} for details.
      
      \emph{Usage:}
      
      \code{   split(flowSet, filter)}
      
      \code{   split(flowSet, filterResult)}
      
      \code{   split(flowSet, list(filters))}
      
      \code{   split(flowSet, factor)}
      
    }
    
    \item{Subset}{Returns a \code{flowSet} of
      \code{\linkS4class{flowFrame}}s that have been subset according
      to a \code{\linkS4class{filter}} or
      \code{\linkS4class{filterResult}}, or according to a list of such
      items of equal length as the \code{flowSet}.
      
      \emph{Usage:}
      
      \code{   Subset(flowSet, filter)}
      
      \code{   Subset(flowSet, filterResult)}
      
      \code{   Subset(flowSet, list(filters))}
      
    }
    
    
    \item{rbind2}{Combine two \code{flowSet} objects, or one
      \code{flowSet} and one \code{\linkS4class{flowFrame}} object.
      
      \emph{Usage:}
      
      \code{   rbind2(flowSet, flowSet)}
      
      \code{   rbind2(flowSet, flowFrame)}
      
    }
    
    \item{spillover}{Compute spillover matrix from a compensation
      set. See separate documentation for details.
    }
  }
}



\section{Important note on storage and performance}{
  
  The bulk of the data in a \code{flowSet} object is stored in an
  \code{\link[base:environment]{environment}}, and is therefore not
  automatically copied when the \code{flowSet} object is copied. If
  \code{x} is an object of class \code{flowSet}, then the code
  \preformatted{y <- x} will create an object \code{y} that contains
  copies of the \code{phenoData} and administrative data in \code{x},
  but refers to the \emph{same} environment with the actual fluorescence
  data. See below for how to create proper copies.

  The reason for this is performance. The pass-by-value semantics of
  function calls in \code{R} can result in numerous copies of the same
  data object being made in the course of a series of nested function
  calls. If the data object is large, this can result in considerable
  cost of memory and performance. \code{flowSet} objects are intended to
  contain experimental data in the order of hundreds of Megabytes, which
  can effectively be treated as read-only: typical tasks are the
  extraction of subsets and the calculation of summary statistics.  This
  is afforded by the design of the \code{flowSet} class: an object of
  that class contains a \code{phenoData} slot, some administrative
  information, and a \emph{reference} to an environment with the
  fluorescence data; when it is copied, only the reference is copied,
  but not the potentially large set of fluorescence data themselves.

  However, note that subsetting operations, such as \code{y <- x[i]} do
  create proper copies, including a copy of the appropriate part of the
  fluorescence data, as it should be expected. Thus, to make a proper
  copy of a \code{flowSet} \code{x}, use \code{y <- x[seq(along=x)]}

}


\author{

  F. Hahne, B. Ellis, P. Haaland and N. Le Meur

}

\seealso{
  
  \code{\linkS4class{flowFrame}}, \code{\link{read.flowSet}}
  
}

\examples{

## load example data and object creation
data(GvHD)

## subsetting to flowSet
set <- GvHD[1:4]
GvHD[1:4,1:2]
sel <- sampleNames(GvHD)[1:2]
GvHD[sel, "FSC-H"]
GvHD[sampleNames(GvHD) == sel[1], colnames(GvHD[1]) == "SSC-H"]

## subsetting to flowFrame
GvHD[[1]]
GvHD[[1, 1:3]]
GvHD[[1, "FSC-H"]]
GvHD[[1, colnames(GvHD[1]) == "SSC-H"]]
GvHD$s5a02

## constructor
flowSet(GvHD[[1]], GvHD[[2]])
pd <- phenoData(GvHD)[1:2,]
flowSet(s5a01=GvHD[[1]], s5a02=GvHD[[2]],phenoData=pd)

## colnames
colnames(set)
colnames(set) <- make.names(colnames(set))

## object name
identifier(set)
identifier(set) <- "test"

## phenoData
pd <- phenoData(set)
pd
pd$test <- "test"
phenoData(set) <- pd
pData(set)
varLabels(set)
varLabels(set)[6] <- "Foo"
varLabels(set)

## sampleNames
sampleNames(set)
sampleNames(set) <- LETTERS[1:length(set)]
sampleNames(set)

## keywords
keyword(set, list("transformation"))

## length
length(set)

## compensation
samp <- read.flowSet(path=system.file("extdata","compdata","data",
package="flowCore"))
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
comp.mat
summary(samp[[1]])
samp <- compensate(samp, as.matrix(comp.mat))
summary(samp[[1]])

## transformation
opar <- par(mfcol=c(1:2))
plot(set[[1]], c("FL1.H", "FL2.H"))
set <- transform(set, FL1.H=log(FL1.H), FL2.H=log(FL2.H))
plot(set[[1]], c("FL1.H", "FL2.H"))
par(opar)

## filtering of flowSets
rectGate <- rectangleGate(filterId="nonDebris", FSC.H=c(200,Inf))
fres <- filter(set, rectGate)
class(fres)
summary(fres[[1]])
rectGate2 <- rectangleGate(filterId="nonDebris2", SSC.H=c(300,Inf))
fres2 <- filter(set, list(A=rectGate, B=rectGate2, C=rectGate, D=rectGate2))

## Splitting frames of a flowSet
split(set, rectGate)
split(set[1:2], rectGate, populatiuon="nonDebris2+")
split(set, c(1,1,2,2))

## subsetting according to filters and filter results
Subset(set, rectGate)
Subset(set, filter(set, rectGate))
Subset(set, list(A=rectGate, B=rectGate2, C=rectGate, D=rectGate2))

## combining flowSets
rbind2(set[1:2], set[3:4])
rbind2(set[1:3], set[[4]])
rbind2(set[[4]], set[1:2])

}

\keyword{classes}
