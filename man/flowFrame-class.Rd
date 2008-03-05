\name{flowFrame-class}
\docType{class}
\alias{flowFrame-class}
\alias{flowFrame}
\alias{[,flowFrame,ANY-method}
\alias{[,flowFrame,filter-method}
\alias{[,flowFrame,filterResult-method}
\alias{==,flowFrame,filterResult-method}
\alias{==,flowFrame,flowFrame-method}
\alias{Subset,flowFrame,ANY-method}
\alias{Subset,flowFrame,logical-method}
\alias{description,flowFrame-method}
\alias{description<-,flowFrame,character-method}
\alias{exprs,flowFrame-method}
\alias{exprs<-,flowFrame,matrix-method}
\alias{colnames,flowFrame-method}
\alias{colnames<-,flowFrame-method}
\alias{names,flowFrame-method}
\alias{parameters,flowFrame-method}
\alias{keyword,flowFrame,function-method}
\alias{keyword,flowFrame,character-method}
\alias{keyword,flowFrame,list-method}
\alias{show,flowFrame-method}
\alias{plot,flowFrame,missing-method}
\alias{plot,flowFrame,character-method}
\alias{plot,flowFrame,filter-method}
\alias{$.flowFrame}
\alias{filter,flowFrame}
\alias{compensate,flowFrame,matrix-method}
\alias{spillover,flowFrame-method}
\alias{transform,flowFrame-method}
\alias{split,flowFrame,ANY,ANY-method}
\alias{split,flowFrame,ANY-method}
\alias{split,flowFrame,filter,ANY-method}
\alias{split,flowFrame,filterSet,ANY-method}
\alias{split,flowFrame,logicalFilterResult,ANY-method}
\alias{split,flowFrame,multipleFilterResult,ANY-method}
\alias{split,flowFrame,manyFilterResult,ANY-method}
\alias{ncol,flowFrame-method}
\alias{nrow,flowFrame-method}
\alias{dim,flowFrame-method}
\alias{identifier,flowFrame-method}
\alias{summary,flowFrame-method}
\alias{featureNames,flowFrame-method}
\alias{each_row,flowFrame-method}
\alias{each_col,flowFrame-method}
\alias{range,flowFrame-method}

  
\title{'flowFrame': a class for storing observed quantitative
  properties for a population of cells from a FACS run
}

\description{This class represents the data contained in a \acronym{FCS}
  file or similar data structure. There are three parts of the data:
  \enumerate{

    \item a numeric matrix of the raw measurement values with
    \kbd{rows=events} and \kbd{columns=parameters}
    
    \item annotation for the parameters (e.g., the measurement channels,
    stains, dynamic range)
    
    \item additional annotation provided through keywords in the
    \acronym{FCS}file
    
  }
}


\section{Creating Objects}{
  Objects can be created using\cr
  \code{  new("flowFrame",}\cr
  \code{    exprs       = ....,  Object of class matrix}\cr
  \code{    parameters  = ....,  Object of class AnnotatedDataFrame}\cr
  \code{    description = ....,  Object of class list}\cr
  \code{  )}\cr
  
  or the function \code{\link[flowCore]{read.FCS}}, which will read the
  content of a \acronym{FCS} file. This is the recommended and safest
  way of object creation, since \code{read.FCS} will perform basic data
  quality checks upon import. Unless you know exactly what you are
  doing, creating objects using \code{new("flowFrame")} is discouraged.

}

\section{Slots}{
  \describe{
    
    \item{\code{exprs}:}{Object of class \code{matrix} containing the
      measured intensities. Rows correspond to cells, columns to the
      different measurement channels. The \code{colnames} attribute of
      the matrix is supposed to hold the names or identifiers for the
      channels. The \code{rownames} attribute would usually not be set.
      }
    
    \item{\code{parameters}:}{An
      \code{\link[Biobase]{AnnotatedDataFrame-class}} containing information
      about each column of the flowFrame. This will generally be filled
      in by \code{read.FCS} or similar functions using data from the
      \acronym{FCS} keywords describing the parameters.}
    
    \item{\code{description}:}{A list containing the meta data included
      in the FCS file.}
    
  }
}

\section{Methods}{
  \describe{
    
    \item{[}{Subsetting. Returns an object of class \code{flowFrame}.
      The subsetting is applied to the \code{exprs} slot, while the
      \code{description} slot is unchanged. The syntax for subsetting is
      similar to that of matrices. In addition to index vectors
      \code{flowFrames} can be subset via \code{\link{filterResult}} and
      \code{\link{filter}} objects.

      \emph{Usage:}
      
      \code{   flowFrame[x,y]}
      
      \code{   flowFrame[filter,]}

      \code{   flowFrame[filterResult,]}

    }
    
    \item{\$}{Subsetting by channel name. This is similar to subsetting
      of columns of data frames, i.e., \code{frame$FSC.H} is equivalent
      to \code{frame[, "FSC.H"]}. Note that column names may have to be
      quoted if they are no valid R symbols
      (e.g. \code{frame$"FSC-H"}).}
  
    \item{exprs, exprs<-}{Extract or replace the raw data
      intensities. The replacement value must be a numeric matrix with
      colnames matching the parameter definitions. Implicit subsetting
      is allowed (i.e. less columns in the replacement value compared to
      the original flowFrame, but all have to be defined there).

      \emph{Usage:}
      
      \code{   exprs(flowFrame)}
      
      \code{   exprs(flowFrame) <- value}

    }
    
    \item{description}{Extract the whole list of annotation
      keywords. Usually one would only be interested in a subset of
      keywords, in which case the \code{keyword} method is more
      appropriate
      
      \emph{Usage:}
      
      \code{   description(flowFrame)}
      
    }

    \item{keyword}{Extract one or more entries from the
      \code{description} slot by keyword. Methods are defined for
      character vectors (select a keyword by name), functions (select a
      keyword by evaluating a function on their content) and for lists
      (a combination of the above).
      
      \emph{Usage:}

      \code{   keyword(flowFrame)}

      \code{   keyword(flowFrame, character)}
      
    }
    
    \item{parameters}{Extract parameters and return an object of class
      \code{AnnotatedDataFrame}. To access the actual parameter
      annotation, use \code{pData(parameters(frame))}.
      
      \emph{Usage:}
      
      \code{   parameters(flowFrame)}

    }

    \item{show}{Display details about the \code{flowFrame} object.}

    \item{summary}{Return descriptive statistical summary (min, max,
      mean and quantile) for each channel

      \emph{Usage:}
      
      \code{   summary(flowFrame)}
      
    }
    
    \item{plot}{Basic plots for \code{flowFrame} objects. If the object
      has only a single parameter this produces a
      \code{\link[graphics::hist]{histogram}}. For exactly two
      parameters we plot a bivariate density map (see
      \code{\link[geneplotter]{smootScatter}}) and for more than two
      parameters we produce a simple \code{\link[graphics]{pairs}}
      plot. To select specific parameters from a \code{flowFrame} for
      plotting, either subset the object or specifiy the parameters as a
      character vector in the second argument to \code{plot}.

      \emph{Usage:}
      
      \code{   plot(flowFrame, ...)}
      
      \code{   plot(flowFrame, character, ...)}
      
    }

    \item{ncol,nrow,dim}{Extract the dimensions of the data matrix.

      \emph{Usage:}
      
      \code{   ncol(flowFrame)}
      
      \code{   nrow(flowFrame)}

      \code{   dim(flowFrame)}

    }

    \item{featureNames, colnames}{Extract parameter names (i.e., the
      colnames of the data matrix). For \code{colnames} there is also a
      replacement method.

      \emph{Usage:}
      
      \code{   featureNames(flowFrame)}
      
      \code{   colnames(flowFrame)}

      \code{   colnames(flowFrame) <- value}
      
    }
    
    \item{names}{Extract pretty formated names of the parameters
      including parameter descriptions.
      
      \emph{Usage:}
      
      \code{   names(flowFrame)}

    }

    \item{identifier}{Extract GUID. Returns file name if not available. See
      \code{\link[flowCore]{identifier}} for details.

      \emph{Usage:}
      
      \code{   identifier(flowFrame)}
    }
    
    \item{range}{Get dynamic range of the flowFame. Note that this is
      not necessarily the range of the actual data values, but the
      theoretical range of values the measurement instrument was able to
      capture. The values of the dynamic range will be transformed when
      using the out-of-line transformation method for
      \code{flowFrames}. Additional character arguments are evaluated as
      parameter names for which to return the dynamic range.
      
      \emph{Usage:}

      \code{   range(flowFrame, ...)}

    }

    \item{each\_row, each\_col}{Apply functions over rows or columns of the
      data matrix. These are convenince methods. To have more control
      over the iterator use the \code{\link{fsApply}} function.

      \emph{Usage:}
      
      \code{   each_row(flowFrame, function, ...)}

      \code{   each_col(flowFrame, function, ...)}
    }
    
    \item{transform}{Apply a transformation function on a
      \code{flowFrame} object. This uses R's
      \code{\link[base]{transfom}} function by treating the
      \code{flowFrame} like a regular \code{data.frame}. \code{flowCore}
      provides an additional inline mechanism for transformations (see
      \code{\link{filter-on-methods}}) which is strictly more limited
      than the out-of-line transformation described here. Note that only
      the out-of-line transformation makes sure that the values for the
      dynamic range in the \code{parameters} slot are also transformed
      accordingly.
      
      \emph{Usage:}
      
      \code{   transform(flowFrame, ...)}

    }
    
    \item{filter}{Apply a \code{filter} object on a \code{flowSet}
      object. This returns an object of class
      \code{\link{filterResult}}, which could then be used for
      subsetting of the data or to calculate summary statistics. See
      \code{\link{filter}} for details.

      \emph{Usage:}
      
      \code{   filter(flowFrame, filter)}

    }
    
    \item{split}{Split \code{flowFrame} object according to a
      \code{\link{filter}}, a \code{\link{filterResult}} or a
      \code{factor}. For most types of filters, an optional
      \code{flowSet=TRUE} parameter will create a \code{flowSet} rather
      than a simple list.

      \emph{Usage:}
      
      \code{   split(flowFrame, filter, flowSet=FALSE, ...)}

      \code{   split(flowFrame, filterResult, flowSet=FALSE, ...)}

      \code{   split(flowFrame, factor, flowSet=FALSE, ...)}

    }
    
    \item{Subset}{Subset a flowFrame according to a \code{filter} or a
      logical vector. The same can be done using the standard subsetting
      operator with a \code{filter}, \code{filterResult}, or a logical
      vector as first argument.

      \emph{Usage:}
      
      \code{   Subset(flowFrame, filter)}

      \code{   Subset(flowFrame, logical)}
      
    }

    \item{compensate}{Apply a compensation matrix on a \code{flowFrame}
      object.
      
      \emph{Usage:}
      
      \code{   compensate(flowFrame, matrix)}

    }

    \item{spillover}{Extract spillover matrix from description slot if
      present.
      
      \emph{Usage:}
      
      \code{   spillover(flowFrame, matrix)}

    }
    
    \item{==}{Test equality between \code{flowFrames}}
  }
}

\details{
  Objects of class \code{flowFrame} can be used to hold arbitrary
  data of cell populations, acquired in flow-cytometry.
  
  \acronym{FCS} is the Data File Standard for Flow Cytometry, the
  current version is FCS 3.0. See the vignette of this package for
  additional information on using the object system for handling of
  flow-cytometry data.

}

\author{F. Hahne, B. Ellis, P. Haaland and N. Le Meur}

\examples{
## load example data
data(GvHD)
frame <- GvHD[[1]]

## ubsetting
frame[1:4,]
frame[,3]
frame[,"FSC-H"]
frame$"SSC-H"

## accessing and replacing raw values
head(exprs(frame))
exprs(frame) <- exprs(frame)[1:3000,]
frame
exprs(frame) <- exprs(frame)[,1:6]
frame

## access FCS keywords
head(description(frame))
keyword(frame, c("FILENAME", "$FIL"))

## parameter annotation
parameters(frame)
pData(parameters(frame))

## summarize frame data
summary(frame)

## plotting
plot(frame)
plot(frame, c("FSC-H", "SSC-H"))
plot(frame[,1])

## frame dimensions
ncol(frame)
nrow(frame)
dim(frame)

## accessing and replacing parameter names
featureNames(frame)
all(featureNames(frame) == colnames(frame))
colnames(frame) <- make.names(colnames(frame))
colnames(frame)
parameters(frame)$name
names(frame)

## accessing a GUID
identifier(frame)

## dynamic range of a frame
range(frame)
range(frame, "FSC.H", "FL1.H")
range(frame)$FSC.H

## iterators
head(each_row(frame, mean))
head(each_col(frame, mean))

## transformation
opar <- par(mfcol=c(1:2))
plot(frame, c("FL1.H", "FL2.H"))
frame <- transform(frame, FL1.H=log(`FL1.H`), FL2.H=log(`FL2.H`))
plot(frame, c("FL1.H", "FL2.H"))
par(opar)
range(frame)

## filtering of flowFrames
rectGate <- rectangleGate(filterId="nonDebris","FSC.H"=c(200,Inf))
fres <- filter(frame, rectGate)
summary(fres)

## splitting of flowFrames
split(frame, rectGate)
split(frame, rectGate, flowSet=TRUE)
split(frame, fres)
f <- cut(exprs(frame$FSC.H), 3)
split(frame, f)

## subsetting according to filters and filter results
Subset(frame, rectGate)
Subset(frame, fres)
Subset(frame, as.logical(exprs(frame$FSC.H) < 300))
frame[rectGate,]
frame[fres,]

## accessing the spillover matrix
try(spillover(frame))

## check equality
frame2 <- frame
frame == frame2
exprs(frame2) <- exprs(frame)*2
frame == frame2

}

\keyword{classes}

