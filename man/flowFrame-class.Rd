\name{flowFrame-class}
\docType{class}
\alias{flowFrame-class}
\alias{flowFrame}
\alias{[,flowFrame,ANY-method}
\alias{[,flowFrame,filter-method}
\alias{[,flowFrame,filterResult-method}
\alias{$.flowFrame}
\alias{exprs}
\alias{exprs<-}
\alias{exprs,flowFrame-method}
\alias{exprs<-,flowFrame,matrix-method}
\alias{exprs<-,flowFrame,ANY-method}
\alias{initialize,flowFrame-method}
\alias{head,flowFrame-method}
\alias{tail,flowFrame-method}
\alias{description}
\alias{description,flowFrame-method}
\alias{description<-,flowFrame,list-method}
\alias{description<-,flowFrame,ANY-method}
\alias{show,flowFrame-method}
\alias{plot,flowFrame,ANY-method}
\alias{plot,flowFrame-method}
\alias{summary,flowFrame-method}
\alias{ncol}
\alias{ncol,flowFrame-method}
\alias{nrow}
\alias{nrow,flowFrame-method}
\alias{dim}
\alias{dim,flowFrame-method}
\alias{featureNames}
\alias{featureNames,flowFrame-method}
\alias{colnames}
\alias{colnames,flowFrame-method}
\alias{colnames<-}
\alias{colnames<-,flowFrame-method}
\alias{names}
\alias{names,flowFrame-method}
\alias{range}
\alias{range,flowFrame-method}
\alias{cbind2,flowFrame,matrix-method}
\alias{cbind2,flowFrame,numeric-method}
\alias{transform,flowFrame-method}
\alias{compensate,flowFrame,matrix-method}
\alias{compensate,flowFrame,data.frame-method}
\alias{compensate,flowFrame,compensation-method}
\alias{==,flowFrame,filterResult-method}
\alias{==,flowFrame,flowFrame-method}
\alias{<,flowFrame,ANY-method}
\alias{<=,flowFrame,ANY-method}
\alias{>,flowFrame,ANY-method}
\alias{>=,flowFrame,ANY-method}
\alias{spillover,flowFrame-method}


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
    \acronym{FCS} file
    
  }
}


\section{Creating Objects}{
  Objects can be created using\cr
  \code{  new("flowFrame",}\cr
  \code{    exprs       = ....,  Object of class matrix}\cr
  \code{    parameters  = ....,  Object of class AnnotatedDataFrame}\cr
  \code{    description = ....,  Object of class list}\cr
  \code{  )}\cr

  or the constructor \code{flowFrame}, with mandatory arguments
  \code{exprs} and optional arguments \code{parameters} and
  \code{description}.

  \code{flowFrame(exprs, parameters, description=list())}
  
  To create a \code{flowFrame} directly from an \acronym{FCS} file, use
  function \code{\link[flowCore]{read.FCS}}. This is the recommended and
  safest way of object creation, since \code{read.FCS} will perform
  basic data quality checks upon import. Unless you know exactly what
  you are doing, creating objects using \code{new} or the constructor is
  discouraged.

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
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
      containing information about each column of the
      \code{flowFrame}. This will generally be filled in by
      \code{read.FCS} or similar functions using data from the
      \acronym{FCS} keywords describing the parameters.}
    
    \item{\code{description}:}{A list containing the meta data included
      in the FCS file.}
    
  }
}

\section{Methods}{
  \describe{

    There are separate documentation pages for most of the methods
    listed here which should be consulted for more details.
    
    \item{[}{Subsetting. Returns an object of class \code{flowFrame}.
      The subsetting is applied to the \code{exprs} slot, while the
      \code{description} slot is unchanged. The syntax for subsetting is
      similar to that of \code{\link[=data.frame]{data.frames}}. In
      addition to the usual index vectors (integer and logical by
      position, character by parameter names), \code{flowFrames} can be
      subset via \code{\link{filterResult}} and
      \code{\linkS4class{filter}} objects.

      \emph{Usage:}
      
      \code{   flowFrame[i,j]}
      
      \code{   flowFrame[filter,]}

      \code{   flowFrame[filterResult,]}

      Note that the value of argument \code{drop} is ignored when
      subsetting \code{flowFrames}.
      
    }
    
    \item{\$}{Subsetting by channel name. This is similar to subsetting
      of columns of \code{\link[=data.frame]{data.frames}}, i.e.,
      \code{frame$FSC.H} is equivalent to \code{frame[, "FSC.H"]}. Note
      that column names may have to be quoted if they are no valid R
      symbols (e.g. \code{frame$"FSC-H"}).

    }
  
    \item{exprs, exprs<-}{Extract or replace the raw data
      intensities. The replacement value must be a numeric matrix with
      colnames matching the parameter definitions. Implicit subsetting
      is allowed (i.e. less columns in the replacement value compared to
      the original \code{flowFrame}, but all have to be defined there).

      \emph{Usage:}
      
      \code{   exprs(flowFrame)}
      
      \code{   exprs(flowFrame) <- value}

    }

    \item{head, tail}{Show first/last elements of the raw data matrix

      \emph{Usage:}
      
      \code{   head(flowFrame)}
      
      \code{   tail(flowFrame)}
      
    }
    
    \item{description, description<-}{Extract or replace the whole list
      of annotation keywords. Usually one would only be interested in a
      subset of keywords, in which case the \code{keyword} method is
      more appropriate. The optional \code{hideInternal} parameter can
      be used to exclude internal FCS parameters starting
      with\code{\\$}.
      
      \emph{Usage:}
      
      \code{   description(flowFrame)}

      \code{   description(flowFrame) <- value}
      
    }

    \item{keyword, keyword<-}{Extract ore replace one or more entries
      from the \code{description} slot by keyword. Methods are defined
      for character vectors (select a keyword by name), functions
      (select a keyword by evaluating a function on their content) and
      for lists (a combination of the above). See \code{\link{keyword}}
      for details.
      
      \emph{Usage:}

      \code{   keyword(flowFrame)}

      \code{   keyword(flowFrame, character)}

      \code{   keyword(flowFrame, list)}
      
      \code{   keyword(flowFrame) <- list(value) }
      
    }
    
    \item{parameters, parameters<-}{Extract parameters and return an
      object of class
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}},
      or replace such an object. To access the actual parameter
      annotation, use \code{pData(parameters(frame))}. Replacement is
      only valid with
      \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrames}}
      containing all varLabels \code{name}, \code{desc}, \code{range},
      \code{minRange} and \code{maxRange}, and matching entries in the
      \code{name} column to the colnames of the \code{exprs} matrix. See
      \code{\link{parameters}} for more details.
      
      \emph{Usage:}
      
      \code{   parameters(flowFrame)}

      \code{   parameters(flowFrame) <- value}

    }

    \item{show}{

      Display details about the \code{flowFrame} object.
      
    }

    \item{summary}{Return descriptive statistical summary (min, max,
      mean and quantile) for each channel

      \emph{Usage:}
      
      \code{   summary(flowFrame)}
      
    }
    
    \item{plot}{Basic plots for \code{flowFrame} objects. If the object
      has only a single parameter this produces a
      \code{\link[graphics:hist]{histogram}}. For exactly two parameters
      we plot a bivariate density map (see
      \code{\link[graphics]{smoothScatter}}
      and for more than two parameters we produce a simple
      \code{\link[lattice]{splom}} plot. To select specific parameters
      from a \code{flowFrame} for plotting, either subset the object or
      specify the parameters as a character vector in the second
      argument to \code{plot}. The smooth parameters lets you toggle
      between density-type
      \code{\link[graphics]{smoothScatter}}
      plots and regular scatterplots.  For far more sophisticated
      plotting of flow cytometry data, see the
      \code{\link[flowViz:flowViz-package]{flowViz}} package.

      \emph{Usage:}
      
      \code{   plot(flowFrame, ...)}
      
      \code{   plot(flowFrame, character, ...)}

      \code{   plot(flowFrame, smooth=FALSE, ...)}
      
    }

    \item{ncol, nrow, dim}{Extract the dimensions of the data matrix.

      \emph{Usage:}
      
      \code{   ncol(flowFrame)}
      
      \code{   nrow(flowFrame)}

      \code{   dim(flowFrame)}

    }

    \item{featureNames, colnames, colnames<-}{Extract parameter names (i.e., the
      colnames of the data matrix). \code{colnames} and
      \code{featureNames} are synonyms. For \code{colnames} there is
      also a replacement method. This will update the \code{name} column
      in the \code{parameters} slot as well.

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

    \item{identifier}{Extract GUID of a \code{flowFrame}. Returns the
      file name if no GUID is available. See \code{\link{identifier}}
      for details.

      \emph{Usage:}
      
      \code{   identifier(flowFrame)}
    }
    
    \item{range}{Get dynamic range of the \code{flowFame}. Note that
      this is not necessarily the range of the actual data values, but
      the theoretical range of values the measurement instrument was
      able to capture. The values of the dynamic range will be
      transformed when using the transformation methods for
      \code{flowFrames}. Additional character arguments are evaluated as
      parameter names for which to return the dynamic range.
      
      \emph{Usage:}

      \code{   range(flowFrame, ...)}

    }

    \item{each\_row, each\_col}{Apply functions over rows or columns of
      the data matrix. These are convenience methods. See
      \code{\link{each_col}} for details.

      \emph{Usage:}
      
      \code{   each_row(flowFrame, function, ...)}

      \code{   each_col(flowFrame, function, ...)}
    }
    
    \item{transform}{Apply a transformation function on a
      \code{flowFrame} object. This uses R's
      \code{\link[base]{transform}} function by treating the
      \code{flowFrame} like a regular \code{data.frame}. \code{flowCore}
      provides an additional inline mechanism for transformations (see
      \code{\link{\%on\%}}) which is strictly more limited
      than the out-of-line transformation described here.
      
      \emph{Usage:}
      
      \code{   transform(flowFrame, ...)}

    }
    
    \item{filter}{Apply a \code{\linkS4class{filter}} object on a
      \code{flowFrame} object. This returns an object of class
      \code{\link{filterResult}}, which could then be used for
      subsetting of the data or to calculate summary statistics. See
      \code{\link{filter}} for details.

      \emph{Usage:}
      
      \code{   filter(flowFrame, filter)}

    }
    
    \item{split}{Split \code{flowFrame} object according to a
      \code{\link{filter}}, a \code{\link{filterResult}} or a
      \code{factor}. For most types of filters, an optional
      \code{flowSet=TRUE} parameter will create a
      \code{\linkS4class{flowSet}} rather than a simple list. See
      \code{\link{split}} for details.

      \emph{Usage:}
      
      \code{   split(flowFrame, filter, flowSet=FALSE, ...)}

      \code{   split(flowFrame, filterResult, flowSet=FALSE, ...)}

      \code{   split(flowFrame, factor, flowSet=FALSE, ...)}

    }
    
    \item{Subset}{Subset a \code{flowFrame} according to a \code{filter}
      or a logical vector. The same can be done using the standard
      subsetting operator with a \code{filter}, \code{filterResult}, or
      a logical vector as first argument.

      \emph{Usage:}
      
      \code{   Subset(flowFrame, filter)}

      \code{   Subset(flowFrame, logical)}
      
    }

    \item{cbind2}{Expand a \code{flowFrame} by the data in a
      \code{numeric matrix} of the same length. The \code{matrix} must
      have column names different from those of the
      \code{flowFrame}. The additional method for \code{numerics} only
      raises a useful error message.

      \emph{Usage:}
      
      \code{   cbind2(flowFrame, matrix)}

      \code{   cbind2(flowFrame, numeric)}
      
    }
    

    \item{compensate}{Apply a compensation matrix (or a
      \code{\linkS4class{compensation}} object) on a \code{flowFrame}
      object. This returns a compensated \code{flowFrame}.
      
      \emph{Usage:}
      
      \code{   compensate(flowFrame, matrix)}
      \code{   compensate(flowFrame, data.frame)}

    }

    \item{spillover}{Extract spillover matrix from description slot if
      present.
      
      \emph{Usage:}
      
      \code{   spillover(flowFrame, matrix)}

    }
    
    \item{==}{Test equality between two \code{flowFrames}}

    \item{<, >, <=, >=}{These operators basically treat the
      \code{flowFrame} as a numeric matrix.}
    
    \item{\code{initialize(flowFrame)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
  }
}



\details{
  
  Objects of class \code{flowFrame} can be used to hold arbitrary data
  of cell populations, acquired in flow-cytometry.
  
  \acronym{FCS} is the Data File Standard for Flow Cytometry, the
  current version is FCS 3.0. See the vignette of this package for
  additional information on using the object system for handling of
  flow-cytometry data.

}

\author{

  F. Hahne, B. Ellis, P. Haaland and N. Le Meur

}

\seealso{
  
  \code{\linkS4class{flowSet}},  \code{\link{read.FCS}}
  
}

\examples{

## load example data
data(GvHD)
frame <- GvHD[[1]]

## subsetting
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
if(require(flowViz)){
plot(frame)
plot(frame, c("FSC-H", "SSC-H"))
plot(frame[,1])
plot(frame, c("FSC-H", "SSC-H"), smooth=FALSE)
}

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
identifier(frame) <- "test"

## dynamic range of a frame
range(frame)
range(frame, "FSC.H", "FL1.H")
range(frame)$FSC.H

## iterators
head(each_row(frame, mean))
head(each_col(frame, mean))

## transformation
opar <- par(mfcol=c(1:2))
if(require(flowViz))
plot(frame, c("FL1.H", "FL2.H"))
frame <- transform(frame, FL1.H=log(`FL1.H`), FL2.H=log(`FL2.H`))
if(require(flowViz))
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

