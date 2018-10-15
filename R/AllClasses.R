## =========================================================================##
## =========================================================================##
##                    Class definitions and contructors                     ##
## =========================================================================##
## =========================================================================##






## ===========================================================================
##  Some helpers
## ---------------------------------------------------------------------------
## Check for the class of object x and its length and cast error if wrong
checkClass <- function(x, class, length=NULL, verbose=FALSE,
                       mandatory=TRUE)
{
    if(mandatory && missing(x))
        stop("Argument '", substitute(x), "' missing with no default",
             call.=verbose)
    msg <- paste("'", substitute(x), "' must be object of class ",
                 paste("'", class, "'", sep="", collapse=" or "), sep="")
    fail <- !any(sapply(class, function(c, y) is(y, c), x))
    if(!is.null(length) && length(x) != length)
    {
        if(!is.null(x))
        {
            fail <- TRUE
            msg <- paste(msg, "of length", length)
        }
    }
    if(fail) stop(msg, call.=verbose) else invisible(NULL)     
}



## ===========================================================================
##  flowFrame
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, parameters
## and description. exprs contains measurement values, description contains 
## information from file headers of FCS file and parameters contains
## information about the FCS measurement parameters (i.e. channels) available.
## Exprs is a matrix (values are stored in internal memory) 
## ---------------------------------------------------------------------------
#' 'flowFrame': a class for storing observed quantitative properties for a
#' population of cells from a FACS run
#' 
#' This class represents the data contained in a \acronym{FCS} file or similar
#' data structure. There are three parts of the data: \enumerate{
#' \item a numeric matrix of the raw measurement values with \kbd{rows=events}
#' and \kbd{columns=parameters}
#' \item annotation for the parameters (e.g., the measurement channels, stains,
#' dynamic range)
#' \item additional annotation provided through keywords in the \acronym{FCS}
#' file
#' }
#' 
#' 
#' 
#' Objects of class \code{flowFrame} can be used to hold arbitrary data of cell
#' populations, acquired in flow-cytometry.
#' 
#' \acronym{FCS} is the Data File Standard for Flow Cytometry, the current
#' version is FCS 3.0. See the vignette of this package for additional
#' information on using the object system for handling of flow-cytometry data.
#' 
#' @name flowFrame-class
#' @aliases flowFrame-class flowFrame [,flowFrame,ANY-method
#' [,flowFrame,filter-method [,flowFrame,filterResult-method $.flowFrame exprs
#' exprs<- exprs,flowFrame-method exprs<-,flowFrame,matrix-method
#' exprs<-,flowFrame,ANY-method initialize,flowFrame-method
#' head,flowFrame-method tail,flowFrame-method description
#' description,flowFrame-method description<-,flowFrame,list-method
#' description<-,flowFrame,ANY-method show,flowFrame-method
#' plot,flowFrame,ANY-method plot,flowFrame-method summary,flowFrame-method
#' ncol,flowFrame-method nrow,flowFrame-method dim dim,flowFrame-method
#' featureNames featureNames,flowFrame-method colnames,flowFrame-method
#' colnames<- colnames<-,flowFrame-method names names,flowFrame-method range
#' range,flowFrame-method cbind2,flowFrame,matrix-method
#' cbind2,flowFrame,numeric-method transform,flowFrame-method
#' compensate,flowFrame,matrix-method compensate,flowFrame,data.frame-method
#' compensate,flowFrame,compensation-method ==,flowFrame,filterResult-method
#' ==,flowFrame,flowFrame-method <,flowFrame,ANY-method <=,flowFrame,ANY-method
#' >,flowFrame,ANY-method >=,flowFrame,ANY-method spillover,flowFrame-method
#' @docType class
#' 
#' @slot exprs {Object of class \code{matrix} containing the
#' measured intensities. Rows correspond to cells, columns to the
#' different measurement channels. The \code{colnames} attribute of
#' the matrix is supposed to hold the names or identifiers for the
#' channels. The \code{rownames} attribute would usually not be set.
#' }
#' @slot parameters {An
#' \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#' containing information about each column of the
#' \code{flowFrame}. This will generally be filled in by
#' \code{read.FCS} or similar functions using data from the
#' \acronym{FCS} keywords describing the parameters.}
#' @slot description {A list containing the meta data included
#' in the FCS file.}
#' 
#' @section Creating Objects: 
#' Objects can be created using\cr \code{
#' new("flowFrame",}\cr \code{ exprs = ...., Object of class matrix}\cr \code{
#' parameters = ...., Object of class AnnotatedDataFrame}\cr \code{ description
#' = ...., Object of class list}\cr \code{ )}\cr
#' 
#' or the constructor \code{flowFrame}, with mandatory arguments \code{exprs}
#' and optional arguments \code{parameters} and \code{description}.
#' 
#' \code{flowFrame(exprs, parameters, description=list())}
#' 
#' To create a \code{flowFrame} directly from an \acronym{FCS} file, use
#' function \code{\link[flowCore]{read.FCS}}. This is the recommended and
#' safest way of object creation, since \code{read.FCS} will perform basic data
#' quality checks upon import. Unless you know exactly what you are doing,
#' creating objects using \code{new} or the constructor is discouraged. 
#' 
#' @section Methods:
#'   There are separate documentation pages for most of the methods
#'   listed here which should be consulted for more details.
#'   \describe{
#'   \item{[}{Subsetting. Returns an object of class \code{flowFrame}.
#'     The subsetting is applied to the \code{exprs} slot, while the
#'     \code{description} slot is unchanged. The syntax for subsetting is
#'     similar to that of \code{\link[=data.frame]{data.frames}}. In
#'     addition to the usual index vectors (integer and logical by
#'                                          position, character by parameter names), \code{flowFrames} can be
#'     subset via \code{\link{filterResult}} and
#'     \code{\linkS4class{filter}} objects.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   flowFrame[i,j]}
#'     
#'     \code{   flowFrame[filter,]}
#'     
#'     \code{   flowFrame[filterResult,]}
#'     
#'     Note that the value of argument \code{drop} is ignored when
#'     subsetting \code{flowFrames}.
#'     
#'   }
#'   \item{$}{Subsetting by channel name. This is similar to subsetting
#'     of columns of \code{\link[=data.frame]{data.frames}}, i.e.,
#'     \code{frame$FSC.H} is equivalent to \code{frame[, "FSC.H"]}. Note
#'     that column names may have to be quoted if they are no valid R
#'     symbols (e.g. \code{frame$"FSC-H"}).
#'     
#'   }
#'   \item{exprs, exprs<-}{Extract or replace the raw data
#'     intensities. The replacement value must be a numeric matrix with
#'     colnames matching the parameter definitions. Implicit subsetting
#'     is allowed (i.e. less columns in the replacement value compared to
#'                 the original \code{flowFrame}, but all have to be defined there).
#'     
#'     \emph{Usage:}
#'     
#'     \code{   exprs(flowFrame)}
#'     
#'     \code{   exprs(flowFrame) <- value}
#'     
#'   }
#'   \item{head, tail}{Show first/last elements of the raw data matrix
#'     
#'     \emph{Usage:}
#'     
#'     \code{   head(flowFrame)}
#'     
#'     \code{   tail(flowFrame)}
#'     
#'   }
#'   \item{description, description<-}{Extract or replace the whole list
#'     of annotation keywords. Usually one would only be interested in a
#'     subset of keywords, in which case the \code{keyword} method is
#'     more appropriate. The optional \code{hideInternal} parameter can
#'     be used to exclude internal FCS parameters starting
#'     with\code{\\$}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   description(flowFrame)}
#'     
#'     \code{   description(flowFrame) <- value}
#'     
#'   }
#'   \item{keyword, keyword<-}{Extract ore replace one or more entries
#'     from the \code{description} slot by keyword. Methods are defined
#'     for character vectors (select a keyword by name), functions
#'     (select a keyword by evaluating a function on their content) and
#'     for lists (a combination of the above). See \code{\link{keyword}}
#'     for details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   keyword(flowFrame)}
#'     
#'     \code{   keyword(flowFrame, character)}
#'     
#'     \code{   keyword(flowFrame, list)}
#'     
#'     \code{   keyword(flowFrame) <- list(value) }
#'     
#'   }
#'   \item{parameters, parameters<-}{Extract parameters and return an
#'     object of class
#'     \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}},
#'     or replace such an object. To access the actual parameter
#'     annotation, use \code{pData(parameters(frame))}. Replacement is
#'     only valid with
#'     \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrames}}
#'     containing all varLabels \code{name}, \code{desc}, \code{range},
#'     \code{minRange} and \code{maxRange}, and matching entries in the
#'     \code{name} column to the colnames of the \code{exprs} matrix. See
#'     \code{\link{parameters}} for more details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   parameters(flowFrame)}
#'     
#'     \code{   parameters(flowFrame) <- value}
#'     
#'   }
#'   \item{show}{
#'     
#'     Display details about the \code{flowFrame} object.
#'     
#'   }
#'   \item{summary}{Return descriptive statistical summary (min, max,
#'                                                          mean and quantile) for each channel
#'     
#'     \emph{Usage:}
#'     
#'     \code{   summary(flowFrame)}
#'     
#'   }
#'   \item{plot}{Basic plots for \code{flowFrame} objects. If the object
#'     has only a single parameter this produces a
#'     \code{\link[graphics:hist]{histogram}}. For exactly two parameters
#'     we plot a bivariate density map (see
#'                                      \code{\link[graphics]{smoothScatter}}
#'                                      and for more than two parameters we produce a simple
#'                                      \code{\link[lattice]{splom}} plot. To select specific parameters
#'                                      from a \code{flowFrame} for plotting, either subset the object or
#'                                      specify the parameters as a character vector in the second
#'                                      argument to \code{plot}. The smooth parameters lets you toggle
#'                                      between density-type
#'                                      \code{\link[graphics]{smoothScatter}}
#'                                      plots and regular scatterplots.  For far more sophisticated
#'                                      plotting of flow cytometry data, see the
#'                                      \code{\link[flowViz:flowViz-package]{flowViz}} package.
#'                                      
#'                                      \emph{Usage:}
#'                                      
#'                                      \code{   plot(flowFrame, ...)}
#'                                      
#'                                      \code{   plot(flowFrame, character, ...)}
#'                                      
#'                                      \code{   plot(flowFrame, smooth=FALSE, ...)}
#'                                      
#'   }
#'   \item{ncol, nrow, dim}{Extract the dimensions of the data matrix.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   ncol(flowFrame)}
#'     
#'     \code{   nrow(flowFrame)}
#'     
#'     \code{   dim(flowFrame)}
#'     
#'   }
#'   \item{featureNames, colnames, colnames<-}{. \code{colnames} and
#'     \code{featureNames} are synonyms, they extract parameter names (i.e., the
#'                                                                     colnames of the data matrix) .
#'     For \code{colnames} there is
#'     also a replacement method. This will update the \code{name} column
#'     in the \code{parameters} slot as well.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   featureNames(flowFrame)}
#'     
#'     \code{   colnames(flowFrame)}
#'     
#'     \code{   colnames(flowFrame) <- value}
#'     
#'   }
#'   \item{names}{Extract pretty formated names of the parameters
#'     including parameter descriptions.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   names(flowFrame)}
#'     
#'   }
#'   \item{identifier}{Extract GUID of a \code{flowFrame}. Returns the
#'     file name if no GUID is available. See \code{\link{identifier}}
#'     for details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   identifier(flowFrame)}
#'   }
#'   \item{range}{Get instrument or actual data range of the \code{flowFame}. Note that
#'     instrument dynamic range is not necessarily the same as the range of the actual data values, but
#'     the theoretical range of values the measurement instrument was
#'     able to capture. The values of the dynamic range will be
#'     transformed when using the transformation methods for\code{flowFrames}.
#'     
#'     parameters:
#'       
#'       x: flowFrame object.
#'     
#'     type: Range type. either "instrument" or "data". Default is "instrument"
#'     
#'     \emph{Usage:}
#'     
#'     \code{   range(x, type = "data")}
#'     
#'   }
#'   \item{each_row, each_col}{Apply functions over rows or columns of
#'     the data matrix. These are convenience methods. See
#'     \code{\link{each_col}} for details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   each_row(flowFrame, function, ...)}
#'     
#'     \code{   each_col(flowFrame, function, ...)}
#'   }
#'   \item{transform}{Apply a transformation function on a
#'     \code{flowFrame} object. This uses R's
#'     \code{\link[base]{transform}} function by treating the
#'     \code{flowFrame} like a regular \code{data.frame}. \code{flowCore}
#'     provides an additional inline mechanism for transformations (see
#'     \code{\link{\%on\%}}) which is strictly more limited
#'     than the out-of-line transformation described here.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   transform(flowFrame, translist, ...)}
#'     
#'   }
#'   \item{filter}{Apply a \code{\linkS4class{filter}} object on a
#'     \code{flowFrame} object. This returns an object of class
#'     \code{\link{filterResult}}, which could then be used for
#'     subsetting of the data or to calculate summary statistics. See
#'     \code{\link{filter}} for details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   filter(flowFrame, filter)}
#'     
#'     }
#'   \item{split}{Split \code{flowFrame} object according to a
#'     \code{\link{filter}}, a \code{\link{filterResult}} or a
#'     \code{factor}. For most types of filters, an optional
#'     \code{flowSet=TRUE} parameter will create a
#'     \code{\linkS4class{flowSet}} rather than a simple list. See
#'     \code{\link{split}} for details.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   split(flowFrame, filter, flowSet=FALSE, ...)}
#'     
#'     \code{   split(flowFrame, filterResult, flowSet=FALSE, ...)}
#'     
#'     \code{   split(flowFrame, factor, flowSet=FALSE, ...)}
#'     
#'     }
#'   \item{Subset}{Subset a \code{flowFrame} according to a \code{filter}
#'     or a logical vector. The same can be done using the standard
#'     subsetting operator with a \code{filter}, \code{filterResult}, or
#'     a logical vector as first argument.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   Subset(flowFrame, filter)}
#'     
#'     \code{   Subset(flowFrame, logical)}
#'     
#'     }
#'   \item{cbind2}{Expand a \code{flowFrame} by the data in a
#'     \code{numeric matrix} of the same length. The \code{matrix} must
#'     have column names different from those of the
#'     \code{flowFrame}. The additional method for \code{numerics} only
#'     raises a useful error message.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   cbind2(flowFrame, matrix)}
#'     
#'     \code{   cbind2(flowFrame, numeric)}
#'      
#'     }
#'   \item{compensate}{Apply a compensation matrix (or a
#'     \code{\linkS4class{compensation}} object) on a \code{flowFrame}
#'     object. This returns a compensated \code{flowFrame}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   compensate(flowFrame, matrix)}
#'     \code{   compensate(flowFrame, data.frame)}
#'     
#'     }
#'   \item{spillover}{Extract spillover matrix from description slot if
#'     present. It is equivalent to 
#'     \code{keyword(x, c("spillover", "SPILL"))}
#'     Thus will simply return a list of keywords value for "spillover" and "SPILL".
#'     
#'     \emph{Usage:}
#'     
#'     \code{   spillover(flowFrame)}
#'     
#'     }
#'   \item{==}{Test equality between two \code{flowFrames}}
#'   \item{<, >, <=, >=}{These operators basically treat the
#'     \code{flowFrame} as a numeric matrix.}
#'   \item{\code{initialize(flowFrame)}:}{Object instantiation, used
#'     by \code{new}; not to be called directly by the user.}
#' }
#' 
#' @author
#' 
#' F. Hahne, B. Ellis, P. Haaland and N. Le Meur
#' @seealso
#' 
#' \code{\linkS4class{flowSet}}, \code{\link{read.FCS}}
#' @keywords classes
#' @examples
#' 
#' ## load example data
#' data(GvHD)
#' frame <- GvHD[[1]]
#' 
#' ## subsetting
#' frame[1:4,]
#' frame[,3]
#' frame[,"FSC-H"]
#' frame$"SSC-H"
#' 
#' ## accessing and replacing raw values
#' head(exprs(frame))
#' exprs(frame) <- exprs(frame)[1:3000,]
#' frame
#' exprs(frame) <- exprs(frame)[,1:6]
#' frame
#' 
#' ## access FCS keywords
#' head(description(frame))
#' keyword(frame, c("FILENAME", "$FIL"))
#' 
#' ## parameter annotation
#' parameters(frame)
#' pData(parameters(frame))
#' 
#' ## summarize frame data
#' summary(frame)
#' 
#' ## plotting
#' plot(frame)
#' if(require(flowViz)){
#' plot(frame)
#' plot(frame, c("FSC-H", "SSC-H"))
#' plot(frame[,1])
#' plot(frame, c("FSC-H", "SSC-H"), smooth=FALSE)
#' }
#' 
#' ## frame dimensions
#' ncol(frame)
#' nrow(frame)
#' dim(frame)
#' 
#' ## accessing and replacing parameter names
#' featureNames(frame)
#' all(featureNames(frame) == colnames(frame))
#' colnames(frame) <- make.names(colnames(frame))
#' colnames(frame)
#' parameters(frame)$name
#' names(frame)
#' 
#' ## accessing a GUID
#' identifier(frame)
#' identifier(frame) <- "test"
#' 
#' ##  range of a frame
#' range(frame) #instrument range
#' range(frame, type = "data") #actual data range
#' range(frame)$FSC.H
#' 
#' ## iterators
#' head(each_row(frame, mean))
#' head(each_col(frame, mean))
#' 
#' ## transformation
#' opar <- par(mfcol=c(1:2))
#' if(require(flowViz))
#' plot(frame, c("FL1.H", "FL2.H"))
#' frame <- transform(frame, transformList(c("FL1.H", "FL2.H"), log))
#' if(require(flowViz))
#' plot(frame, c("FL1.H", "FL2.H"))
#' par(opar)
#' range(frame)
#' 
#' ## filtering of flowFrames
#' rectGate <- rectangleGate(filterId="nonDebris","FSC.H"=c(200,Inf))
#' fres <- filter(frame, rectGate)
#' summary(fres)
#' 
#' ## splitting of flowFrames
#' split(frame, rectGate)
#' split(frame, rectGate, flowSet=TRUE)
#' split(frame, fres)
#' f <- cut(exprs(frame$FSC.H), 3)
#' split(frame, f)
#' 
#' ## subsetting according to filters and filter results
#' Subset(frame, rectGate)
#' Subset(frame, fres)
#' Subset(frame, as.logical(exprs(frame$FSC.H) < 300))
#' frame[rectGate,]
#' frame[fres,]
#' 
#' ## accessing the spillover matrix
#' try(spillover(frame))
#' 
#' ## check equality
#' frame2 <- frame
#' frame == frame2
#' exprs(frame2) <- exprs(frame)*2
#' frame == frame2
#' 
#' 
setClass("flowFrame",                
         representation=representation(exprs="matrix",
         parameters="AnnotatedDataFrame",
         description="list"),
         prototype=list(exprs=matrix(numeric(0),
                        nrow=0,
                        ncol=0),
         parameters=new("AnnotatedDataFrame"),
         description=list(note="empty")))

## helper function to create empty AnnotatedDataFrame for the parameters slot
parDefault <- function(exp)
{
    vm <- data.frame(labelDescription=c(name="Name of Parameter",
                     desc="Description of Parameter",
                     range="Range of Parameter",
                     minRange="Minimum Parameter Value after Transformation",
                     maxRange="Maximum Parameter Value after Transformation"))
    cols <- colnames(exp)
    pd <- data.frame(name=cols, desc=cols,
                     range=apply(exp, 2, max, na.rm=TRUE),
                     minRange=apply(exp, 2, min, na.rm=TRUE),
                     maxRange=apply(exp, 2, max, na.rm=TRUE)
                     , row.names = paste0("$P", seq_along(cols)))
    new("AnnotatedDataFrame", pd, vm)
}

## check parameter AnnotatedDataFrame for validity
isValidParameters <- function(parameters, exprs)
{
    checkClass(parameters, "AnnotatedDataFrame")
    if(!all(c("name", "desc", "range", "minRange", "maxRange")
            %in% varLabels(parameters)))
        stop("The following columns are mandatory:\n  'name', 'desc',",
             "'range', 'minRange', 'maxRange'", call.=FALSE)
    if(!missing(exprs))
        if(!all(colnames(exprs) %in% parameters$name))
            stop("parameter description doesn't match colnames of the ",
                 "data matrix", call.=FALSE)
    return(TRUE)
}

## constructor
flowFrame <- function(exprs, parameters, description=list())
{
    if(!is.matrix(exprs) || !is.numeric(exprs) || is.null(colnames(exprs)))
        stop("Argument 'exprs' must be numeric matrix with colnames ",
             "attribute set", call.=FALSE)
    if(missing(parameters))
        parameters <- parDefault(exprs)
    else
        isValidParameters(parameters, exprs)
    checkClass(description, "list")
    new("flowFrame", exprs=exprs, parameters=parameters,
        description=description)
}



## ===========================================================================
##  flowSet
## ---------------------------------------------------------------------------
## A collection of several cytoFrames making up one experiment. Slots 
## frames, phenoData, colnames. Frames contains the cytoFrame objects,
## phenoData the experiment meta data and colnames the channel names.
## An additional character scalar _.name._ is stored in the environment
## which holds a name of the object that will be used in the workFlows.
## By storing it in the environment we don't have to add an additional
## slot and defunct old serialized flowSet objects.
## ---------------------------------------------------------------------------
#' 'flowSet': a class for storing flow cytometry raw data from quantitative
#' cell-based assays
#' 
#' This class is a container for a set of \code{\linkS4class{flowFrame}}
#' objects
#' 
#' 
#' @name flowSet-class
#' @aliases flowSet-class flowSet [,flowSet-method [,flowSet,ANY-method
#' $,flowSet-method [[,flowSet-method [[,flowSet,ANY-method [[<-,flowSet-method
#' [[<-,flowSet,ANY,ANY,flowFrame-method [[<-,flowFrame-method
#' fsApply,flowSet-method show,flowSet-method length,flowSet-method
#' colnames,flowSet-method colnames<-,flowSet-method identifier,flowSet-method
#' identifier<-,flowSet,ANY-method sampleNames,flowSet-method
#' sampleNames<-,flowSet,ANY-method phenoData,flowSet-method
#' phenoData<-,flowSet,ANY-method phenoData<-,flowSet,phenoData-method
#' pData,flowSet-method pData<-,flowSet,data.frame-method
#' plot,flowSet,ANY-method plot,flowSet-method varLabels,flowSet-method
#' varLabels<-,flowSet-method varLabels<-,flowSet,ANY-method
#' varMetadata,flowSet-method varMetadata<-,flowSet,ANY-method
#' compensate,flowSet,ANY-method compensate,flowSet,list-method
#' compensate,flowSet,data.frame-method transform,flowSet-method
#' rbind2,flowSet,missing rbind2,flowSet,flowSet-method
#' rbind2,flowSet,flowSet,missing-method rbind2,flowSet,flowFrame-method
#' rbind2,flowFrame,flowSet-method rbind2,flowSet,missing-method
#' summary,flowSet-method
#' @docType class
#' 
#' @slot frames An \code{\link[base:environment]{environment}}
#' containing one or more \code{\linkS4class{flowFrame}} objects.
#' @slot phenoData An
#' \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#' containing the phenotypic data for the whole data set. Each row
#' corresponds to one of the \code{\linkS4class{flowFrame}}s in the
#' \code{frames} slot.  The \code{sampleNames} of \code{phenoData}
#' (see below) must match the names of the
#' \code{\linkS4class{flowFrame}} in the \code{frames} environment.
#' @slot colnames A \code{character} object with the (common)
#' column names of all the data matrices in the
#' \code{\linkS4class{flowFrame}}s.
#' 
#' @section Creating Objects:
#' 
#' Objects can be created using\cr \code{ new('flowSet',}\cr \code{ frames =
#' ...., # environment with flowFrames}\cr \code{ phenoData = .... # object of
#' class AnnotatedDataFrame}\cr \code{ colnames = ....  # object of class
#' character}\cr \code{ )}\cr
#' 
#' or via the constructor \code{flowSet}, which takes arbitrary numbers of
#' flowFrames, either as a list or directly as arguments, along with an
#' optional \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#' for the \code{phenoData} slot and a \code{character} scalar for the
#' \code{name} by which the object can be referenced.
#' 
#' \code{flowSet(..., phenoData)}
#' 
#' Alternatively, \code{flowSets} can be coerced from \code{list} and
#' \code{environment} objects.
#' 
#' \code{as(list("A"=frameA,"B"=frameB),"flowSet")}
#' 
#' The safest and easiest way to create \code{flowSet}s directly from
#' \acronym{FCS} files is via the \code{\link{read.flowSet}} function, and
#' there are alternative ways to specify the files to read. See the separate
#' documentation for details.
#' 
#' @section Methods:
#'   \describe{
#' 
#' \item{[, [[}{Subsetting. \code{x[i]} where \code{i} is a scalar,
#'   returns a \code{flowSet} object, and \code{x[[i]]} a
#'   \code{\linkS4class{flowFrame}} object. In this respect the
#'   semantics are similar to the behavior of the subsetting operators
#'   for lists. \code{x[i, j]} returns a \code{flowSet} for which the
#'   parameters of each \code{\linkS4class{flowFrame}} have been subset
#'   according to \code{j}, \code{x[[i,j]]} returns the subset of a
#'   single \code{\linkS4class{flowFrame}} for all parameters in
#'   \code{j}. Similar to data frames, valid values for \code{i} and
#'   \code{j} are logicals, integers and characters.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   flowSet[i]}
#'   
#'   \code{   flowSet[i,j]}
#'   
#'   \code{   flowSet[[i]]}
#'   
#' }
#' 
#' \item{$}{Subsetting by frame name. This will return a single
#'   \code{\linkS4class{flowFrame}} object. Note that names may have to
#'   be quoted if they are no valid R symbols
#'   (e.g. \code{flowSet$"sample 1"}}
#' 
#' \item{colnames, colnames<-}{Extract or replace the \code{colnames}
#'   slot.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   colnames(flowSet)}
#'   
#'   \code{   colnames(flowSet) <- value}
#'   
#' }
#' 
#' \item{identifier, identifier<-}{Extract or replace the \code{name}
#'   item from the environment.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   identifier(flowSet)}
#'   
#'   \code{   identifier(flowSet) <- value}
#'   
#' }
#' 
#' 
#' \item{phenoData, phenoData<-}{Extract or replace the
#'   \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#'   from the \code{phenoData} slot.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   phenoData(flowSet)}
#'   
#'   \code{   phenoData(flowSet) <- value}
#'   
#' }
#' 
#' \item{pData, pData<-}{Extract or replace the data frame (or columns
#'                                                          thereof) containing actual phenotypic information from the
#'   \code{phenoData} slot.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   pData(flowSet)}
#'   
#'   \code{   pData(flowSet)$someColumn <- value}
#'   
#' }
#' 
#' \item{varLabels, varLabels<-}{ Extract and set varLabels in the
#'   \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#'   of the \code{phenoData} slot.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   varLabels(flowSet)}
#'   
#'   \code{   varLabels(flowSet) <- value}
#'   
#' }
#' 
#' \item{sampleNames}{Extract and replace sample names from the
#'   \code{phenoData} object. Sample names correspond to frame
#'   identifiers, and replacing them will also replace the \code{GUID}
#'   slot for each frame. Note that \code{sampleName} need to be
#'   unique.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   sampleNames(flowSet)}
#'   
#'   \code{   sampleNames(flowSet) <- value}
#'   
#' }
#' 
#' \item{keyword}{Extract or replace keywords specified in a character
#'   vector or a list from the \code{description} slot of each
#'   frame. See \code{\link{keyword}} for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   keyword(flowSet, list(keywords))}
#'   
#'   \code{   keyword(flowSet, keywords)}
#'   
#'   \code{   keyword(flowSet) <- list(foo="bar") }
#'   
#' }
#' 
#' \item{length}{number of \code{\linkS4class{flowFrame}} objects in
#'   the set.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   length(flowSet)}
#'   
#' }
#' 
#' \item{show}{display object summary.}
#' 
#' \item{summary}{Return descriptive statistical summary (min, max,
#'                                                        mean and quantile) for each channel of each
#'   \code{\linkS4class{flowFrame}}
#'   
#'   \emph{Usage:}
#'   
#'   \code{   summary(flowSet)}
#'   
#' }
#' 
#' 
#' \item{fsApply}{Apply a function on all frames in a \code{flowSet}
#'   object. Similar to \code{\link{sapply}}, but with additional
#'   parameters. See separate documentation for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   fsApply(flowSet, function, ...)}
#'   
#'   \code{   fsApply(flowSet, function, use.exprs=TRUE, ...)}
#'   
#' }
#' 
#' \item{compensate}{Apply a compensation matrix on all frames in a
#'   \code{flowSet} object. See separate documentation for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   compensate(flowSet, matrix)}
#'   
#' }
#' 
#' \item{transform}{Apply a transformation function on all frames of a
#'   \code{flowSet} object. See separate documentation for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   transform(flowSet, ...)}
#'   
#' }
#' 
#' \item{filter}{Apply a filter object on a \code{flowSet}
#'   object. There are methods for \code{\linkS4class{filter}}s,
#'   \code{\link{filterSet}}s and lists of filters. The latter has to
#'   be a named list, where names of the list items are matching
#'   sampleNames of the \code{flowSet}. See \code{\linkS4class{filter}}
#'   for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   filter(flowSet, filter)}
#'   
#'   \code{   filter(flowSet, list(filters))}
#'   
#' }
#' 
#' \item{split}{Split all \code{flowSet} objects according to a
#'   \code{\link{filter}}, \code{\link{filterResult}} or a list of such
#'   objects, where the length of the list has to be the same as the
#'   length of the \code{flowSet}. This returns a list of
#'   \code{\linkS4class{flowFrame}}s or an object of class
#'   \code{flowSet} if the \code{flowSet} argument is set to
#'   \code{TRUE}. Alternatively, a \code{flowSet} can be split into
#'   separate subsets according to a factor (or any vector that can be
#'                                           coerced into factors), similar to the behaviour of
#'   \code{\link[base]{split}} for lists. This will return a list of
#'   \code{flowSet}s. See \code{\link{split}} for details.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   split(flowSet, filter)}
#'   
#'   \code{   split(flowSet, filterResult)}
#'   
#'   \code{   split(flowSet, list(filters))}
#'   
#'   \code{   split(flowSet, factor)}
#'   
#' }
#' 
#' \item{Subset}{Returns a \code{flowSet} of
#'   \code{\linkS4class{flowFrame}}s that have been subset according
#'   to a \code{\linkS4class{filter}} or
#'   \code{\linkS4class{filterResult}}, or according to a list of such
#'   items of equal length as the \code{flowSet}.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   Subset(flowSet, filter)}
#'   
#'   \code{   Subset(flowSet, filterResult)}
#'   
#'   \code{   Subset(flowSet, list(filters))}
#'   
#' }
#' 
#' 
#' \item{rbind2}{Combine two \code{flowSet} objects, or one
#'   \code{flowSet} and one \code{\linkS4class{flowFrame}} object.
#'   
#'   \emph{Usage:}
#'   
#'   \code{   rbind2(flowSet, flowSet)}
#'   
#'   \code{   rbind2(flowSet, flowFrame)}
#'   
#' }
#' 
#' \item{spillover}{Compute spillover matrix from a compensation
#'   set. See separate documentation for details.
#' }
#' }
#' 
#' @section Important note on storage and performance:
#' The bulk of the data in a \code{flowSet} object is stored in an
#' \code{\link[base:environment]{environment}}, and is therefore not
#' automatically copied when the \code{flowSet} object is copied. If
#' \code{x} is an object of class \code{flowSet}, then the code
#' \preformatted{y <- x} will create an object \code{y} that contains
#' copies of the \code{phenoData} and administrative data in \code{x},
#' but refers to the \emph{same} environment with the actual fluorescence
#' data. See below for how to create proper copies.
#' 
#' The reason for this is performance. The pass-by-value semantics of
#' function calls in \code{R} can result in numerous copies of the same
#' data object being made in the course of a series of nested function
#' calls. If the data object is large, this can result in considerable
#' cost of memory and performance. \code{flowSet} objects are intended to
#' contain experimental data in the order of hundreds of Megabytes, which
#' can effectively be treated as read-only: typical tasks are the
#' extraction of subsets and the calculation of summary statistics.  This
#' is afforded by the design of the \code{flowSet} class: an object of
#' that class contains a \code{phenoData} slot, some administrative
#' information, and a \emph{reference} to an environment with the
#' fluorescence data; when it is copied, only the reference is copied,
#' but not the potentially large set of fluorescence data themselves.
#' 
#' However, note that subsetting operations, such as \code{y <- x[i]} do
#' create proper copies, including a copy of the appropriate part of the
#' fluorescence data, as it should be expected. Thus, to make a proper
#' copy of a \code{flowSet} \code{x}, use \code{y <- x[seq(along=x)]}
#' 
#' @author
#' 
#' F. Hahne, B. Ellis, P. Haaland and N. Le Meur
#' @seealso
#' 
#' \code{\linkS4class{flowFrame}}, \code{\link{read.flowSet}}
#' @keywords classes
#' @examples
#' 
#' 
#' ## load example data and object creation
#' data(GvHD)
#' 
#' ## subsetting to flowSet
#' set <- GvHD[1:4]
#' GvHD[1:4,1:2]
#' sel <- sampleNames(GvHD)[1:2]
#' GvHD[sel, "FSC-H"]
#' GvHD[sampleNames(GvHD) == sel[1], colnames(GvHD[1]) == "SSC-H"]
#' 
#' ## subsetting to flowFrame
#' GvHD[[1]]
#' GvHD[[1, 1:3]]
#' GvHD[[1, "FSC-H"]]
#' GvHD[[1, colnames(GvHD[1]) == "SSC-H"]]
#' GvHD$s5a02
#' 
#' ## constructor
#' flowSet(GvHD[[1]], GvHD[[2]])
#' pd <- phenoData(GvHD)[1:2,]
#' flowSet(s5a01=GvHD[[1]], s5a02=GvHD[[2]],phenoData=pd)
#' 
#' ## colnames
#' colnames(set)
#' colnames(set) <- make.names(colnames(set))
#' 
#' ## object name
#' identifier(set)
#' identifier(set) <- "test"
#' 
#' ## phenoData
#' pd <- phenoData(set)
#' pd
#' pd$test <- "test"
#' phenoData(set) <- pd
#' pData(set)
#' varLabels(set)
#' varLabels(set)[6] <- "Foo"
#' varLabels(set)
#' 
#' ## sampleNames
#' sampleNames(set)
#' sampleNames(set) <- LETTERS[1:length(set)]
#' sampleNames(set)
#' 
#' ## keywords
#' keyword(set, list("transformation"))
#' 
#' ## length
#' length(set)
#' 
#' ## compensation
#' samp <- read.flowSet(path=system.file("extdata","compdata","data",
#' package="flowCore"))
#' cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
#' comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
#' comp.mat
#' summary(samp[[1]])
#' samp <- compensate(samp, as.matrix(comp.mat))
#' summary(samp[[1]])
#' 
#' ## transformation
#' opar <- par(mfcol=c(1:2))
#' plot(set[[1]], c("FL1.H", "FL2.H"))
#' set <- transform(set, transformList(c("FL1.H", "FL2.H"), log))
#' plot(set[[1]], c("FL1.H", "FL2.H"))
#' par(opar)
#' 
#' ## filtering of flowSets
#' rectGate <- rectangleGate(filterId="nonDebris", FSC.H=c(200,Inf))
#' fres <- filter(set, rectGate)
#' class(fres)
#' summary(fres[[1]])
#' rectGate2 <- rectangleGate(filterId="nonDebris2", SSC.H=c(300,Inf))
#' fres2 <- filter(set, list(A=rectGate, B=rectGate2, C=rectGate, D=rectGate2))
#' 
#' ## Splitting frames of a flowSet
#' split(set, rectGate)
#' split(set[1:2], rectGate, populatiuon="nonDebris2+")
#' split(set, c(1,1,2,2))
#' 
#' ## subsetting according to filters and filter results
#' Subset(set, rectGate)
#' Subset(set, filter(set, rectGate))
#' Subset(set, list(A=rectGate, B=rectGate2, C=rectGate, D=rectGate2))
#' 
#' ## combining flowSets
#' rbind2(set[1:2], set[3:4])
#' rbind2(set[1:3], set[[4]])
#' rbind2(set[[4]], set[1:2])
#' 
#' 
setClass("flowSet",                   
         representation=representation(frames="environment",
         phenoData="AnnotatedDataFrame",
         colnames="character"),
         prototype=list(frames=new.env(hash=TRUE, parent=emptyenv()),
         phenoData=new("AnnotatedDataFrame",
         data=data.frame(),
         varMetadata=data.frame()),
         colnames=character(0)),
         validity=function(object){
             nc <- length(colnames(object))
             ## Make sure that all of our samples list
             name.check <- is.na(match(sampleNames(object), ls(object@frames,
                                                               all.names=TRUE)))
             if(any(name.check)) {
                 name.list <- paste(sampleNames(object)[name.check], sep=",")
                 return(paste("These objects are not in the data environment:",
                              name.list))
             }
             
             ##Ensure that all frames match our colnames
             if(!all(sapply(sampleNames(object), function(i) {
                 x <- get(i, env=object@frames)
                 
                 if(all(object@colnames == colnames(x))){
                     TRUE
                 }else{ 
                     message(i, " doesn't have the identical colnames as the other samples!")
                   FALSE
                 }
             }))){
                 return(paste("Some items identified in the data environment",
                              "either have the wrong dimension or type."))
             }
             return(TRUE)
         })

## constructor
flowSet <- function(..., phenoData, name)
{
    x <- list(...)
    if(length(x) == 1 && is.list(x[[1]]))
        x <- x[[1]]
    if(!all(sapply(x, is, "flowFrame")))
        stop("All additional arguments must be flowFrames")
    f <- as(x, "flowSet")
    if(!missing(phenoData))
        phenoData(f) <- phenoData
    if(!missing(name))
        identifier(f) <- name
    f
}



## ===========================================================================
## transform parent class and parameters
## ---------------------------------------------------------------------------
## Parameterize transforms so that we can describe them.
## ---------------------------------------------------------------------------
#' 'transform': a class for transforming flow-cytometry data by applying scale
#' factors.
#' 
#' Transform objects are simply functions that have been extended to allow for
#' specialized dispatch. All of the ``...Transform'' constructors return
#' functions of this type for use in one of the transformation modalities.
#' 
#' 
#' @name transform-class
#' @aliases transform transform,missing-method transform-class
#' summary,transform-method show,transform-method
#' @docType class
#' 
#' @slot .Data Object of class \code{"function"}
#' @slot transformationId A name for the transformation
#' object
#' 
#' @section Methods:
#' \describe{
#' \item{\code{summary}}{Return the parameters}
#' }
#' 
#' @author N LeMeur
#' @seealso \code{\link[flowCore]{linearTransform}},
#' \code{\link[flowCore]{lnTransform}},
#' \code{\link[flowCore]{logicleTransform}},
#' \code{\link[flowCore]{biexponentialTransform}},
#' \code{\link[flowCore]{arcsinhTransform}},
#' \code{\link[flowCore]{quadraticTransform}},
#' \code{\link[flowCore]{logTransform}}
#' @keywords classes
#' @examples
#' 
#' 
#' cosTransform <- function(transformId, a=1, b=1){
#'   t = new("transform", .Data = function(x) cos(a*x+b))
#'   t@transformationId = transformId
#'   t
#' }
#' 
#' cosT <- cosTransform(transformId="CosT",a=2,b=1)
#' 
#' summary(cosT)
#' 
setClass("transform",
         representation=representation(transformationId="character",
                                       .Data="function"),
         prototype=prototype(transformationId=""))

setClass("parameters", contains="list")

setClassUnion("transformation", "transform")

setClassUnion("characterOrTransformation", c("character","transformation"))

setClassUnion("characterOrParameters", c("character","parameters"))

setClass("singleParameterTransform",
         representation=representation(parameters="transformation"),
         contains="transform")

setClass("nullParameter",
         representation=representation(dummy="numeric"))



## ===========================================================================
## Virtual filter and derived concreteFilter and parameterFilter
## ---------------------------------------------------------------------------
## A class describing a selection applied to a flow data matrix. Consist of
## a filterId and the names of the parameters to operate on (for parameter
## filters only). More specific filters all inherit from either of these two
## classes.
## ---------------------------------------------------------------------------
#' A class for representing filtering operations to be applied to flow data.
#' 
#' The \code{filter} class is the virtual base class for all filter/gating
#' objects in \code{flowCore}. In general you will want to subclass or create a
#' more specific filter.
#' 
#' 
#' @name filter-class
#' @aliases filter-class filtergate,filter-class rectangleGate,filter-class
#' polygonGate,filter-class ellipsoidGate,filter-class norm2Filter,filter-class
#' decisionTreeGate,filter-class booleanGate,filter-class filter,filter-method
#' |,filter,filter-method !,filter-method |,filter,list-method
#' |,list,filter-method
#' @docType class
#' 
#' @slot filterId A character vector that identifies this \code{filter}. 
#' This is typically user specified but can be automatically deduced by 
#' certain filter operations, particularly boolean and
#' set operations.
#' 
#' @section Objects from the Class:
#' 
#' All \code{\link[flowCore:filter-class]{filter}} objects in \code{flowCore}
#' should be instantiated through their constructors. These are functions
#' that share the same name with the respective \code{filter}
#' classes. E.g.,
#' \code{\link[flowCore:rectangleGate]{rectangleGate()}} is the 
#' constructor function for rectangular gates, and
#' \code{\link[flowCore:kmeansFilter]{kmeansFilter()}} creates
#' objects of class \code{\link{kmeansFilter}}. Usually these
#' constructors can deal with various different inputs, allowing to
#' utilize the same function in different programmatic or interactive
#' settings. For all \code{filters} that operate on specific flow
#' parameters (i.e., those inheriting from 
#'             \code{\link[flowCore:parameterFilter-class]{parameterFilter}}), the parameters
#' need to be passed to the constructor, either as names or colnames of
#' additional input arguments or explicitly as separate arguments.  See
#' the documentation of the respective \code{filter} classes for
#' details. If parameters are explicitly defined as separate arguments,
#' they may be of class \code{character}, in which case they will be
#' evaluated literally as colnames in a \code{\link{flowFrame}}, or of
#' class \code{\link[flowCore:transform-class]{transform}}, in which case the
#' filtering is performed on a temporarily transformed copy of the input
#' data. See \code{\link[flowCore:parameterFilter-class]{here}} for details.
#' 
#' @section Methods:
#' \describe{
#' \item{\code{\%in\%}}{Used in the usual way this returns a vector of
#'   values that identify which events were accepted by the filter. A
#'   single filter may encode several populations so this can return
#'   either a \code{logical} vector, a \code{factor} vector or a
#'   \code{numeric} vector of probabilities that the event is accepted
#'   by the filter. Minimally, you must implement this method when
#'   creating a new type of filter}
#' 
#' \item{\code{&}, \code{|}, \code{!}}{Two filters can be composed
#'   using the usual boolean operations returning a \code{filter} class
#'   of a type appropriate for handling the operation. These methods
#'   attempt to guess an appropriate \code{filterId} for the new
#'   \code{filter}}
#' 
#' \item{\code{\%subset\%}, \code{\%&\%}}{Defines a filter as being a
#'   subset of another filter. For deterministic filters the results
#'   will typically be equivalent to using an \code{\&} operation to
#'   compose the two filters, though summary methods will use subset
#'   semantics when calculating proportions. Additionally, when the
#'   filter is data driven, such as
#'   \code{\link[flowCore:norm2Filter-class]{norm2Filter}}, the subset
#'   semantics are 
#'   applied to the data used to fit the filter possibly resulting in
#'   quite different, and usually more desirable, results.}
#' 
#' \item{\code{\%on\%}}{Used in conjunction with a
#'   \code{\link[flowCore:transformList-class]{transformList}} to create a
#'   \code{transformFilter}. This filter is similar to the subset
#'   filter in that the filtering operation takes place on transformed
#'   values rather than the original values.}
#' 
#' \item{\code{filter}}{A more formal version of \code{\%in\%}, this
#'   method returns a
#'   \code{\link[flowCore:filterResult-class]{filterResult}} object
#'   that can be used in subsequent filter operations as well as providing
#'   more metadata about the results of the filtering operation}
#' 
#' \item{\code{summarizeFilter}}{When implementing a new filter this
#'   method is used to update the \code{filterDetails} slot of a
#'   \code{filterResult}. It is optional and typically only needs to be
#'   implemented for data-driven filters.}
#' 
#' }
#' 
#' @author B. Ellis, P.D. Haaland and N. LeMeur
#' @seealso \code{\link[flowCore:transform-class]{transform}},
#' \code{\link[flowCore:filter-methods]{filter}}
#' @keywords classes
setClass("filter", 
         representation=representation("VIRTUAL",
         filterId="character"),
         prototype=prototype(filterId=""))

setClass("concreteFilter",
         contains="filter")

                                        # setClass("parameterFilter",
                                        #          representation=representation(parameters="character"),
                                        #          contains="concreteFilter",
                                        #          prototype=prototype(parameters=""))

setClass("parameterFilter", 
         representation(parameters="parameters"),
         contains="concreteFilter",
         prototype=prototype(parameters=new("parameters",.Data="NULL"))
         )

#########################################################################
#filters is a list of filters for the same flowFrame
#thus is different from filerList which is for a flowSet
#---------------------------------------------------------------------
#which are supposed to be gated on the same parent population.
#It is mainly for plotting multiple gates per flowFramein flowViz::xyplot.
#These gates should have the same parameters(channels)
###########################################################################
#' Class "filters" and "filtersList"
#' 
#' The \code{filters} class is the container for a list of
#' \code{\link[flowCore:filter-methods]{filter}} objects.\cr\cr
#' The \code{filtersList}
#' class is the container for a list of \code{filters} objects. 
#' 
#' The \code{filters} class mainly
#' exists for displaying multiple filters/gates on one single panel(flowFrame)
#' of \code{\link[flowViz:xyplot]{xyplot}}. Note that it is different from
#' \code{\link[flowCore:filterList]{filterList}} class which is to be applied to
#' a flowSet. In other words, \code{filter} objects of a \code{fliterList} are
#' to be applied to different flowFrames. However,all of \code{filter} objects
#' of a \code{filters} object are for one single flowFrame, more specifically for one
#' pair of projections(parameters).So these filters should share the common
#' parameters.\cr\cr
#' And \code{filtersList} is a list of \code{filters} objects, which are to be
#' applied to a flowSet.
#' 
#' 
#' @name filters-class
#' @aliases filters-class filters filtersList-class filtersList
#' show,filters-method show,filtersList-method
#' @docType class
#' 
#' @usage 
#' \code{filters(x)}\cr
#' \code{filtersList(x)}
#' 
#' @param   x A list of \code{filter} or \code{filters} objects.
#' 
#' @return  A \code{filters} or \code{filtersList} object from the constructor 
#' 
#' @slot .Data Object of class
#' \code{"list"}. The class directly extends \code{list}, and this slot holds
#' the list data.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{list}"}
#' 
#' @section Objects from the Class:
#' Objects are created from regular lists using the constructors 
#' \code{filters} and \code{filtersList}:
#' 
#' \code{filters(x)}
#' 
#' \code{filtersList(x)}
#' 
#' @author Mike Jiang
#' @seealso \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:filterList-class]{filterList}}
#' @keywords classes
setClass("filters",
		 contains="list"
		 )
 ## Constructor
 filters <- function(x)
 {
	 checkClass(x, "list")
	 x <- new("filters", .Data=x)
	 validFilters(x)
	 return(x)
 }
 #' Check if all filters in a filters matches same paramters
 #' @param flist a filters object
 #' @return TRUE or FALSE
 validFilters<- function(flist)
 {
	 res <- TRUE
	 checkClass(flist, "filters")
		
     
     fParams <- lapply(flist, function(x) sort(parameters(x)))
     nParam <- length(unique(fParams))
     
     valid <- FALSE
     #validity check for 1d gate (nParam up to 2 is allowed)
     if(all(sapply(fParams, length) == 1)){
       valid <- nParam <= 2
     }else{
       #otherwise consider them as 2d
      valid <- nParam == 1 
     }
     if(!valid)
     {
       stop("Not all filter objects in the list have the same paramters", call.=FALSE)
       res <- FALSE
     }  
     
	 
	 if(any(sapply(flist, is, "filterResult")))
	 {
		 stop("filterResults are not allowed in a filterList") 
		 res <- FALSE
	 }
	 return(res)
 
 }
		 
		
		 
#########################################################################
#filtersList is a list filters to be applied to a flowSet
#---------------------------------------------------------------------
 
setClass("filtersList",
		 contains="list"
		 )
		 
 ## Check if a filtersList matches a flowSet.
 validFiltersList <- function(flist, set, strict=TRUE)
 {
	 res <- TRUE
	 checkClass(flist, "filtersList")
	 checkClass(strict, "logical", 1)
	 if(!missing(set)){
		 checkClass(set, "flowSet")
		 if(res <- !all(names(flist) == sampleNames(set)))
			 warning("Sample names don't match between flowSet and ",
					 "filterResultList", call.=FALSE)
	 }
	 
	 if(strict){
		 fTypes <- unname(sapply(flist, class,simplify=F))
		 if(length(unique(fTypes)) != 1)
		 {
			 warning("Not all filter objects in the list are of equal",
					 " type.", call.=FALSE)
			 res <- FALSE
		 }
		 if(any(sapply(flist, is, "filterResult")))
		 {
			 stop("filterResults are not allowed in a filterList") 
				 res <- FALSE
		 }
		 return(res)
	}
 }
 
 ## Constructor
 filtersList <- function(x)
 {
	 checkClass(x, "list")
	 
	 if(is.null(names(x)))
		 stop("Names missing in input list.")
	 x <- new("filtersList", .Data=x)
	 validFiltersList(x)
	 return(x)
 }
## ===========================================================================
## Rectangular gate
## ---------------------------------------------------------------------------
## A class describing a 2D rectangular region in the parameter space. Slots
## min and max hold the boundaries in the two dimensions.
## ---------------------------------------------------------------------------
#' Class "rectangleGate"
#' 
#' 
#' Class and constructor for n-dimensional rectangular
#' \code{\linkS4class{filter}} objects.
#' 
#' 
#' This class describes a rectangular region in n dimensions, which is a
#' Cartesian product of \code{n} orthogonal intervals in these dimensions.
#' \code{n=1} corresponds to a range gate, \code{n=2} to a rectangle gate,
#' \code{n=3} corresponds to a box region and \code{n>3} to a hyper-rectangular
#' regions. Intervals may be open on one side, in which case the value for the
#' boundary is supposed to be \code{Inf} or \code{-Inf}, respectively.
#' \code{rectangleGates} are inclusive, that means that events on the
#' boundaries are considered to be in the gate.
#' 
#' The constructor is designed to be useful in both direct and programmatic
#' usage. To use it programmatically, you may either construct a named list or
#' you may construct a matrix with \code{n} columns and \code{2} rows.  The
#' first row corresponds to the minimal value for each parameter while the
#' second row corresponds to the maximal value for each parameter.  The names
#' of the parameters are taken from the column names or from the list names,
#' respectively. Alternatively, the boundaries of the \code{rectangleGate} can
#' be given as additional named arguments, where each of these arguments should
#' be a numeric vector of length \code{2}; the function tries to collapse these
#' boundary values into a matrix.
#' 
#' Note that boundaries of \code{rectangleGates} where \code{min > max} are
#' syntactically valid, however when evaluated they will always be empty.
#' 
#' \code{rectangleGate} objects can also be multiplied using the \code{*}
#' operator, provided that both gates have orthogonal axes. This results in
#' higher-dimensional \code{rectangleGates}. The inverse operation of
#' subsetting by parameter name(s) is also available.
#' 
#' Evaluating a \code{rectangleGate} generates an object of class
#' \code{\linkS4class{logicalFilterResult}}. Accordingly, \code{rectangleGates}
#' can be used to subset and to split flow cytometry data sets.
#' 
#' @name rectangleGate-class
#' @aliases rectangleGate-class rectangleGate summary,rectangleGate-method
#' show,rectangleGate-method [,rectangleGate,character-method
#' [,rectangleGate,ANY-method *,rectangleGate,rectangleGate-method
#' @docType class
#' 
#' 
#' @usage rectangleGate(\dots, .gate, filterId="defaultRectangleGate")
#' 
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' gate. The object can later be identified by this name.
#' @param .gate A definition of the gate. This can be either a list, or a
#' matrix, as described below.
#' @param \dots You can also directly provide the boundaries of a
#' \code{rectangleGate} as additional named arguments, as described below.
#' @return
#' 
#' Returns a \code{\link{rectangleGate}} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @note
#' 
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for details on plotting of \code{rectangleGates}.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{parameterFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, by class
#' \code{parameterFilter}, distance 2.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
#' distance 3.
#' 
#' @slot min,max Objects of class \code{"numeric"}. The
#' minimum and maximum values of the n-dimensional rectangular
#' region.
#' 
#' @slot parameters Object of class \code{"character"},
#' indicating the parameters for which the \code{rectangleGate} is
#' defined.
#' 
#' @slot filterId Object of class \code{"character"},
#' referencing the filter.
#' 
#' @section Objects from the Class:
#' 
#' Objects can be created by calls of the form \code{new("rectangleGate",
#' ...)}, by using the constructor \code{rectangleGate} or by combining
#' existing \code{rectangleGates} using the \code{*} method.  Using the
#' constructor is the recommended way of object instantiation.
#' 
#' @section Methods:
#' \describe{
#'    \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'      "rectangleGate")}: The workhorse used to evaluate the filter on
#'      data. This is usually not called directly by the user, but
#'      internally by calls to the \code{\link{filter}} methods. }
#'    
#'    \item{show}{\code{signature(object = "rectangleGate")}: Print
#'      information about the filter. }
#'    
#'    \item{*}{\code{signature(e1 = "rectangleGate", e2 =
#'      "rectangleGate")}: combining two \code{rectangleGates} into one
#'      higher dimensional representation. }
#'    
#'    \item{[}{\code{signature(x = "rectangleGate", i = "character")}:
#'      Subsetting of a \code{rectangleGate} by parameter name(s). This
#'      is essentially the inverse to \code{*}. }
#' }
#' 
#' @author F.Hahne, B. Ellis N. Le Meur
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{polygonGate}},
#' \code{\link{ellipsoidGate}}, \code{\link{polytopeGate}},
#' \code{\link{filter}} for evaluation of \code{rectangleGates} and
#' \code{\link{split}} and \code{\link{Subset}}for splitting and subsetting of
#' flow cytometry data sets based on that.
#' @keywords methods classes
#' @examples
#' 
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' #Create directly. Most likely from a command line
#' rectangleGate(filterId="myRectGate", "FSC-H"=c(200, 600), "SSC-H"=c(0, 400))
#' 
#' #To facilitate programmatic construction we also have the following
#' rg <- rectangleGate(filterId="myRectGate", list("FSC-H"=c(200, 600),
#' "SSC-H"=c(0, 400)))
#' mat <- matrix(c(200, 600, 0, 400), ncol=2, dimnames=list(c("min", "max"),
#' c("FSC-H", "SSC-H")))
#' rg <- rectangleGate(filterId="myRectGate", .gate=mat)
#' 
#' ## Filtering using rectangleGates
#' fres <- filter(dat, rg)
#' fres
#' summary(fres)
#' 
#' ## The result of rectangle filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' ## Multiply rectangle gates
#' rg1 <- rectangleGate(filterId="FSC-", "FSC-H"=c(-Inf, 50))
#' rg2 <- rectangleGate(filterId="SSC+", "SSC-H"=c(50, Inf))
#' rg1 * rg2
#' 
#' ## Subset rectangle gates
#' rg["FSC-H"]
#' 
#' ##2d rectangleGate can be coerced to polygonGate
#' as(rg, "polygonGate")
#' 
#' 
setClass("rectangleGate",
         representation=representation(min="numeric",
         max="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultRectangleGate",
         min=-Inf,
         max=Inf)
         )

## parse '...' argument of a gate constructor. The return value is a list
## with parameters (as transforms) and values.
parseDots <- function(dl, collapseFirst=TRUE, len=NULL){
    parseItem <- function(i, x, len){
        ## We can return transforms directly
        y <- x[[i]]
        if(is(y, "transform")){
            dl[[i]] <<- NA
            y
        }else{
            li <- length(y)
            if(!is.character(y) && !is.null(len) && li!=len)
                stop("All additional arguments must be of length ",
                     len, call.=FALSE)
            if(!is.character(y) && li!=allLen)
                stop("All additional arguments must be of equal length ",
                     call.=FALSE)     
            if(!is.character(y) && is.null(names(x)[i]))
                stop("Additional arguments have to be named.",
                     call.=FALSE)
            if(is.character(y)){
                ## We return character scalars as unitytransforms
                dl[[i]] <<- NA
                unitytransform(y) 
            }else{
                ## For eerything else we make unitytransforms from the
                ## argument names
                unitytransform(names(x)[i])
            }
        }
    }
    ## We only parse ..1 if it is a list and drop all other arguments
    if(collapseFirst && length(dl) && is.list(dl[[1]]))
        dl <- dl[[1]]
    if(length(dl)){
        ## If ..1 is a character vector we return unitytransforms only
        if(is.character(dl[[1]]) && length(dl[[1]])>1)
            return(list(parameters=sapply(dl[[1]], unitytransform,
                        simplify=FALSE),
                        values=as.list(rep(NA, length(dl[[1]])))))
        ## If ..1 is a matrix we return unitytransforms and the matrix
        if(is.matrix(dl[[1]])){
            if(is.null(colnames(dl[[1]])))
                stop("Matrix of gate boundaries must have colnames.",
                     call.=FALSE)
            return(list(parameters=sapply(colnames(dl[[1]]), unitytransform,
                        simplify=FALSE), values=dl[[1]]))
        }
        ## All items in dl must be of equal length
        allLen <- if(is.character(dl[[1]])) length(dl[[min(length(dl), 2)]]) 
        else length(dl[[1]])
       
    }
    parms <- sapply(seq_along(dl), parseItem, dl, len, simplify=FALSE)
    return(list(parameters=parms, values=dl))
}

## Further process the output of parseDots to collapse individual arguments
prepareInputs <- function(parsed, .gate, ...)
{
    parms <- parsed$parameters
    values <- parsed$values
    if(missing(.gate)){
        if(any(sapply(values, is.na)))
            stop("The gate boundaries has to be provides as argument",
                 " '.gate'", call.=FALSE)
        if(!is.matrix(values)){
            sel <- sapply(values, is, "numeric")
            if(any(sel)){
                values <- matrix(sapply(values[sel], function(x){
                    if(length(x) ==2)
                        x <- sort(x)
                    x}), ncol=length(parms))
                parms <- parms[sel]
                colnames(values) <- sapply(parms, parameters)
            }
            return(list(parameters=parms, values=values))
        }
        if(!length(parms))
            stop("No arguments provided.", call.=FALSE)
        return(parsed)
    }else{
        if(is.matrix(.gate) && !is.null(colnames(.gate)))
            return(parseDots(list(.gate), ...))
        if(any(sapply(values, is.na))){
            if(ncol(.gate) != length(parms))
                stop("Number of parameters and dimensions of supplied",
                     " gate boundaries don't match.", call.=FALSE)
            return(list(parameters=parms, values=.gate))
        }
        if(!length(parms) || !all(sapply(parms, is, "unityTranform"))){
            return(prepareInputs(parseDots(list(.gate), ...)))
        }else{
            return(parsed)
        }
    }
}


## Constructor. We allow for the following inputs:
##  ... are named numerics, each of length 2
##  ... are transforms or a mix of transforms and characters, .gate is
##      the associated matrix of min and max values
##  ..1 is a named list of numerics
##  ..1 is a list of transformations or characters and .gate is the
##      associated matrix of min and max values, each of length 2
##  .gate is a matrix of min and max values with colnames = parameters
##  .gate is a named list of numerics, each of lenght 2
rectangleGate <- function(..., .gate, filterId="defaultRectangleGate")
{
    checkClass(filterId, "character", 1)
    parms <- parseDots(list(...), len=2)
    parms <- prepareInputs(parms, .gate, len=2)
    parms$values <- apply(parms$values, 2, sort)
    new("rectangleGate", filterId = filterId, parameters=parms$parameters,
        min=parms$value[1, ], max=parms$value[2, ])
}



## ===========================================================================
## Quadrant gate
## ---------------------------------------------------------------------------
## A class describing a gate which separates a 2D parameter space into
## four quadrants. Slot boundary holds a vector of length two indicating
## the quadrant boundaries in each of the two dimensions.
## ---------------------------------------------------------------------------
#' Class "quadGate"
#' 
#' 
#' Class and constructors for quadrant-type \code{\link{filter}} objects.
#' 
#' 
#' \code{quadGates} are defined by two parameters, which specify a separation
#' of a two-dimensional parameter space into four quadrants. The
#' \code{quadGate} function is designed to be useful in both direct and
#' programmatic usage.
#' 
#' For the interactive use, these parameters can be given as additional named
#' function arguments, where the names correspond to valid parameter names in a
#' \code{\link{flowFrame}} or \code{\link{flowSet}}. For a more programmatic
#' approach, a named list or numeric vector of the gate boundaries can be
#' passed on to the function as argument \code{.gate}.
#' 
#' Evaluating a \code{quadGate} results in four sub-populations, and hence in
#' an object of class \code{\link{multipleFilterResult}}. Accordingly,
#' \code{quadGates} can be used to split flow cytometry data sets.
#' 
#' @name quadGate-class
#' @aliases quadGate-class quadGate show,quadGate-method
#' @docType class
#' 
#' @usage quadGate(\dots, .gate, filterId="defaultQuadGate")
#' 
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' \code{\link{filter}}. The object can later be identified by this name.
#' @param .gate A definition of the gate for programmatic access. This can be
#' either a named list or a named numeric vector, as described below.
#' @param \dots The parameters of \code{quadGates} can also be directly
#' described using named function arguments, as described below.
#' @return
#' 
#' Returns a \code{quadGate} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @note
#' 
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{quadGates}.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{parameterFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, by class
#' \code{parameterFilter}, distance 2.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
#' distance 3.
#' 
#' @slot boundary Object of class \code{"numeric"}, length
#' 2. The boundaries of the quadrant regions.
#' @slot parameters Object of class \code{"character"},
#' describing the parameter used to filter the \code{flowFrame}.
#' @slot filterId Object of class \code{"character"},
#' referencing the gate.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("quadGate",
#' ...)} or using the constructor \code{quadGate}. The latter is the
#' recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "quadGate")}: The workhorse used to evaluate the gate on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "quadGate")}: Print
#'     information about the gate. }
#'   
#' }
#' 
#' @author F.Hahne, B. Ellis N. Le Meur
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{flowSet}}, \code{\link{filter}} for
#' evaluation of \code{quadGates} and \code{\link{split}} for splitting of flow
#' cytometry data sets based on that.
#' @keywords classes methods
#' @examples
#' 
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Create directly. Most likely from a command line
#' quadGate(filterId="myQuadGate1", "FSC-H"=100, "SSC-H"=400)
#' 
#' ## To facilitate programmatic construction we also have the following
#' quadGate(filterId="myQuadGate2", list("FSC-H"=100, "SSC-H"=400))
#' ## FIXME: Do we want this?
#' ##quadGate(filterId="myQuadGate3", .gate=c("FSC-H"=100, "SSC-H"=400))
#' 
#' ## Filtering using quadGates
#' qg <- quadGate(filterId="quad", "FSC-H"=600, "SSC-H"=400)
#' fres <- filter(dat, qg)
#' fres
#' summary(fres)
#' names(fres)
#' 
#' ## The result of quadGate filtering are multiple sub-populations
#' ## and we can split our data set accordingly
#' split(dat, fres)
#' 
#' ## We can limit the splitting to one or several sub-populations
#' split(dat, fres, population="FSC-H-SSC-H-")
#' split(dat, fres, population=list(keep=c("FSC-H-SSC-H-",
#' "FSC-H-SSC-H+")))
#' 
#' 
setClass("quadGate",
         representation=representation(boundary="numeric"),        
         contains="parameterFilter",
         prototype=list(filterId="defaultQuadGate",
         boundary=c(Inf, Inf)))

## Constructor. We allow for the following inputs:
##  ..1 and ..2 are named numerics of length 1
##  ..1 and ..2 are transforms or a mix of transforms and characters, .gate
##      is the associated numeric vector of boundary values of length 2
##  ..1 is a named list of numerics of length 1
##  ..1 is a list of transformations or characters and .gate is the
##      associated numeric vector of boundary values of length 2
##  .gate is a named list of numerics, each of lenght 1
quadGate <- function(..., .gate, filterId="defaultQuadGate")
{
    checkClass(filterId, "character", 1)
    if(!missing(.gate) && !is.list(.gate) && !is.matrix(.gate))
        .gate <- matrix(.gate, nrow=1, dimnames=list(NULL, names(.gate)))
    parms <- prepareInputs(parseDots(list(...), len=1), .gate, len=1)
    p <- as.numeric(parms$values)
    names(p) <- colnames(parms$values)
    if(length(parms$parameters) !=2 || nrow(parms$value)!=1)
        stop("Expecting two named arguments or a single named vector\n",
             "of length 2 as input for gate boundaries.", call.=FALSE)
    new("quadGate", filterId=filterId, parameters=parms$parameters,
        boundary=p)
}



## ===========================================================================
## Polygon gate
## ---------------------------------------------------------------------------
## A class describing a 2D polygonal region in the parameter space. Slot
## boundary holds the vertices of the polygon in a 2 colum matrix.
## ---------------------------------------------------------------------------
#' Class "polygonGate"
#' 
#' 
#' Class and constructor for 2-dimensional polygonal \code{\link{filter}}
#' objects.
#' 
#' 
#' Polygons are specified by the coordinates of their vertices in two
#' dimensions. The constructor is designed to be useful in both direct and
#' programmatic usage. It takes either a list or a named matrix with \code{2}
#' columns and at least \code{3} rows containing these coordinates.
#' Alternatively, vertices can be given as named arguments, in which case the
#' function tries to convert the values into a matrix.
#' 
#' @name polygonGate-class
#' @aliases polygonGate-class polygonGate show,polygonGate-method
#' @docType class
#' 
#' @usage polygonGate(\dots, .gate, boundaries, filterId="defaultPolygonGate")
#' 
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' gate.
#' @param .gate,boundaries A definition of the gate. This can be either a list
#' or a named matrix as described below. Note the argument boundaries is
#' deprecated and will go away in the next release.
#' @param \dots You can also directly describe a gate without wrapping it in a
#' list or matrix, as described below.
#' @return
#' 
#' Returns a \code{\link{polygonGate}} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @note
#' 
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{polygonGates}.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{parameterFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, by class
#' \code{parameterFilter}, distance 2.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
#' distance 3.
#' 
#' @slot boundaries Object of class \code{"matrix"}. The
#' vertices of the polygon in two dimensions. There need to be at
#' least 3 vertices specified for a valid polygon.
#' @slot parameters Object of class \code{"character"},
#' describing the parameter used to filter the \code{flowFrame}.
#' @slot filterId Object of class \code{"character"},
#' referencing the filter.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("polygonGate",
#' ...)} or by using the constructor \code{polygonGate}. Using the
#' constructor is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "polygonGate")}: The workhorse used to evaluate the filter on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "polygonGate")}: Print
#'     information about the filter. }
#'   
#' }
#' 
#' @author F.Hahne, B. Ellis N. Le Meur
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{rectangleGate}},
#' \code{\link{ellipsoidGate}}, \code{\link{polytopeGate}},
#' \code{\link{filter}} for evaluation of \code{rectangleGates} and
#' \code{\link{split}} and \code{\link{Subset}}for splitting and subsetting of
#' flow cytometry data sets based on that.
#' @keywords methods
#' @examples
#' 
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Defining the gate
#' sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
#' colnames(sqrcut) <- c("FSC-H","SSC-H")
#' pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
#' pg
#' 
#' ## Filtering using polygonGates
#' fres <- filter(dat, pg)
#' fres
#' summary(fres)
#' 
#' ## The result of polygon filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' 
setClass("polygonGate",
         representation(boundaries="matrix"),
         contains="parameterFilter",
         prototype=list(filterId="defaultPolygonGate",
         boundaries=matrix(ncol=2, nrow=3)),
         validity=function(object)
     {
         msg <- TRUE
         if(!is.matrix(object@boundaries) || nrow(object@boundaries)<3 ||
            ncol(object@boundaries)!=2
            )
             msg <- paste("\nslot 'boundaries' must be a numeric matrix",
                          "of at least 3 rows and exactly 2 columns")
         return(msg)
     })

## Constructor. We allow for the following inputs:
##  ..1 and ..2 are named numerics, each of the same length
##  ..1 and ..2  are transforms or a mix of transforms and characters, .gate is
##      the associated matrix of polygon vertices of ncol=2
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and .gate is the
##      associated matrix of polygon vertices of ncol=2
##  .gate is a matrix of polygon vertices of ncol=2, colnames = parameters
##  .gate is a named list of two numerics, both of the same length
polygonGate <- function(..., .gate, boundaries, filterId="defaultPolygonGate")
{
    checkClass(filterId, "character", 1)
    if(missing(.gate))
        if(!missing(boundaries)){
            .Deprecated(msg=paste("The 'boundaries' argument is deprecated,",
                        "please use '.gate' instead."))
            .gate=boundaries
        }     
    parms <- prepareInputs(parseDots(list(...)), .gate)
    if(length(parms$parameters) !=2)
        stop("Polygon gates are only defined in two dimensions.",
             call.=FALSE)
    new("polygonGate", filterId=filterId, parameters=parms$parameters,
        boundaries=parms$values)
}



## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
## A class describing a nD polytope region in the parameter space. Slot a
## holds the coefficients of the linear equations for m halfspaces in n
## dimensions and b is a vector of m intercepts.
## ---------------------------------------------------------------------------
#' Define filter boundaries
#' 
#' 
#' Convenience methods to facilitate the construction of \code{\link{filter}}
#' objects
#' 
#' 
#' These functions are designed to be useful in both direct and programmatic
#' usage.
#' 
#' For rectangle gate in n dimensions, if n=1 the gate correspond to a range
#' gate. If n=2, the gate is a rectangle gate. To use this function
#' programmatically, you may either construct a list or you may construct a
#' matrix with \code{n} columns and \code{2} rows.  The first row corresponds
#' to the minimal value for each parameter while the second row corresponds to
#' the maximal value for each parameter.  The names of the parameters are taken
#' from the column names as in the third example.
#' 
#' Rectangle gate objects can also be multiplied together using the \code{*}
#' operator, provided that both gate have orthogonal axes.
#' 
#' For polygon gate, the boundaries are specified as vertices in 2 dimensions,
#' for polytope gate objects as vertices in n dimensions. 
#' 
#' Polytope gate objects will represent the convex polytope determined
#' by the vertices and parameter b which together specify the polytope as 
#' an intersection of half-spaces represented as a system of linear inequalities,
#' \eqn{Ax\le b}
#' 
#' For quadrant gates, the boundaries are specified as a named list or vector
#' of length two.
#' 
#' 
#' @name polytopeGate-class
#' @aliases polytopeGate-class polytopeGate show,polytopeGate-method
#' @docType class
#' 
#' @usage polytopeGate(\dots, .gate, b, filterId="defaultPolytopeGate")
#' 
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' gate.
#' @param .gate A definition of the gate. This can be either a list, vector or
#' matrix, described below.
#' @param b Need documentation
#' @param \dots You can also directly describe a gate without wrapping it in a
#' list or matrix, as described below.
#' @return
#' 
#' Returns a \code{\link{rectangleGate}} or \code{\link{polygonGate}} object
#' for use in filtering \code{\link{flowFrame}}s or other flow cytometry
#' objects.
#' @author F.Hahne, B. Ellis N. Le Meur
#' @seealso \code{\link{flowFrame}}, \code{\link{filter}}
#' @keywords methods
setClass("polytopeGate",
         representation(a="matrix",b="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultPolytopeGate", a=matrix(), b=1))

## Constructor. We allow for the following inputs:
##  b is always a numeric of length = ncol(a)
##  ... are named numerics, each of the same length
##  ...  are transforms or a mix of transforms and characters, a is
##      the associated matrix of coefficients
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and a is the
##      associated matrix of coefficients
##  .gate is a matrix of coefficients , colnames = parameters
##  .gate is a named list of numerics, all of the same length
polytopeGate <- function(..., .gate, b, filterId="defaultPolytopeGate")
{
    checkClass(filterId, "character", 1)
    checkClass(b, "numeric")
    parms <- prepareInputs(parseDots(list(...)), .gate)
    colnames(parms$values) <- sapply(parms$parameters, parameters)
    new("polytopeGate", filterId=filterId, parameters=parms$parameters,
        a=parms$values, b=b)
}



## ===========================================================================
## Ellipsoid gate
## ---------------------------------------------------------------------------
## A class describing an ellipsoid region in the parameter space. Slots
## mean and cov contain the mean values and the covariance matrix describing
## the ellipse, slot distance holds a scaling factor, i.e., the Mahalanobis
## distance.
## ---------------------------------------------------------------------------
#' Class "ellipsoidGate"
#' 
#' 
#' Class and constructor for n-dimensional ellipsoidal \code{\link{filter}}
#' objects.
#' 
#' 
#' A convenience method to facilitate the construction of a ellipsoid
#' \code{\link{filter}} objects. Ellipsoid gates in n dimensions (n >= 2) are
#' specified by a a covarinace matrix and a vector of mean values giving the
#' center of the ellipse.
#' 
#' This function is designed to be useful in both direct and programmatic
#' usage. In the first case, simply describe the covariance matrix through
#' named arguments. To use this function programmatically, you may pass a
#' covarince matrix and a mean vector directly, in which case the parameter
#' names are the colnames of the matrix.
#' 
#' @name ellipsoidGate-class
#' @aliases ellipsoidGate-class ellipsoidGate show,ellipsoidGate-method
#' @docType class
#' @usage
#' ellipsoidGate(\dots, .gate, mean, distance=1, filterId="defaultEllipsoidGate")
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' gate.
#' @param .gate A definition of the gate via a covariance matrix.
#' @param mean Numeric vector of equal length as dimensions in \code{.gate}.
#' @param distance Numeric scalar giving the Mahalanobis distance defining the
#' size of the ellipse. This mostly exists for compliance reasons to the
#' gatingML standard as \code{mean} and \code{gate} should already uniquely
#' define the ellipse. Essentially, \code{distance} is merely a factor that
#' gets applied to the values in the covariance matrix.
#' @param \dots You can also directly describe the covariance matrix through
#' named arguments, as described below.
#' @return
#' Returns a \code{\link{ellipsoidGate}} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{parameterFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, by class
#' \code{parameterFilter}, distance 2.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
#' distance 3.
#' 
#' @slot mean Objects of class \code{"numeric"}. Vector giving
#' the location of the center of the ellipse in n dimensions.
#' @slot cov Objects of class \code{"matrix"}. The covariance
#' matrix defining the shape of the ellipse.
#' @slot distance Objects of class \code{"numeric"}. The
#' Mahalanobis distance defining the size of the ellipse.
#' @slot parameters Object of class \code{"character"},
#' describing the parameter used to filter the \code{flowFrame}.
#' @slot filterId Object of class \code{"character"},
#' referencing the filter.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("ellipsoidGate",
#' ...)} or by using the constructor \code{ellipsoidGate}.  Using the
#' constructor is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'     \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                     "ellipsoidGate")}: The workhorse used to evaluate the filter on
#'       data. This is usually not called directly by the user, but
#'       internally by calls to the \code{\link{filter}} methods. }
#' 
#'     \item{show}{\code{signature(object = "ellipsoidGate")}: Print
#'      information about the filter. }
#' }
#' @note
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{ellipsoidGates}.
#' 
#' @author F.Hahne, B. Ellis, N. LeMeur
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{polygonGate}},
#' \code{\link{rectangleGate}}, \code{\link{polytopeGate}},
#' \code{\link{filter}} for evaluation of \code{rectangleGates} and
#' \code{\link{split}} and \code{\link{Subset}}for splitting and subsetting of
#' flow cytometry data sets based on that.
#' @keywords methods
#' @examples
#' 
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Defining the gate
#' cov <- matrix(c(6879, 3612, 3612, 5215), ncol=2,
#' dimnames=list(c("FSC-H", "SSC-H"), c("FSC-H", "SSC-H")))
#' mean <- c("FSC-H"=430, "SSC-H"=175)
#' eg <- ellipsoidGate(filterId= "myEllipsoidGate", .gate=cov, mean=mean)
#' 
#' ## Filtering using ellipsoidGates
#' fres <- filter(dat, eg)
#' fres
#' summary(fres)
#' 
#' ## The result of ellipsoid filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' ##ellipsoidGate can be converted to polygonGate by interpolation
#' pg <- as(eg, "polygonGate")
#' pg
#' 
#' 
#' 
setClass("ellipsoidGate",
         representation(mean="numeric",
                        cov="matrix",
			distance="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultEllipsoidGate",
         mean=numeric(), cov=matrix(), distance=1),
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@cov) ||
                nrow(object@cov) != ncol(object@cov) ||
                nrow(object@cov) < 2) 
                 msg <- "\nslot 'cov' must be a symmetric matrix of at least 2 rows"
             if(!is.numeric(object@mean) ||
                length(object@mean) != nrow(object@cov))
                 msg <- paste("\nslot 'mean' must be numeric vector of",
                              "same length as dimensions in 'cov'")
             if(!is.numeric(object@distance) ||	length(object@distance)!=1)
                 msg <- "'distance' must be numeric of length 1"      
             return(msg)
         })

## Constructor. We allow for the following inputs:
##  mean always is a numeric of the same length as number of dimensions,
##  distance is always a vector of length 1
##  ... are named numerics, each of the same length
##  ...  are transforms or a mix of transforms and characters, .gate is
##      the associated covariance matrix
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and .gate is the
##      associated covariance matrix
##  .gate is the covariance matrix, colnames=parameters
##  .gate is a named list of numerics, each of the same length
ellipsoidGate <- function(..., .gate, mean, distance=1,
                          filterId="defaultEllipsoidGate")
{
    checkClass(filterId, "character", 1)
    checkClass(mean, "numeric")
    checkClass(distance, "numeric", 1)
    parms <- prepareInputs(parseDots(list(...)), .gate)
    names(mean) <- sapply(parms$parameters, parameters)
    new("ellipsoidGate", filterId=filterId, parameters=parms$parameters,
        cov=parms$values, mean=mean, distance=distance)
}



## ===========================================================================
## norm2Filter
## ---------------------------------------------------------------------------
## A class to describe the fit of a bivariate normal distribution.
## Slot method is a character describing the method used to compute the
## covariance matrix, slot scale.factor holds a numeric representing the
## Mahalanobis distance. Slot transformation holds a list of length
## giving transformations, if applicable that are applied to the data
## before gating. n is the number of points used in the subsampling step.
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         representation=representation(method="character",
         scale.factor="numeric",
         n="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultNorm2Filter",
         scale.factor=1,
         transformation=list(),
         method="covMcd",
         n=50000))

## Constructor. We allow for the following inputs:
##  method is always a character and scale.factor and n both are always
##     numerics, all of length 1
##  x and y are characters of length 1 or a mix of characters and
##     transformations
##  x is a character of length 2 and y is missing
##  x is a list of characters and/or transformations, y is missing
norm2Filter <- function(x, y, method="covMcd", scale.factor=1,
                        n=50000, filterId="defaultNorm2Filter")
{
    checkClass(method, "character", 1)
    checkClass(scale.factor, "numeric", 1)
    checkClass(n, "numeric", 1)
    checkClass(filterId, "character", 1)
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a norm2 gate.")
        if(length(x)>2)
            warning("Only the first two parameters will be used.")
        y=x[[2]]
        x=x[[1]]
    }
    new("norm2Filter", parameters=c(x, y), method=method,
        scale.factor=scale.factor, filterId=filterId, n=n)
}



## ===========================================================================
## kmeansFilter
## ---------------------------------------------------------------------------
## Apply kmeans clustering on a single parameter. The number k of clusters
## is given by the length of the 'populations' slot. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("kmeansFilter",
         representation=representation(populations="character"),
         prototype=list(filterId="defaultKmeansFilter"),
         contains="parameterFilter")

## Constructor. We allow for the following inputs:
##  ..1 is transform and .2 is some vector that can be coerced to character
##  ..1 is some vector that can be coerced to character
kmeansFilter <- function(..., filterId="defaultKmeansFilter")
{
    checkClass(filterId, "character", 1)
    ll <- list(...)
    if(length(ll)){
        n <- names(ll)[1]
        if(is.list(ll[[1]])){
            n <- names(ll[[1]])[1]
            ll[[1]] <- unlist(ll[[1]], recursive=FALSE)
        }
        parameter <- if(is(ll[[1]], "transform")) ll[[1]] else n
        populations <- if(is(parameter, "transform")){
            if(length(ll)==1)
                stop("List of populations needs to be provided as ",
                     "an additional argument.", call.=FALSE) 
            as.character(unlist(ll[[2]]))} else as.character(unlist(ll[[1]]))
    }else{
        stop("No arguments provided.", .call=FALSE)
    }
    new("kmeansFilter", parameters=parameter,
        populations=populations, filterId=filterId)
}


## ===========================================================================
## sampleFilter
## ---------------------------------------------------------------------------
## Sample 'size' rows from a flowFrame. 
## ---------------------------------------------------------------------------
setClass("sampleFilter",
         representation=representation(size="numeric"),
         contains="concreteFilter",
         prototype=list(size=10000, filterId="defaultSampleFilter"))

##Constructor: We allow for the following inputs:
##  size is always a numeric of length 1
sampleFilter <- function(size, filterId="defaultSampleFilter")
{
    checkClass(filterId, "character", 1)
    checkClass(size, "numeric", 1)
    new("sampleFilter", filterId=filterId, size=size)
}


## ===========================================================================
## boundaryFilter
## ---------------------------------------------------------------------------
## Remove events piled up on the margins of a particular channel
## ---------------------------------------------------------------------------
setClass("boundaryFilter",
         representation=representation(tolerance="numeric", side="character"),
         contains="parameterFilter",
         prototype=list(tolerance=.Machine$double.eps, filterId="defaultBoundaryFilter",
         side="both"))

##Constructor: We allow for the following inputs:
##  tolerance is always a numeric of length 1
boundaryFilter <- function(x, tolerance=.Machine$double.eps, side=c("both", "lower", "upper"),
                           filterId="defaultBoundaryFilter")
{
    checkClass(filterId, "character")
    checkClass(tolerance, "numeric")
    side <- rep(match.arg(side), length(x))[1:length(x)]
    tolerance <- rep(tolerance, length(x))[1:length(x)]
    names(tolerance) <- names(side) <- as.character(x)
    new("boundaryFilter", parameters=x, filterId=filterId, tolerance=tolerance,
        side=side)
}




## ===========================================================================
## expressionFilter
## ---------------------------------------------------------------------------
## Let's us encapsulate an expression as a gate. There also is a constructor
## to create the filter from a character representation of the expression
## which is helpful for programmatic use. The args slot can contain additional
## arguments that are passed on to the evaluation environment. deparse stores
## a deparsed version of the expression.
## ---------------------------------------------------------------------------
setClass("expressionFilter",
         representation=representation(expr="expression",
         args="list",
         deparse="character"),
         contains="concreteFilter",
         prototype=list(filterId="defaultExpressionFilter",
         exprs=expression(rep(TRUE, length(get(ls()[1])))),
         args=list(),
         deparse="default"))

## Constructor: We allow for the following inputs:
##  expr is always an expression
##  ... are further arguments to the expression
expressionFilter <- function(expr, ..., filterId="defaultExpressionFilter")
{
    subs <- substitute(expr)
    if(missing(filterId)){
        filterId <- deparse(subs)
        if(length(filterId)>1)
            filterId <- paste(gsub("^ *", "", filterId[2]), "...", sep="")
    }else checkClass(filterId, "character", 1)
    new("expressionFilter", filterId=filterId, expr=as.expression(subs),
        args=list(...), deparse=deparse(subs))
}

## Constructor from a character string: We allow for the following inputs:
##  expr is always a character string
char2ExpressionFilter <- function(expr, ...,
                                  filterId="defaultExpressionFilter")
{
    checkClass(expr, "character", 1)
    subs <- parse(text=expr)
    if(missing(filterId))
        filterId <- expr
    else
        checkClass(filterId, "character", 1)
    new("expressionFilter", filterId=filterId, expr=subs,
        args=list(...), deparse=expr)
}



## ===========================================================================
## timeFilter
## ---------------------------------------------------------------------------
## Detect turbulences and abnormalities in the aquisition of flow data over
## time and gate them out. Argument 'bandwidth' sets the sensitivity, i.e.,
## the amount of local variance of the signal we want to allow. 'binSize'
## controls the size of the bins for the local variance and location
## estimation, 'timeParameter' can be used to explicitely give the paramter
## name of the time parameter (we will make an educated guess if this is not
## given).
## ---------------------------------------------------------------------------
setClass("timeFilter",
         representation=representation(bandwidth="numeric",
         binSize="numeric",
         timeParameter="character"),
         contains="parameterFilter",
         prototype=list(filterId="defaultTimeFilter",
         bandwidth=0.75,
         binSize=NULL,
         timeParameter=NULL))

## Constructor: We allow for the following inputs:
##  bandwidth and binSize are always numerics of lenght 1, timeParameter
##      is always a character of length 1
##  ..1 is a character
##  ..1 is a list of character and/or transformations
##  ... are characters and/or transformations
timeFilter <- function(..., bandwidth=0.75, binSize, timeParameter,
                       filterId="defaultTimeFilter")
{
    checkClass(bandwidth, "numeric", 1)
    
    if(!missing(binSize))
        checkClass(binSize, "numeric", 1)
    else
        binSize <- NULL
    if(!missing(timeParameter))
        checkClass(timeParameter, "character", 1)
    else
        timeParameter <- NULL
    checkClass(filterId, "character", 1)
    parms <- parseDots(list(...))
    new("timeFilter", parameters=parms$parameters,
        bandwidth=bandwidth, binSize=as.numeric(binSize),
        timeParameter=as.character(timeParameter), filterId=filterId)
}



## ===========================================================================
## filterSet
## ---------------------------------------------------------------------------
## Stores a list of filters from a gating sequence as an environment. Although
## this will be kept as part of flowCore we encourage to use the new workflow
## infrastructure instead.
## ---------------------------------------------------------------------------
setClass("filterSet",
         representation=representation(env="environment",
         name="character"),
         prototype=prototype(env=new.env(hash=TRUE, parent=emptyenv()),
         name="Filter Set"))

## constructor
filterSet <- function(..., name="default") {
    .Deprecated("flowWorkspace::GatingSet")
    filters <- list(...)
    ## Allow the list(x, y, z) format as well.
    if(length(filters)==1 && is.list(filters[[1]]))
        filters <- filters[[1]]
    if(length(filters) == 0)
        new("filterSet", env=new.env(parent=emptyenv()), name=name)
    else{
        tmp <- as(filters, "filterSet")
        tmp@name <- name
        tmp
    }
    
}



## ===========================================================================
## filterReference
## ---------------------------------------------------------------------------
## References a filter (contained within a filterSet). Everything is just
## passed to the referenced filter. This may be better handled by the type
## system by having "real" filters inherit from concreteFilter (or something)
## and then simply having a setAs(), but I think that will be too much work
## for filter authors.
## ---------------------------------------------------------------------------
setClass("filterReference",
         representation=representation(name="character",
         env="environment"),
         contains="filter")

## Constructor from an environment
setMethod("filterReference",
          signature("environment", "character"),
          function(from, name) {
              new("filterReference", name=name, env=from)
          })

## Constructor from another filterSet
setMethod("filterReference",
          signature("filterSet", "character"),
          function(from,name)
      {
          new("filterReference", env=from@env,
              name=name)
      })



## ===========================================================================
## setOperationFilter
## ---------------------------------------------------------------------------
## Superclass for union intersect, complement and subset filter, which all
## consist of two or more component filters
## ---------------------------------------------------------------------------
setClass("setOperationFilter",
         representation=representation(filters="list"),
         contains="concreteFilter")



## ===========================================================================
## unionFilter 
## ---------------------------------------------------------------------------
## The union of two filters, .i.e, the logical | operation.
## A simple optimization would be to linearize the union of a filter and
## another union filter.
## ---------------------------------------------------------------------------
setClass("unionFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
setMethod("|",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("unionFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "or", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
setMethod("|",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "|", e2=e2))
setMethod("|",
          signature=signature(e1="filter",
          e2="list"),
          definition=function(e1, e2) lapply(e2, "|", e1=e1))



## ===========================================================================
## intersectFilter 
## ---------------------------------------------------------------------------
## The intersection of two filters, i.e, the logical & operation.
## This is somewhat different from the %subset% operation because
## some filters depend on the data and would return different results
## when applied to the full dataset.
## --------------------------------------------------------------------------
setClass("intersectFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
setMethod("&",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("intersectFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "and", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
setMethod("&",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "&", e2=e2))
setMethod("&",
          signature=signature(e1="filter",
          e2="list"),
          definition=function(e1, e2) lapply(e2, "&", e1=e1))



## ===========================================================================
## complementFilter 
## ---------------------------------------------------------------------------
## The complement of a filters, i.e, the logical ! operation.
## ---------------------------------------------------------------------------
setClass("complementFilter",
         representation=representation("setOperationFilter"),
         validity=function(object)
     { 
         if(length(object@filters) != 1) {
             warning("Complement filters can only operate on a ",
                     "single filter")
             return(FALSE)
         }
         TRUE
     })


## constructor
setMethod("!",
          signature=signature(x="filter"),
          definition=function(x)
      {
          new("complementFilter",filters=list(x),
              filterId=paste("not",identifier(x)))
      })



## ===========================================================================
## subsetFilter 
## ---------------------------------------------------------------------------
## Combining two filters in a way that the RHS filter  takes the subset
## of the LHS filter as input. For many cases this is equivalent to an
## intersection filter, the only differnce is in data-driven filters.
## ---------------------------------------------------------------------------
setClass("subsetFilter",
         representation=representation("setOperationFilter"),
         validity=function(object)
     {
         if(length(object@filters) != 2) {
             warning("Subset filters are only defined as binary operators")
             return(FALSE)
         }
         TRUE
     })

## constructor from two filters. %&% is an alias for %subset%
setMethod("%subset%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("subsetFilter",
              filters=list(e1, e2), filterId=paste(identifier(e1),"in",
                                    identifier(e2)))
      })
setMethod("%&%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2) e1 %subset% e2)

## constructor from a list of filters and a filter
setMethod("%subset%",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "%subset%", e2=e2))

## constructor from a filterSet and a filter
setMethod("%subset%",
          signature=signature(e1="filterSet",
          e2="filter"),
          definition=function(e1,e2)
      {
          ## Make a copy of the filterSet, preserving R semantics
          x <- as(as(e1, "list"), "filterSet")
          n <- names(e1)
          x[[""]] <- e2
          target <- as.symbol(identifier(e2))
          for(i in n){
              x[[""]] <- substitute(~ a %subset% b, list(a=as.symbol(i),
                                                         b=target))
          }
          x
      })



## ===========================================================================
## filterResult
## ---------------------------------------------------------------------------
## A container for the results after applying a filter to flow cytometry
## data with slots frameId (identifier of the object) and filterDetails,
## which is a list containing and further describing the input filter.
## ---------------------------------------------------------------------------
setClass("filterResult",
         representation=representation(frameId="character",
         filterDetails="list"),
         contains="concreteFilter",
         prototype=list(frameId="Filter Result",
         filterDetails=list()))



## ===========================================================================
## logicalFilterResult
## ---------------------------------------------------------------------------
## Resuls from a filtering operation that only produces a single population.
## Slot subet is a logical vector indicating the population membership of the
## data in the gated flowFrame.
## ---------------------------------------------------------------------------
setClass("logicalFilterResult",
         representation=representation(subSet="logical"),
         contains="filterResult")



## ===========================================================================
## multipleFilterResult
## ---------------------------------------------------------------------------
## Resuls from a filtering operation that produces multiple populations.
## Slot subet is a factor vector indicating the population membership of the
## data in the gated flowFrame. Factor names are used as population names.
## ---------------------------------------------------------------------------
setClass("multipleFilterResult",
         representation=representation(subSet="factor"),
         contains="filterResult")



## ===========================================================================
## manyFilterResult
## ---------------------------------------------------------------------------
## A special case of multipleFilterResult that arises when there are
## overlapping sets. The subset indices are stored as a matrix, where
## each row contains the results of a single filtering operation.
## ---------------------------------------------------------------------------
setClass("manyFilterResult",
         representation=representation(subSet="matrix",
         dependency="ANY"),
         contains="filterResult")

##constructor
manyFilterResult <- function(filters, frameId, dependency=NULL)
{
    q <- new("manyFilterResult",
             filterDetails=lapply(filters, slot, "filterDetails"),
             subSet=do.call(cbind, lapply(filters, as, "logical")),
             dependency=dependency)
    colnames(q@subSet) <- sapply(filters, slot, "filterId")
    q
}



## ===========================================================================
## randomFilterResult
## ---------------------------------------------------------------------------
## A result of a filtering operation where the population membership is
## considered to be stochastic rather than absolute. Currently there is no
## implementation of a filter that produces such a filterResult, although
## norm2Filter, curvFilters and the t-mixture filters in flowClust are
## obvious candidates.
## ---------------------------------------------------------------------------
setClass("randomFilterResult",
         representation=representation(subSet="numeric"),
         contains="filterResult")



## ===========================================================================
## filterResultList
## ---------------------------------------------------------------------------
## A list of filterResults which typically is generated when applying a
## filter to a whole flowSet. This is a class union of list and filterResult
## and mainly exists to allow for method dispatch and sanity checking.
## FIXME: Do we want to allow for mixed filter classes in the list?
## ---------------------------------------------------------------------------
setClass("filterResultList",
         contains=c("list", "filterResult"))

## Check if a filterResultList matches a flowSet. If strict=TRUE, the
## function will also check whether all items in the filterResultSet
## are of equal type and produce the same number of populations.
validFilterResultList <- function(fres, set, strict=TRUE)
{
    res <- TRUE
    checkClass(fres, "filterResultList")
    checkClass(strict, "logical", 1)
    if(!missing(set)){
        checkClass(set, "flowSet")
        if(res <- !all(names(fres) == sampleNames(set)))
            warning("Sample names don't match between flowSet and ",
                    "filterResultList", call.=FALSE)
    }
    if(strict){
        fTypes <- sapply(fres, function(x) class(x))
        if(length(unique(fTypes)) != 1){
            warning("Not all filterResults in the list are of equal",
                    " type.", call.=FALSE)
            res <- FALSE
        }
        nrPops <- sapply(fres, function(x) length(x))
        if(length(unique(nrPops)) != 1){
            warning("Not all filterResults in the list share the",
                    " same number of sub-populations.", call.=FALSE)
            res <- FALSE
        }
        return(res)
    }
}


## ---------------------------------------------------------------------------
## A list of filters serving as input for a filtering operation of whole
## flowSets or for generating workflow gateActionItems. This directly extends
## class 'list' and mainly exists to allow for method dispatch and sanity
## checking. The filterId slot is supposed to contain a unique identifier for all
## individual filter objects in the list. Names of the list items should
## always correspond to sampleNames of the flowSet.
## ---------------------------------------------------------------------------
setClass("filterList",
         contains="list",
         representation=representation(filterId="character"))

## Check if a filteList matches a flowSet. If strict=TRUE, the
## function will also check whether all items in the filterResultSet
## are of equal type and produce the same number of populations.
validFilterList <- function(flist, set, strict=TRUE)
{
    res <- TRUE
    checkClass(flist, "filterList")
    checkClass(strict, "logical", 1)
    if(!missing(set)){
        checkClass(set, "flowSet")
        if(res <- !all(names(flist) == sampleNames(set)))
            warning("Sample names don't match between flowSet and ",
                    "filterResultList", call.=FALSE)
    }
    if(strict){
        fTypes <- sapply(flist, function(x) class(x))
        if(length(unique(fTypes)) != 1)
        {
            warning("Not all filter objects in the list are of equal",
                    " type.", call.=FALSE)
            res <- FALSE
        }
        if(any(sapply(flist, is, "filterResult")))
        {
            stop("filterResults are not allowed in a filterList") 
            res <- FALSE
        }
        return(res)
    }
}

## Constructor
filterList <- function(x, filterId=identifier(x[[1]]))
{
    checkClass(x, "list")
    checkClass(filterId, "character", 1)
    if(is.null(names(x)))
        stop("Names missing in input list.")
    x <- new("filterList", .Data=x, filterId=filterId)
    validFilterList(x)
    return(x)
}


## ===========================================================================
## filterSummary
## ---------------------------------------------------------------------------
## A class containing the results of calling summary methods on filterResult.
## In the case of multipleFilterrResults, the individual slots(except 'count')
## will be vectors.
## Slots are:
##   - name:  The name of the summary, usually this will be set to be the
##            identifier of the filterResult, or the names of the individual
##            populations for a multipleFilterResult
##   - true:  The number of events in the filter (or the individual
##            populations)
##   - count: The total number of events the filter was applied on
##   - p:     The ratio of events within the filter (i.e., true/count)
## ---------------------------------------------------------------------------
setClass("filterSummary",
         representation=representation(name="character",
         true="numeric",
         count="numeric",
         p="numeric"))



## ===========================================================================
## filterSummaryList
## ---------------------------------------------------------------------------
## A list of filterSummaries which typically is generated when summarizing a
## filterResultList. This directly extends the list class  and mainly exists
## to allow for method dispatch.
## ---------------------------------------------------------------------------
setClass("filterSummaryList",
         contains="list")



## ===========================================================================
## transform functions
## ---------------------------------------------------------------------------
## Constructors for the different varieties of transforms. All of these
## create objects of the basic class 'transform', unless stated otherwise.
## ---------------------------------------------------------------------------
## linear (polynomial) transform constructor
linearTransform <- function(transformationId="defaultLinearTransform",
                            a=1, b=0)
{
    checkClass(a, "numeric")
    checkClass(b, "numeric")
    t <- new("transform", .Data=function(x)  x <- a*x+b)
    t@transformationId <- transformationId
    t
}

## Quadratic transformation constructor
quadraticTransform <- function(transformationId="defaultQuadraticTransform",
                               a=1, b=1, c=0)
{
    if(!is.double(a)) 
        stop("a must be numeric")
    if(!is.double(b))
        stop("b must be numeric")
    if(!is.double(c))
        stop("c must be numeric")
    t <- new("transform", .Data=function(x) x <- a*x^2 + b*x + c)
    t@transformationId <- transformationId
    t
}

## Natural logarithm transformation constructor
lnTransform <- function(transformationId="defaultLnTransform",
                        r=1, d=1)
{
    if(!is.double(r) || r <= 0)
        stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
        stop("d must be numeric")
    t <- new("transform", .Data=function(x)
             x<-log(x)*(r/d))
    t@transformationId <- transformationId
    t
}

## Logarithm transformation constructor
logTransform <- function(transformationId="defaultLogTransform",
                         logbase=10, r=1, d=1)
{
    if(!is.double(r) || r <= 0)
        stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
        stop("d must be numeric")
    if(!is.double(r) || r <=0)
        stop("r must be numeric and positive")
    if(!is.double(logbase) || logbase <= 1)
        stop("logabse must be a pnumeric greater than 1")
    t <- new("transform", .Data=function(x) x <- log(x, logbase)*(r/d))
    t@transformationId <- transformationId
    t
}


## General biexponential transformation constructor
biexponentialTransform <-
    function(transformationId="defaultBiexponentialTransform",
             a=.5, b=1, c=.5, d=1, f=0, w=0,
             tol=.Machine$double.eps^0.25, maxit=as.integer(5000))
{
    t <- new("transform", .Data=function(x)
             x <- biexponential_transform(x, a, b, c, d, f, w, tol, maxit))
    t@transformationId <- transformationId
    t
}

## Logicle transformation constructor
## Input parameters are to be provided in decades
logicleTransform <- function(transformationId="defaultLogicleTransform", 
        w = 0.5, t = 262144, m = 4.5, a = 0) {

    k <- new("transform", .Data=function(x) 
            x <- logicle_transform(as.double(x), as.double(t),as.double(w), as.double(m), as.double(a), FALSE)
            )            
    k@transformationId <- transformationId
    k
}

### Inverse logicle transformation constructor
inverseLogicleTransform <- function(trans, transformationId, ...)UseMethod("inverseLogicleTransform")
inverseLogicleTransform.default <- function(trans, transformationId, ...) {
  
    stop("trans has to be an object of class \"transform\"
            created using the \"logicleTransform\" function\n
         or a 'transformList' created by 'estimateLogicle'\n")
}
inverseLogicleTransform.transform <- function(trans, transformationId, ...) {
    k <- .inverseLogicleTransform(trans@.Data)
   if(missing(transformationId))
    k@transformationId <- paste( "inverse", trans@transformationId, sep ="_")
    k
}
.inverseLogicleTransform <- function(func){
  pars <- c("w", "t", "m", "a")
  vals <- ls(environment(func))
  if(!all(pars %in% vals))
    stop("\"trans\" is not a valid object produced using the
           \"logicle\" function")
  
  w = environment(func)[["w"]] 
  t = environment(func)[["t"]] 
  m = environment(func)[["m"]]
  a = environment(func)[["a"]]
  k <- new("transform", .Data=function(x)
    x <- logicle_transform(as.double(x), as.double(t),as.double(w), as.double(m), as.double(a), TRUE)
  )
  
}
inverseLogicleTransform.transformList <- function(trans, transformationId, ...) {
  invs <- sapply(trans@transforms, function(obj){
    .inverseLogicleTransform(obj@f)
  })
  channels <- names(invs)
  if(missing(transformationId))
    transformationId <- paste( "inverse", trans@transformationId, sep ="_")
  
  transformList(channels, invs, transformationId = transformationId)
}
#' It is mainly trying to estimate w (linearization width in asymptotic decades) value based on given m and data range
#' @param dat flowFrame
#' @param p channel name
#' @param m full length of transformed display in decodes
#' @param t top of the scale of data value
#' @param a additional negative range to be included in display in decades
#' @param q quantile of negative data value (used to adjust w calculation)
#' @param type character either "instrument" or "data". The data range.
#' @noRd
.lgclTrans  <- function(dat, p, t , m, a = 0, q = 0.05, type = "instrument") {
    type <- match.arg(type, c("instrument", "data"))
    transId <- paste(p,"logicleTransform", sep = "_")
    
    rng <- range(dat)
    dat <- exprs(dat)[,p]
    
    if(missing(t)){
      if(type == "instrument")
        t <- rng[,p][2]
      else
        t <- max(dat)
    }
    
    if(missing(m)){
      if(type == "instrument")
        m <- 4.5#hardcoded value to keep consistency with the legacy behavior
      else
        m <- log10(t) + 1 
    }
      
    dat <- dat[dat<0]
    w <- 0
    if(length(dat)) {
        r <- .Machine$double.eps + quantile(dat, q)
        w=(m-log10(t/abs(r))) / 2
        if(w<0)
          stop("w is negative!Try to increase 'm'")
    } 
    logicleTransform( transformationId = transId, w=w, t = t, m = m, a = a)
}

estimateLogicle <- function(x, channels, ...)UseMethod("estimateLogicle")
estimateLogicle.flowFrame <- function(x, channels, ...){
  trans <- .estimateLogicle(x, channels, ...)
  channels <- names(trans)
  transformList(channels, trans)
}
.estimateLogicle <- function(x, channels,...){
            if(!is(x,"flowFrame"))
                stop("x has to be an object of class \"flowFrame\"")
            if(missing(channels))
                stop("Please specify the channels to be logicle transformed");
#            indx <- channels %in% colnames(x)
#            if(!all(indx))
#                stop(paste("Channels", channels[!indx] , "were not found in x ",
#                            sep = " "))
            channels <- sapply(channels, function(channel)getChannelMarker(x, channel)[["name"]], USE.NAMES = FALSE)
            
            sapply(channels, function(p) {
                        .lgclTrans(x, p, ...)               
                    })
              
        }

## Truncation transformation constructor
truncateTransform <- function(transformationId="defaultTruncateTransform",
                              a=1)
{
    t <- new("transform", .Data=function(x){
        x[x<=a] <- a
        x
    })
    t@transformationId <- transformationId
    t
}

## Scale transformation constructor
scaleTransform <- function(transformationId="defaultScaleTransform",
                           a=1, b=10^4)
{
    t <- new("transform", .Data=function(x) (x-a)/(b-a))
    t@transformationId <- transformationId
    t
}

## Split-scale transformation constructor
splitScaleTransform <- function(transformationId="defaultSplitscaleTransform",
                                maxValue=1023,
                                transitionChannel=64, r=192)
{
    maxChannel <- r + transitionChannel
    b <- transitionChannel/2
    d <- 2*log10(exp(1))*r/transitionChannel
    logt <- -2*log10(exp(1))*r/transitionChannel + log10(maxValue)
    t <- 10^logt
    a <- transitionChannel/(2*t)
    logCT <- (a*t+b)*d/r
    c <- 10^logCT/t
    tr <- new("transform", .Data= function(x){
        idx <- which(x <= t)
        idx2 <- which(x > t)
        if(length(idx2)>0)
            x[idx2] <- log10(c*x[idx2])*r/d
        if(length(idx)>0)
            x[idx] <- a*x[idx]+b
        x
    })
    tr@transformationId <- transformationId
    tr
}

## Hyperbolic Arcsin transformation constructor
arcsinhTransform <- function(transformationId="defaultArcsinhTransform",
                             a=1, b=1, c=0)
{
    t <- new("transform", .Data=function(x) asinh(a+b*x)+c)
    t@transformationId <- transformationId
    t
}



## ===========================================================================
## parameterTransform
## ---------------------------------------------------------------------------
## A class used to map parameters of a transform during %on% operations.
## ---------------------------------------------------------------------------
setClass("parameterTransform",
         representation=representation(parameters="character"),
         contains="transform")

## constructor
parameterTransform <- function(FUN, params)
    new("parameterTransform", .Data=as.function(FUN),
        parameters=as.character(params))



## ===========================================================================
## transformMap
## ---------------------------------------------------------------------------
## We want to be able to include transforms within a filter. First we need to
## know which parameters should be input filters
## ---------------------------------------------------------------------------
setClass("transformMap",
         representation=representation(output="character",
         input="character",
         f="function"))



## ===========================================================================
## transformList
## ---------------------------------------------------------------------------
## A list of transformMaps
## ---------------------------------------------------------------------------
setClass("transformList",
         representation=representation(transforms="list",
                                       transformationId="character"),
         prototype=prototype(transformationId="defaultTransformation"),
         validity=function(object)
         if(all(sapply(object@transforms, is, "transformMap"))) TRUE else
         stop("All list items of a 'transformList' must be of class ",
              "'transformMap.'", call.=FALSE))

## constructor
transformList <- function(from, tfun, to=from,
                          transformationId="defaultTransformation")
{
    from <- unique(from)
    to <- unique(to)
    if(!is.character(from) || !is.character(to) || length(from) != length(to))
        stop("'from' and 'to' must be character vectors of equal length.",
             call.=FALSE)
    if(is.character(tfun))
        tfun <- lapply(tfun, get)
    if(!is.list(tfun)) tfun <- list(tfun)
    if(!all(sapply(tfun, is, "function") | sapply(tfun, is, "transform")))
        stop("'tfun' must be a list of functions or a character vector ",
             "with the function names.", call.=FALSE)
    tfun <- rep(tfun, length(from))
    tlist <- mapply(function(x, y, z)
                    new("transformMap", input=x, output=y, 
                    f=if(is(z, "transform")) z@.Data else z),
                    from, to, tfun[1:length(from)])
    tlist <- as(tlist, "transformList")
    identifier(tlist) <- transformationId
    return(tlist)
}



## ===========================================================================
## transformFilter
## ---------------------------------------------------------------------------
## FIXME: I have no clue what that is supposed to be but my guess is that it
## can go away once we have the new transformations in place
## ---------------------------------------------------------------------------
setClass("transformFilter",
         representation=representation(transforms="transformList",
         filter="filter"),
         contains="concreteFilter")



## ===========================================================================
## compensation
## ---------------------------------------------------------------------------
## A class to define a compensation operation.
## Slots are:
##   - compensationId: The identifier of the object
##   - spillover:      The spillover matrix
##   - parameters:     The parameters for which the data is to be compensated,
##                     an object of class parameters 
## ---------------------------------------------------------------------------
setClass("compensation",
         representation(spillover="matrix",
                        compensationId="character",
                        parameters="parameters"),
         prototype=prototype(spillover=matrix(),
                             compensationId="default",
                             parameters=new("parameters",.Data="")))

## Constructor: We allow for the following inputs:
##  spillover is always a symmetric numerical matrix with colnames set
## invert is deprecated
##  invert is always a logical of length 1
##  ..1 is a character vector
##  ..1 is a list of character and/or transformations
##  ..1 is a matrix and spillover is missing
##  ... are characters and/or transformations
## If parameters are given explicitely they need to match the colnames
## of the spillover matrix.
compensation <- function(..., spillover, compensationId="defaultCompensation")
{
    parms <- parseDots(list(...))
    if(missing(spillover))
        spillover <- as.matrix(parms$values)

#    J.Spidlen, Oct 23, 2013: Removed check for square matrices
#    We now support non-square matrices as well
#
#    if(!is.matrix(spillover) || !is.numeric(spillover) ||
#       ncol(spillover) != nrow(spillover))
#        stop("'spillover' must be numeric matrix with same number of ",
#             "rows and columns", call.=FALSE)

    if(!is.matrix(spillover) || !is.numeric(spillover))
        stop("'spillover' must be numeric matrix", call.=FALSE)
    if(is.null(colnames(spillover)))
        stop("Spillover matrix must have colnames", call.=FALSE)
    checkClass(compensationId, "character", 1)
#    checkClass(inv, "logical", 1)
    if(!length(parms$parameters))
        parms <- sapply(colnames(spillover), unitytransform)
    if(all(sapply(parms$parameters,function(x) is(x,"unitytransform"))) &&
       !all(sapply(parms$parameters, parameters) %in% colnames(spillover)))
        stop("Parameters and column names of the spillover matrix ",
             "don't match.", call.=FALSE)
#    if(inv)
      ## spillover <- solve(spillover/max(spillover))
#      spillover <- solve(spillover)
    new("compensation", spillover=spillover, 
        compensationId=compensationId,
        parameters=new("parameters", parms$parameters))
}



## ===========================================================================
## compensatedParameter
## ---------------------------------------------------------------------------
## FIXME NG: Please document
## ---------------------------------------------------------------------------
setClass("compensatedParameter",
          contains=c("transform"),
          representation=representation(parameters="character",spillRefId="character",
                                        searchEnv="environment"
                                       )
        )

## Constructor
compensatedParameter <- function(parameters,
                                 spillRefId="defaultCompensatedParameter",
                                 transformationId="defaultTransformationId",
                                 searchEnv)
{
    
    new("compensatedParameter", parameters=parameters, spillRefId=spillRefId,
        transformationId=transformationId,searchEnv=searchEnv)
}
    	 

## ===========================================================================
## fcReference
## ---------------------------------------------------------------------------
## The general definition of a reference: We have an identifier and an
## evaluation environment to do the lookup in. In the context of workflows,
## the environment will live inside the workFlow object.
## Rather than storing the actual items, the workflow framework will provide
## a more reference based semantic. To this end, we define a number
## of different reference classes to be used in the subsequent class
## definition. Objects of class "fcNullReference" allow to assign
## empty (hence unresolvable) references without breaking dispatch
## or class validity. The task of creating the correct reference type is
## handled by the workFlow-specific assign methods, however these only call
## the appropriate fcReference constructors, so those need to make sure that
## all necesary side-effects take place.
## ---------------------------------------------------------------------------
## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address or similar.
#guid <- function()
#    as.vector(format.hexmode(as.integer(Sys.time())/
#                             runif(1)*proc.time()["elapsed"]))
guid <- function(len=10){
       ltrs <- c(LETTERS,letters)
       paste(c(sample(ltrs,1),sample(c(ltrs,0:9),len-1,replace=TRUE)),collapse="")
}

setClass("fcReference", 
         representation=representation(ID="character",
         env="environment"),
         prototype=prototype(ID=paste("ref", guid(), sep="_"),
         env=new.env(parent=emptyenv()))
         )

## Set a value in the alias table of the workFlow 
setAlias <- function(alias, value, workflow)
{
    checkClass(alias, "character", 1)
    checkClass(value, "character", 1)
    checkClass(workflow, c("workFlow", "environment"))
    aliasEnv <- alias(workflow)
    aliasEnv[[alias]] <- unique(c(aliasEnv[[alias]], value))
    return(invisible(NULL))
}

## Get a value from the alias table of the workFlow 
getAlias <- function(alias, workflow)
{
    checkClass(alias, "character")
    checkClass(workflow, c("workFlow", "environment"))
    fun <- function(x)
        if(x %in% ls(workflow)) x else alias(workflow)[[x]]
    return(as.vector(sapply(alias, fun)))
}

## Figure out which reference type to create, based on the class of 'value'
refType <- function(value)
{
    if(is(value, "flowFrame") || is(value, "flowSet")) "fcDataReference"
    else if(is(value, "filterResult")) "fcFilterResultReference"
    else if(is(value, "filter")) "fcFilterReference"
    else if(is(value, "filterList")) "fcFilterReference"
    else if(is(value, "actionItem")) "fcActionReference"
    else if(is(value, "view")) "fcViewReference"
    else if(is(value, "compensation")) "fcCompensateReference"
    else if(is(value, "transformList")) "fcTransformReference"
    else if(is(value, "graphNEL")) "fcTreeReference"
    else if(is(value, "environment")) "fcAliasReference"
    else if(is(value, "normalization")) "fcNormalizationReference"
    else if(is(value, "subsetting")) "fcSubsettingReference"
    else if(is(value, "list")) "fcJournalReference"
    else if(is.null(value)) "fcNullReference"
    else "fcReference"
}

## Create useful identifiers for references
refName <- function(value)
{
    prefix <- if(is(value, "flowFrame") || is(value, "flowSet")) "dataRef"
    else if(is(value, "filterResult")) "fresRef"
    else if(is(value, "filter")) "filterRef"
    else if(is(value, "filterList")) "filterRef"
    else if(is(value, "actionItem")) "actionRef"
    else if(is(value, "view")) "viewRef"
    else if(is(value, "compensation")) "compRef"
    else if(is(value, "transformList")) "transRef"
    else if(is(value, "graphNEL")) "treeRef"
    else if(is(value, "environment")) "aliasRef"
    else if(is(value, "normalization")) "normRef"
    else if(is(value, "subsetting")) "subRef"
    else if(is(value, "list")) "journalRef"
    else if(is.null(value)) "nullRef"
    else "genericRef"
    return(paste(prefix, guid(), sep="_"))
}

## constructor
fcReference <- function(ID=paste("genericRef", guid(), sep="_"),
                        env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workflow")) env@env else env
    ref <- new("fcReference", ID=ID, env=envir)
    setAlias(substitute(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcStructureReference
## ---------------------------------------------------------------------------
## This only exists to subclass structural parts of the workFlow like views
## and actionItems to allow for group-wise method dispatch.
## ---------------------------------------------------------------------------
setClass("fcStructureReference",
         contains=list("VIRTUAL",
         "fcReference"))



## ===========================================================================
## fcTreeReference
## ---------------------------------------------------------------------------
## A reference to a graphNEL object representing the workflow tree
## ---------------------------------------------------------------------------
setClass("fcTreeReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("treeRef", guid(), sep="_"))
         )

## constructor
fcTreeReference <- function(ID=paste("treeRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcTreeReference", ID=ID, env=envir)
    setAlias(".__tree", identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcJournalReference
## ---------------------------------------------------------------------------
## A reference to a list representing the journal
## ---------------------------------------------------------------------------
setClass("fcJournalReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("journalRef", guid(), sep="_"))
         )

## constructor
fcJournalReference <- function(ID=paste("journalRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcJournalReference", ID=ID, env=envir)
    setAlias(".__journal", identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcAliasReference
## ---------------------------------------------------------------------------
## A reference to an environment object representing the alias table
## ---------------------------------------------------------------------------
setClass("fcAliasReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("aliasRef", guid(), sep="_"))
         )

## constructor
fcAliasReference <- function(ID=paste("aliasRef", guid(), sep="_"),
                             env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    new("fcAliasReference", ID=ID, env=envir)
}



## ===========================================================================
## fcDataReference
## ---------------------------------------------------------------------------
## A reference to a data type object (a flowSet or flowFrame)
## ---------------------------------------------------------------------------
setClass("fcDataReference",
         contains="fcReference",
         prototype=prototype(ID=paste("dataRef", guid(), sep="_"))
         )

## constructor
fcDataReference <- function(ID=paste("dataRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcDataReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcActionReference
## ---------------------------------------------------------------------------
## A reference to an action item object (a gate, transformation or
## compensation)
## ---------------------------------------------------------------------------
setClass("fcActionReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("actionRef", guid(), sep="_"))
         )

## constructor
fcActionReference <- function(ID=paste("actionRef", guid(), sep="_"),
                              env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcActionReference", ID=ID, env=envir)
    setAlias(names(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcViewReference
## ---------------------------------------------------------------------------
## A reference to a view object. This allows to bind action items to
## particular views on the data
## ---------------------------------------------------------------------------
setClass("fcViewReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("viewRef", guid(), sep="_"))
         )

## constructor
fcViewReference <- function(ID=paste("viewRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcViewReference", ID=ID, env=envir)
    setAlias(names(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcFilterResultReference
## ---------------------------------------------------------------------------
## A reference to a filterResult object. We need this to store filterResults
## along with the views the filtering generated, without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcFilterResultReference",
         contains="fcReference",
         prototype=prototype(ID=paste("fresRef", guid(), sep="_"))
         )

## constructor
fcFilterResultReference <- function(ID=paste("fresRef",
                                    guid(), sep="_"),
                                    env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcFilterResultReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcFilterReference
## ---------------------------------------------------------------------------
## A reference to a filter object. We need this to store filters within 
## a gateActionItem without unnecessarily copying things.
## ---------------------------------------------------------------------------
setClass("fcFilterReference",
         contains="fcReference",
         prototype=prototype(ID=paste("filterRef", guid(), sep="_"))
         )

## constructor
fcFilterReference <- function(ID=paste("filterRef",
                              guid(), sep="_"),
                              env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcFilterReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcCompensateReference
## ---------------------------------------------------------------------------
## A reference to a compensation object. We need this to store a compensation
## within a compensateActionItem without unnecessarily copying things.
## ---------------------------------------------------------------------------
setClass("fcCompensateReference",
         contains="fcReference",
         prototype=prototype(ID=paste("compRef", guid(), sep="_"))
         )

## constructor
fcCompensateReference <- function(ID=paste("compRef",
                                  guid(), sep="_"),
                                  env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcCompensateReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcTransformReference
## ---------------------------------------------------------------------------
## A reference to a transformation object. For now we use transformList
## objects until transformation is a more usefull class. We need this to
## store a trasnforamtion within a transformActionItem without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcTransformReference",
         contains="fcReference",
         prototype=prototype(ID=paste("transRef", guid(), sep="_"))
         )

## constructor
fcTransformReference <- function(ID=paste("transRef",
                                 guid(), sep="_"),
                                 env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcTransformReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcNormalizationReference
## ---------------------------------------------------------------------------
## A reference to a normalization object. We need this to store a
## normalization within a normActionItem without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcNormalizationReference",
         contains="fcReference",
         prototype=prototype(ID=paste("normRef", guid(), sep="_"))
         )

## constructor
fcNormalizationReference <- function(ID=paste("normRef",
                                              guid(), sep="_"),
                                     env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcNormalizationReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcSubsettingReference
## ---------------------------------------------------------------------------
## A reference to a subsetting object. We need this to store a
## subsetting within a normActionItem without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcSubsettingReference",
         contains="fcReference",
         prototype=prototype(ID=paste("subRef", guid(), sep="_"))
         )

## constructor
fcSubsettingReference <- function(ID=paste("subRef",
                                              guid(), sep="_"),
                                     env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, c("workFlow", "environment"))
    envir <- if(is(env, "workFlow")) env@env else env
    ref <- new("fcSubsettingReference", ID=ID, env=envir)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcNullReference
## ---------------------------------------------------------------------------
## A NULL reference to be used whenever a slot is supposed to be empty
## ---------------------------------------------------------------------------
setClass("fcNullReference",
         contains=c("fcDataReference",
                    "fcActionReference",
                    "fcViewReference",
                    "fcFilterResultReference",
                    "fcFilterReference",
                    "fcCompensateReference",
                    "fcTransformReference",
                    "fcNormalizationReference",
                    "fcSubsettingReference", 
                    "fcTreeReference",
                    "fcJournalReference",
                    "fcAliasReference"),
         prototype=prototype(ID=paste("nullRef", guid(), sep="_"))
         )

fcNullReference <- function(...) new("fcNullReference")



## ===========================================================================
## normalization
## ---------------------------------------------------------------------------
## A class to describe normalization operations on a complete flowSet.
## Currently this is only the warping, but more methods may follow. The
## class mainly exists to allow for dispatch in the workFlow system. The
## function 'normFunction' is supposed to take a flowSet, perform an
## operation on 'parameters' and return the altered flowSet. It has two
## mandatory arguments: 'x' and 'parameters'. All additional arguments
## have to be supplied via the list in the 'arguments' slot.
## ---------------------------------------------------------------------------
setClass("normalization",
         representation(parameters="character",
                        normalizationId="character",
                        normFunction="function",
                        arguments="list"),
         prototype=prototype(normalizationId="defaultNormalization",
                             normFunction=function(x) x)
         )

## constructor
normalization <- function(parameters, normalizationId="defaultNormalization",
                          normFunction, arguments=list())
{
    checkClass(normalizationId, "character", 1)
    checkClass(parameters, "character")
    checkClass(normFunction, "function")
    new("normalization", parameters=parameters,
        normalizationId=normalizationId, normFunction=normFunction,
        arguments=arguments)
}




## ===========================================================================
## subsetting
## ---------------------------------------------------------------------------
## A class to describe subsetting operations on a flowSet. The
## class mainly exists to allow for dispatch in the workFlow system.
## ---------------------------------------------------------------------------
setClassUnion("characterOrNumeric", c("character","numeric"))
setClass("subsetting",
         representation(subsettingId="character",
                        indices="characterOrNumeric"),
         prototype=prototype(subsettingId="defaultSubsetting")
         )

## constructor
subsetting <- function(indices, subsettingId="defaultSubsetting")
{
    checkClass(indices, c("character","numeric"))
    checkClass(subsettingId, "character",1)
    new("subsetting", indices=indices,
        subsettingId=subsettingId)
}







## ===========================================================================
## workFlow
## ---------------------------------------------------------------------------
## This class is intended to store all necessary information about a
## workflow. Since the individual bits and pieces (actionItems, views) know
## about their inter-relations, it is only the tree that links the
## individual views and the common evaluation environment containing all
## stored objects. Nodes in the tree correspond to views, edges to
## actionItems, stored as references in the edgeData slot. The tree itself
## is also stored in the environment, which allows for reference-based
## updating without the necessaty of an assignment method or the like.
## In addition to the tree we store an alias table in the environment.
## Internally, all objects are referenced by their guid, but we allow for
## more human readable aliases (usually the "name" slot) which can be used
## to identify objects if they are unique. Whenever possible we try to plot
## these readable names and those are also available for completion.
## Note that the environment in the prototype gets created once and all
## objects created via "new" without explicitely defining "env" will
## essentially share a common environment. This is fixed in the constructor.
## We also keep a journal in the 'journal' slot, which essentially is a list
## holding references to all objects that are created by any operation on the
## workflow. This allows for a simple stepwise undo mechanism and also for
## some clean-up by the error handling in case an operation fails.
## ---------------------------------------------------------------------------
setClass("workFlow",
         representation=representation(name="character",
         tree="fcTreeReference",
         alias="fcAliasReference",
         journal="fcJournalReference",
         env="environment"),
         prototype=prototype(name="default",
         tree=fcNullReference(),
         alias=fcNullReference(),
         journal=fcNullReference(),
         env=new.env(parent=emptyenv())))


## make deep copy of a flowSet
copyFlowSet <- function(x) x[1:length(x)]

## copy a flowFrame
copyFlowFrame <- function(x) x[1:nrow(x)]

## The constructor takes a flow data object (flowFrame or flowSet) and
## makes a copy in the evaluation environment. It also sets up the views
## graph and the alias table in the environment
workFlow <- function(data, name="default", env=new.env(parent=emptyenv()))
{
  .Deprecated("flowWorkspace::GatingSet")
    ## some sanity checks up front
    checkClass(name, "character", 1)
    checkClass(env, "environment")
    wf <-  new("workFlow", name=name, env=env)
    ## set up the alias table as an environment in the workFlow
    aliasTable <- new.env(hash=TRUE, parent=emptyenv())
    id <- refName(aliasTable)
    assign(".__alias", id, aliasTable)
    assign(id, aliasTable, wf@env)
    wf@alias <- new("fcAliasReference", env=wf@env, ID=id)
    ## Set up the views tree and the journal
    tree <- new("graphNEL", edgemode="directed")
    journal <- list()
    ## Assign the data to the workFlow and create a base view
    if(!missing(data)){
        if(!is(data, "flowFrame") && !is(data, "flowSet"))
            stop("'data' must be a flow data structure (flowFrame or flowSet)",
                 call.=FALSE)
	if(is(data,"flowSet"))
		data <- copyFlowSet(data)
	if(is(data,"flowFrame"))
	    data <- copyFlowFrame(data)
	
        dataRef <- assign(value=data, envir=wf)
        viewRef <- view(workflow=wf, name="base view", data=dataRef)
        tree <- addNode(identifier(viewRef), tree)
        journal <- list(c(identifier(viewRef), identifier(dataRef)))
        names(journal) <- ".action_baseView"
    }
    wf@tree <- assign(value=tree, envir=wf)
    wf@journal <- assign(value=journal, envir=wf)
    return(wf)
}


## Create a reference by assigning 'value' to the symbol 'x' in the
## evaluation environment in 'workflow'. If 'x' is not specified, we
## create a reasonable unique identifier. These methods return a reference
## to 'value' for further use.
## Note that creation of a NULL reference does not result in any assignment
## to the environment, but still a fcNullReference object is returned.


## Assign any object to a workFlow object, the symbol (as guid) is created
## automatically
setMethod("assign",
          signature=signature(x="missing",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(value, pos)
      {
          id <- refName(value)
          if(!is.null(value))
              assign(id, value, envir=pos)
          a <- do.call(refType(value), list(ID=id, env=pos))
          return(a)
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("assign",
          signature=signature(x="missing",
          value="ANY",
          pos="missing",
          envir="workFlow",
          inherits="missing",
          immediate="missing"),
          definition=function(value, envir) assign(value=value, pos=envir))

## Assign to a particular symbol (potentially overwriting existing ones)
setMethod("assign",
          signature=signature(x="character",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, pos)
      {
          rmAlias(x, pos)
          if(!is.null(value)){
              if(x %in% ls(pos))
                  warning("Overwriting object in the environment.", call.=FALSE)
              ## add new item to journal
              jid <- identifier(pos@journal)
              journal <- get(pos@journal)
              lj <- length(journal)
              assign(x, value, envir=pos@env)
              if(lj){
                  journal[[length(journal)]] <- c(journal[[length(journal)]], x)
                  assign(jid, journal, pos@env)
              }
          }else{
              suppressWarnings(rm(list=x, envir=pos@env))
          }
          do.call(refType(value), list(ID=x, env=pos))     
      })

## Assign via existing reference.
setMethod("assign",
          signature=signature(x="fcReference",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, pos)
      {
          rmAlias(identifier(x), pos)
          if(is.null(value)){
              Rm(x, rmRef=FALSE)
          }else{
              assign(identifier(x), value, envir=pos@env)  
          }
          do.call(refType(value), list(ID=identifier(x), env=pos))
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("assign",
          signature=signature(x="ANY",
          value="ANY",
          pos="missing",
          envir="workFlow",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, envir)
      {
          assign(x=x, value=value, pos=envir)
      })



## ===========================================================================
## actionItem
## ---------------------------------------------------------------------------
## actionItems are either gates, transformations or compensations that
## work on a particular view of the data, linked to it by the parentView
## slot.
## ---------------------------------------------------------------------------
setClass("actionItem", 
         representation=representation("VIRTUAL",
         ID="character",
         name="character",
         parentView="fcViewReference",
         alias="fcAliasReference",
         env="environment"),
         prototype=prototype(ID=paste("actionRef", guid(), sep="_"),
         name="",
         parentView=fcNullReference(),
         alias=fcNullReference(),
         env=new.env(parent=emptyenv())))



## ===========================================================================
## gateActionItem
## ---------------------------------------------------------------------------
## A subclass od actionItem. This contains the definiton of the filter/gate
## operation
## ---------------------------------------------------------------------------
setClass("gateActionItem",
         contains="actionItem",
         representation=representation(gate="fcFilterReference",
         filterResult="fcFilterResultReference"),
         prototype=prototype(ID=paste("gateActionRef", guid(), sep="_"),
         gate=fcNullReference(),
         filterResult=fcNullReference()))

## The constructor creates the gateActionItem object and directly assigns
## it to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
gateActionItem <- function(ID=paste("gateActionRef", guid(), sep="_"),
                           name=paste("action", identifier(get(gate)), sep="_"),
                           parentView, gate, filterResult, workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(gate, "fcFilterReference")
    checkClass(parentView, "fcViewReference")
    if(missing(filterResult))
        filterResult <-fcNullReference()
    action <- new("gateActionItem", ID=ID, name=name, gate=gate,
                  parentView=parentView, env=workflow@env,
                  filterResult=filterResult, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}



## ===========================================================================
## transformActionItem
## ---------------------------------------------------------------------------
## transformation actionItem. This contains the definiton of the
## transformation, which is not defined yet. For now, we use the
## transformMapList class for this purpose
## ---------------------------------------------------------------------------
setClass("transformActionItem",
         contains="actionItem",
         representation=representation(transform="fcTransformReference"))

## The constructor creates the transformActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
transformActionItem <- function(ID=paste("transActionRef", guid(), sep="_"),
                                name=paste("action", identifier(get(transform)),
                                           sep="_"),
                                parentView, transform, workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(transform, "fcTransformReference")
    checkClass(parentView, "fcViewReference")
    action <- new("transformActionItem", ID=ID, name=name,
                  transform=transform, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}


## ===========================================================================
## compensateActionItem
## ---------------------------------------------------------------------------
## compensation actionItem. This contains the definiton of the
## compensation
## ---------------------------------------------------------------------------
setClass("compensateActionItem",
         contains="actionItem",
         representation=representation(compensate="fcCompensateReference"))

## The constructor creates the compensateActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
compensateActionItem <- function(ID=paste("compActionRef", guid(), sep="_"),
                                 name=paste("action", identifier(get(compensate)),
                                            sep="_"),
                                 parentView, compensate, workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(compensate, "fcCompensateReference")
    checkClass(parentView, "fcViewReference")
    action <- new("compensateActionItem", ID=ID, name=name,
                  compensate=compensate, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}


## ===========================================================================
## normalizeActionItem
## ---------------------------------------------------------------------------
## normalization actionItem. This contains the definiton of the
## normalization
## ---------------------------------------------------------------------------
setClass("normalizeActionItem",
         contains="actionItem",
         representation=representation(normalization="fcNormalizationReference"))

## The constructor creates the normalizeActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
normalizeActionItem <- function(ID=paste("normActionRef", guid(), sep="_"),
                                name=paste("action", identifier(get(normalization)),
                                           sep="_"),
                                 parentView, normalization,
                                 workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(normalization, "fcNormalizationReference")
    checkClass(parentView, "fcViewReference")
    action <- new("normalizeActionItem", ID=ID, name=name,
                  normalization=normalization, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}


## ===========================================================================
## subsettingActionItem
## ---------------------------------------------------------------------------
## This contains the definiton of the subsetting
## ---------------------------------------------------------------------------
setClass("subsettingActionItem",
         contains="actionItem",
         representation=representation(subsetting="fcSubsettingReference"))

## The constructor creates the subsettingActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
subsettingActionItem <- function(ID=paste("subActionRef", guid(), sep="_"),
                                name=paste("action", identifier(get(subsetting)),
                                           sep="_"),
                                 parentView, subsetting,
                                 workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(subsetting, "fcSubsettingReference")
    checkClass(parentView, "fcViewReference")
    action <- new("subsettingActionItem", ID=ID, name=name,
                  subsetting=subsetting, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}



## ===========================================================================
## view
## ---------------------------------------------------------------------------
## The concept of views in this context does not adhere strictly to the
## definition. Whenever possible, we try to provide a "real" view on the
## original (or at least on the parent) data, however operations like
## compensation or transformation will alter the data values, and unless
## we want to add new columns to the original data (or make a copy while
## retaining the original set), there is no simple way around that. Views
## mainly exist in order to provide a uniform organisational structure
## for the data after an actionItem has been applied. The first node in
## the workflow tree is always a view with a NULL action reference, and
## a non-NULL data reference.
## ---------------------------------------------------------------------------
setClass("view", 
         representation=representation(ID="character",
         name="character",
         action="fcActionReference",
         env="environment",
         alias="fcAliasReference",
         data="fcDataReference"),
         prototype=prototype(ID=paste("view", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         alias=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the view object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
view <- function(workflow, ID=paste("viewRef", guid(), sep="_"),
                 name="default", data, action)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    if(missing(data))
        data <- fcNullReference()
    if(missing(action))
        action <- fcNullReference()
    bv <-  new("view", ID=ID, name=name, env=workflow@env,
               action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}



## ===========================================================================
## gateView
## ---------------------------------------------------------------------------
## Gate views store indices of the gating result for further subsetting as a
## list, where each list item contains indices for a single flowFrame. No
## subset will be produced unless another actionItem is applied to the view.
## We need some form of special treatment for gates that produce multiple
## populations. A gateView will always capture the result for only a single
## sub-population, however, the whole filterResult is necessary for plotting
## and a reference to that and the index in the filterResult (i.e., the
## subpopulation) will be stored along with the view . 
## ---------------------------------------------------------------------------
setClass("gateView",
         contains="view",
         representation=representation(indices="list",
         filterResult="fcFilterResultReference",
         frEntry="character"),
         prototype=prototype(ID=paste("gateViewRef", guid(), sep="_"),
         name="",
         filterResult=fcNullReference(),
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the gateView object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
gateView <- function(workflow, ID=paste("gateViewRef", guid(), sep="_"),
                     name="default", action, data, indices, 
                     filterResult, frEntry)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(frEntry, "character", 1)
    checkClass(indices, "list")
    checkClass(action, "fcActionReference")
    checkClass(filterResult, "fcFilterResultReference")
    if(missing(data))
        data <- fcNullReference()
    bv <- new("gateView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, indices=indices,
              filterResult=filterResult, frEntry=frEntry,
              alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}


## Helper function to create gateActionItem from either a filter or filterList
## object
createGateActionItem <- function(wf, action, parent, name)
{
    if(is(action, "filterResult"))
        stop("Don't know how to handle object of class '",
             class(action), "'", call.=FALSE)
    else if(is(action, "filter") || is(action, "filterList")){
        ## get the parentView. If not explicitely specified, use the
        ## root node
        pid <- if(is.null(parent)) views(wf)[[1]] else parent
        pid <- getAlias(pid, wf)
        if(is.null(unlist(pid)))
            stop("'", parent, "' is not a valid view name in this",
                 " workflow.", call.=FALSE)
        ## assign the filter to the evaluation environment and create
        ## a reference to it
        gateRef <- assign(value=action, envir=wf)
        pview <- fcViewReference(ID=pid, env=wf)
        ## create and assign a new gateActionItem
        actionRef <- gateActionItem(parentView=pview, gate=gateRef,
                                    workflow=wf)
        ## add a useful name to the journal entry
        journal <- get(wf@journal)
        names(journal)[length(journal)] <- identifier(actionRef)
        assign(wf@journal, journal, wf)
        ## now evaluate the filter and assign the result
        fres <- filter(Data(get(pview)), action)
        if(is(action, "filterList"))
            identifier(fres) <- identifier(action)
		if(is(fres,"filterResultList")){
        	if(!validFilterResultList(fres))
            	stop("Don't know how to proceed.", call.=FALSE)
		}
        fresRef <-  assign(value=fres, envir=wf)
        gAction <- get(actionRef)
        gAction@filterResult <- fresRef
        assign(actionRef, gAction, wf) 
        ## we need to distinguish between logicalFilterResults and
        ## multipleFilterResults
        nodes <- NULL
        warned <- FALSE
        if(!is(fres, "filterResultList")){
            len <-
                if(is(fres, "logicalFilterResult")) 2
                else length(fres)
            for(i in seq_len(len)){
                vid <- gateView(name=paste(name, names(fres)[i], sep=""),
                                workflow=wf,
                                action=actionRef, filterResult=fresRef,
                                indices=list(fres[[i]]@subSet),
                                frEntry=names(fres)[i])
                if(!warned)
                {
                    isUniqueName(names(get(vid)), wf, TRUE, identifier(vid))
                    warned <- TRUE
                }
                nodes <- c(nodes, identifier(vid))
            }
        }else{
            len <-
                if(is(fres[[1]], "logicalFilterResult")) 2
                else length(fres[[1]])
            for(i in seq_len(len)){
                vid <- gateView(name=paste(name, names(fres[[1]])[i], sep=""),
                                workflow=wf,
                                action=actionRef, filterResult=fresRef,
                                indices=lapply(fres, function(y) y[[i]]@subSet),
                                frEntry=names(fres[[1]])[i])
                if(!warned)
                {
                    isUniqueName(names(get(vid)), wf, TRUE, identifier(vid))
                    warned <- TRUE
                }
                nodes <- c(nodes, identifier(vid))
            }
        }
        ## update the filter and filterResult IDs
        identifier(action) <- paste("filter", identifier(action),
                                    sep="_")
        identifier(fres) <- paste("fres", identifier(fres),
                                  sep="_")
        assign(gateRef, action, wf)
        assign(fresRef, fres, wf)
        ## add new nodes and edges to the workflow tree
        tree <- get(wf@tree)
        tree <- addNode(nodes, tree)
        tree <- addEdge(pview@ID, nodes, tree)
        edgeDataDefaults(tree, "actionItem") <- fcNullReference()
        edgeData(tree, pview@ID, nodes, "actionItem") <-
            actionRef
        assign(x=wf@tree, value=tree, envir=wf)
        return(invisible(wf))
    } else stop("Don't know how to handle object of class '",
                class(action), "'", call.=FALSE)
}

## constructor directly from a filter object. This creates a gateActionItem
## in the workFlow and from that a gateView which is also directly stored in
## the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="concreteFilter"),
          definition=function(wf, action, parent=NULL, name="")
      {
         
          ## create a new journal entry first
          journal <- c(get(wf@journal), list(NULL))
          assign(wf@journal, journal, wf)
          ## now do the actual adding, new objects are captured in the journal
          tryCatch(createGateActionItem(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        message("Error adding item '", substitute(action),
                            "':\n  ", conditionMessage(e))
                        undo(wf)
                    },
                   warning=function(e) message("Warning: ", conditionMessage(e)))
          return(invisible(wf))
      })



## constructor directly from a filterList object. This creates a gateActionItem
## in the workFlow and from that a gateView which is also directly stored in
## the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="filterList"),
          definition=function(wf, action, parent=NULL, name="")
      {
         
          ## create a new journal entry first
          journal <- c(get(wf@journal), list(NULL))
          assign(wf@journal, journal, wf)
          ## now do the actual adding, new objects are captured in the journal
          tryCatch(createGateActionItem(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        cat("Error adding item '", substitute(action),
                            "':\n  ", sep="")
                        cat(gsub("Error.*\\)\\:", "", e))
                        undo(wf)
                    })
          return(invisible(wf))
      })



## ===========================================================================
## transformView
## ---------------------------------------------------------------------------
## A transfomation makes a copy of the data independent of whether it
## is a leaf node or not.
## Do we want to allow to introduce new data columns. My guess is that we
## have to, but for now let's assume we don't.
## ---------------------------------------------------------------------------
setClass("transformView",
         contains="view",
         prototype=prototype(ID=paste("transViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the transformView object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
transformView <- function(workflow, ID=paste("transViewRef", guid(), sep="_"),
                          name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("transformView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a transformation object. This creates a
## transformActionItem in the workFlow and from that a transformView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="transformList"),
          definition=function(wf, action, parent=NULL,
                              name=identifier(action))
      {
          fun <- function(wf, action, parent, name)
          {
              ## get the parentView. If not explicitely specified, use the
              ## root node
              pid <- if(is.null(parent)) views(wf)[[1]] else parent
              pid <- getAlias(pid, wf)
              if(is.null(unlist(pid)))
                  stop("'", parent, "' is not a valid view name in this",
                       " workflow.", call.=FALSE)
              ## assign the transformation to the evaluation environment and
              ## create a reference to it
              transRef <- assign(value=action, envir=wf)
              pview <- fcViewReference(ID=pid, env=wf)
              tree <- get(wf@tree)
              if(length(unlist(adj(tree, pid))))
                  warning("The selected parent view is not a leaf node.\n",
                          "Don't know how to update yet.", call.=FALSE)
              ## create and assign a new transformActionItem
              actionRef <- transformActionItem(parentView=pview,
                                               transform=transRef,
                                               workflow=wf)
              ## add a useful name to the journal entry
              journal <- get(wf@journal)
              names(journal)[length(journal)] <- identifier(actionRef)
              assign(wf@journal, journal, wf)
              ## now transform the data and assign the result
              tData <- action %on% Data(get(pview))
              dataRef <- assign(value=tData, envir=wf)
              vid <- transformView(name=name, workflow=wf,
                                   action=actionRef, data=dataRef)
              ## update the identifier of the transformation object
              identifier(action) <- paste("trans", identifier(action),
                                          sep="_")
              assign(transRef, value=action, envir=wf) 
              ## add new nodes and edges to the workflow tree
              nid <- identifier(vid)
              tree <- addNode(nid, tree)
              tree <- addEdge(pid, identifier(vid), tree)
              edgeDataDefaults(tree, "actionItem") <- fcNullReference()
              edgeData(tree, pid , nid, "actionItem") <- actionRef
              assign(x=wf@tree, value=tree, envir=wf)
              return(invisible(wf))
          }
          ## create a new journal entry first
          journal <- c(get(wf@journal), list(NULL))
          assign(wf@journal, journal, wf)
          ## now do the actual adding, new objects are captured in the journal
          tryCatch(fun(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        cat("Error adding item '", substitute(action),
                            "':\n  ", sep="")
                        cat(gsub("Error.*\\)\\:", "", e))
                        undo(wf)
                    })
          return(invisible(wf))
      })



## ===========================================================================
## compensateView
## ---------------------------------------------------------------------------
## A compensation makes a copy of the data independent of whether it
## is a leaf node or not.
## ---------------------------------------------------------------------------
setClass("compensateView",
         contains="view",
         prototype=prototype(ID=paste("compViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the compensateView object and directly assigns it
## to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
compensateView <- function(workflow, ID=paste("compViewRef", guid(), sep="_"),
                           name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("compensateView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a compensation object. This creates a
## compensateActionItem in the workFlow and from that a compensateView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="compensation"),
          definition=function(wf, action, parent=NULL,
          name=identifier(action))
      {
          fun <- function(wf, action, parent, name)
          {
              ## get the parentView. If not explicitely specified, use the
              ## root node
              pid <- if(is.null(parent)) views(wf)[[1]] else parent
              pid <- getAlias(pid, wf)
              if(is.null(unlist(pid)))
                  stop("'", parent, "' is not a valid view name in this",
                       " workflow.", call.=FALSE)
              if(pid != getAlias(views(wf), wf)[1])
                  warning("The selected parent view is not a root node.\n",
                      "Are you sure this is correct?", call.=FALSE)
              ## assign the compensation to the evaluation environment and
              ## create a reference to it
              compRef <- assign(value=action, envir=wf)
              pview <- fcViewReference(ID=pid, env=wf)
              ## create and assign a new ActionItem
              actionRef <- compensateActionItem(parentView=pview,
                                                compensate=compRef,
                                                workflow=wf)
              ## add a useful name to the journal entry
              journal <- get(wf@journal)
              names(journal)[length(journal)] <- identifier(actionRef)
              assign(wf@journal, journal, wf)
              ## now transform the data and assign the result
              tData <- compensate(Data(get(pview)), action)
              dataRef <- assign(value=tData, envir=wf)
              vid <- compensateView(name=name, workflow=wf,
                                    action=actionRef, data=dataRef)
              ## update the identifier of the compensation object
              identifier(action) <- paste("comp", identifier(action),
                                          sep="_")
              assign(compRef, value=action, envir=wf) 
              ## add new nodes and edges to the workflow tree
              nid <- identifier(vid)
              tree <- get(wf@tree)
              tree <- addNode(nid, tree)
              tree <- addEdge(pid, identifier(vid), tree)
              edgeDataDefaults(tree, "actionItem") <- fcNullReference()
              edgeData(tree, pid , nid, "actionItem") <- actionRef
              assign(x=wf@tree, value=tree, envir=wf)
              return(invisible(wf))
          }
          ## create a new journal entry first
          journal <- c(get(wf@journal), list(NULL))
          assign(wf@journal, journal, wf)
          ## now do the actual adding, new objects are captured in the journal
          tryCatch(fun(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        cat("Error adding item '", substitute(action),
                            "':\n  ", sep="")
                        cat(gsub("Error.*\\)\\:", "", e))
                        undo(wf)
                    })
          return(invisible(wf))
      })



## ===========================================================================
## normalizeView
## ---------------------------------------------------------------------------
## A normalization makes a copy of the data independent of whether it
## is a leaf node or not.
## ---------------------------------------------------------------------------
setClass("normalizeView",
         contains="view",
         prototype=prototype(ID=paste("normViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the normalizeView object and directly assigns it
## to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
normalizeView <- function(workflow, ID=paste("normViewRef", guid(), sep="_"),
                           name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("normalizeView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a normalization object. This creates a
## normalizeActionItem in the workFlow and from that a normalizeView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="normalization"),
          definition=function(wf, action, parent=NULL,
                              name=identifier(action))
      {
           fun <- function(wf, action, parent, name)
           {
               ## get the parentView. If not explicitely specified, use the
               ## root node
               pid <- if(is.null(parent)) views(wf)[[1]] else parent
               pid <- getAlias(pid, wf)
               if(is.null(unlist(pid)))
                   stop("'", parent, "' is not a valid view name in this",
                        " workflow.", call.=FALSE)
               ## assign the normalization to the evaluation environment and
               ## create a reference to it
               normRef <- assign(value=action, envir=wf)
               pview <- fcViewReference(ID=pid, env=wf)
               ## create and assign a new ActionItem
               actionRef <- normalizeActionItem(parentView=pview,
                                                normalization=normRef,
                                                workflow=wf)
               ## add a useful name to the journal entry
               journal <- get(wf@journal)
               names(journal)[length(journal)] <- identifier(actionRef)
               assign(wf@journal, journal, wf)
               ## now normalize the data and assign the result
               tData <- normalize(Data(get(pview)), action)
               dataRef <- assign(value=tData, envir=wf)
               vid <- normalizeView(name=name, workflow=wf,
                                    action=actionRef, data=dataRef)
               ## update the identifier of the normalization object
               identifier(action) <- paste("norm", identifier(action),
                                           sep="_")
               assign(normRef, value=action, envir=wf) 
               ## add new nodes and edges to the workflow tree
               nid <- identifier(vid)
               tree <- get(wf@tree)
               tree <- addNode(nid, tree)
               tree <- addEdge(pid, identifier(vid), tree)
               edgeDataDefaults(tree, "actionItem") <- fcNullReference()
               edgeData(tree, pid , nid, "actionItem") <- actionRef
               assign(x=wf@tree, value=tree, envir=wf)
               return(invisible(wf))
           }
           ## create a new journal entry first
           journal <- c(get(wf@journal), list(NULL))
           assign(wf@journal, journal, wf)
           ## now do the actual adding, new objects are captured in the journal
           tryCatch(fun(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        cat("Error adding item '", substitute(action),
                            "':\n  ", sep="")
                        cat(gsub("Error.*\\)\\:", "", e))
                        undo(wf)
                    })
           return(invisible(wf))
       })



## ===========================================================================
## subsettingView
## ---------------------------------------------------------------------------
## A subsetting makes a copy of the data independent of whether it
## is a leaf node or not.
## ---------------------------------------------------------------------------
setClass("subsettingView",
         contains="view",
         prototype=prototype(ID=paste("subViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the subsettingView object and directly assigns it
## to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
subsettingView <- function(workflow, ID=paste("subViewRef", guid(), sep="_"),
                           name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("subsettingView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a subsetting object. This creates a
## subsettingActionItem in the workFlow and from that a subsettingView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="subsetting"),
          definition=function(wf, action, parent=NULL,
                              name=identifier(action))
      {
           fun <- function(wf, action, parent, name)
           {
               ## get the parentView. If not explicitely specified, use the
               ## root node
               pid <- if(is.null(parent)) views(wf)[[1]] else parent
               pid <- getAlias(pid, wf)
               if(is.null(unlist(pid)))
                   stop("'", parent, "' is not a valid view name in this",
                        " workflow.", call.=FALSE)
               ## assign the subsetting to the evaluation environment and
               ## create a reference to it
               subRef <- assign(value=action, envir=wf)
               pview <- fcViewReference(ID=pid, env=wf)
               ## create and assign a new ActionItem
               actionRef <- subsettingActionItem(parentView=pview,
                                                 subsetting=subRef,
                                                 workflow=wf)
               ## add a useful name to the journal entry
               journal <- get(wf@journal)
               names(journal)[length(journal)] <- identifier(actionRef)
               assign(wf@journal, journal, wf)
               ## now subset the data and assign the result
               pdat <- Data(get(pview))
               if(is.character(action@indices))
               {
                   if(!all(action@indices %in% sampleNames(pdat)))
                       stop("Subset out of bounds.")
                   action@indices <- match(action@indices, sampleNames(pdat))
               }
               tData <- Data(get(pview))[action@indices]
               dataRef <- assign(value=tData, envir=wf)
               vid <- subsettingView(name=name, workflow=wf,
                                     action=actionRef, data=dataRef)
               ## update the identifier of the subsetting object
               identifier(action) <- paste("sub", identifier(action),
                                           sep="_")
               assign(subRef, value=action, envir=wf) 
               ## add new nodes and edges to the workflow tree
               nid <- identifier(vid)
               tree <- get(wf@tree)
               tree <- addNode(nid, tree)
               tree <- addEdge(pid, identifier(vid), tree)
               edgeDataDefaults(tree, "actionItem") <- fcNullReference()
               edgeData(tree, pid , nid, "actionItem") <- actionRef
               assign(x=wf@tree, value=tree, envir=wf)
               return(invisible(wf))
           }
           ## create a new journal entry first
           journal <- c(get(wf@journal), list(NULL))
           assign(wf@journal, journal, wf)
           ## now do the actual adding, new objects are captured in the journal
           tryCatch(fun(wf=wf, action=action, parent=parent, name=name),
                    error=function(e){
                        cat("Error adding item '", substitute(action),
                            "':\n  ", sep="")
                        cat(gsub("Error.*\\)\\:", "", e))
                        undo(wf)
                    })
           return(invisible(wf))
       })


setMethod("add",
          signature=signature(wf="workFlow", action="character"),
          definition=function(wf, action, ...)
      {
          action <- subsetting(action)
          add(wf, action,...)
      })

setMethod("add",
          signature=signature(wf="workFlow", action="numeric"),
          definition=function(wf, action, ...)
      {
          action <- subsetting(action)
          add(wf, action,...)
      })

setMethod("add",
          signature=signature(wf="workFlow", action="logical"),
          definition=function(wf, action, ...)
      {
          action <- subsetting(which(action))
          add(wf, action,...)
      })




## ===========================================================================
## Unity transformation
## ---------------------------------------------------------------------------
## Transforms parameters names provided as characters into unity transform 
## objects which can be evaluated to retrive the corresponding columns from the
## data frame
## ---------------------------------------------------------------------------
setClass("unitytransform",
	 contains="transform",
	 representation=representation(parameters="character"))

unitytransform <- function(parameters,
                           transformationId="defaultUnityTransform")
{
    checkClass(transformationId, "character", 1)
    if(missing(parameters))
        parameters <- character()
    new("unitytransform", parameters=parameters,
        transformationId=transformationId)
}



## ===========================================================================
## Polynomial transformation of degree 1 
## ---------------------------------------------------------------------------
## Allows for scaling ,linear combination and translation within a single 
## transformation
## ---------------------------------------------------------------------------
setClass("dg1polynomial", 		
         contains="transform",
         representation=representation(parameters="parameters",
                                       a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=new("parameters"),
                             a=1,
                             b=1))

dg1polynomial <- function(parameters, a=1, b=1,
                          transformationId="defaultDg1polynomialTransform")
{
    checkClass(a, "numeric", length(parameters))
    checkClass(b, "numeric", 1)
    checkClass(transformationId, "character", 1)
    new("dg1polynomial", parameters=parameters, a=a, b=b,
        transformationId=transformationId)
}



## ===========================================================================
## Ratio transformation
## ---------------------------------------------------------------------------
## Ratio of two arguments defined in the transformation
## ---------------------------------------------------------------------------
setClass("ratio",
         contains="transform",
         representation(numerator="transformation",
                        denominator="transformation"),
	 prototype=prototype(numerator=unitytransform(),
                             denominator=unitytransform()))

ratio <- function(numerator=unitytransform(),
                  denominator=unitytransform(),
                  transformationId="defaultRatioTransform")
{
    if(!is(numerator, "transform")){
        checkClass(numerator, "character", 1)
        numerator <- unitytransform(numerator)
    }
    if(!is(denominator, "transform")){
        checkClass(denominator, "character", 1)
        denominator=unitytransform(denominator)
    }  
    new("ratio", numerator=numerator, denominator=denominator,
        transformationId=transformationId)
}



## ===========================================================================
## Quadratic transformation
## ---------------------------------------------------------------------------
setClass("quadratic", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1),
         validity=function(object) 
     {
         msg<-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Quadratic transform is defined for one parameter")
         if(length(object@a)!=1)
             msg <- c(msg, "Only one coefficient is defined for quadratic transform")
         if(object@a==0)
             msg <- c(msg, "'a' should be non zero")
         msg
     })

quadratic <- function(parameters="NULL", a=1,
                      transformationId="defaultQuadraticTransform")
    new("quadratic",parameters=parameters,a=a,
        transformationId=transformationId)

          

## ===========================================================================
## Squareroot transformation
## ---------------------------------------------------------------------------
setClass("squareroot", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Square root transform is defined for one parameter")
         if(length(object@a)!=1)
             msg <- c(msg, "Only one coefficient is defined for quadratic transform")
         if(object@a==0)
             msg <- c(msg, "Coefficient should be non zero")
         msg
     })

squareroot <- function(parameters, a=1,
                       transformationId="defaultSquarerootTransform")
    new("squareroot", parameters=parameters, a=a,
        transformationId=transformationId)



## ===========================================================================
##  Loarithmical Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("logarithm",
         contains="singleParameterTransform",
         representation=representation(a="numeric", b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Logarithm transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

logarithm <- function(parameters, a=1, b=1,
                      transformationId="defaultLogarithmTransform")
    new("logarithm", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Exponential Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------

setClass("exponential", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1,
                             b=1),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Exponential transform is defined for one parameter")
         if(object@a==0)
             msg<-c(msg,"'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg,"'b' should be a non zero number")
         msg  
     })

exponential <- function(parameters, a=1, b=1,
                        transformationId="defaultExponentialTransformation")
    new("exponential", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Inverse hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("asinht", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Inverse hypberbolic transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

asinht <- function(parameters="NULL", a=1, b=1,
                   transformationId="defaultAsinhTransform")
    new("asinht", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Inverse hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("sinht", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1,
                             b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Hypberbolic transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

sinht <- function(parameters, a=1, b=1,
                  transformationId="defaultSinhtTransform")
    new("sinht", parameters=parameters, a=a, b=b,
        transformationId=transformationId)




## ================================================================================
## Inverse hyperbolic sin transformation parametrized according to Gating-ML 2.0 
## --------------------------------------------------------------------------------
## Inputs T, M, A of type numeric and parameter of type transformation or character
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## --------------------------------------------------------------------------------
setClass(
    "asinhtGml2", 		
    contains = "singleParameterTransform",
    representation = representation(T = "numeric", M = "numeric", A = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        parameters = unitytransform(),
        T = 262144,
        M = 4.5,
        A = 0,
        boundMin = -Inf,
        boundMax = Inf),
    validity = function(object)
    {
        msg <- NULL
        if (length(object@parameters) != 1)
            msg <- c(msg, "Inverse hyperbolic sin transformation is defined for one parameter.")
        if (object@T <= 0)
            msg <- c(msg, "'T' should be greater than zero.")
        if (object@M <= 0)
            msg <- c(msg, "'M' should be greater than zero.")
        if (object@A < 0)
            msg <- c(msg, "'A' should be greater than or equal to zero.")
        if (object@A > object@M)
            msg <- c(msg, "'A' should be less than or equal to 'M'.")
        if (object@boundMin > object@boundMax)
          msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
        msg
    }
)

asinhtGml2 <- function(
        parameters, 
        T = 262144, 
        M = 4.5, 
        A = 0, 
        transformationId = "defaultAsinhGml2Transform",
        boundMin = -Inf,
        boundMax = Inf)
    new("asinhtGml2", parameters = parameters, 
        T = T, M = M, A = A, transformationId = transformationId, boundMin = boundMin, boundMax = boundMax)


## ===================================================================================
## Logicle transformation parametrized according to Gating-ML 2.0
## -----------------------------------------------------------------------------------
## Inputs T, M, W, A of type numeric and parameter of type transformation or character
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## -----------------------------------------------------------------------------------
setClass(
    "logicletGml2", 		
    contains = "singleParameterTransform",
    representation = representation(T = "numeric", M = "numeric", W = "numeric", A = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        parameters = unitytransform(),
            T = 262144,
            M = 4.5,
            W = 0.5,
            A = 0,
            boundMin = -Inf,
            boundMax = Inf),
    validity = function(object)
    {
        msg <- NULL
        if (length(object@parameters) != 1)
            msg <- c(msg, "Logicle transformation is defined for one parameter.")
        if (object@T <= 0)
            msg <- c(msg, "'T' should be greater than zero.")
        if (object@M <= 0)
            msg <- c(msg, "'M' should be greater than zero.")
        if (object@W < 0)
            msg <- c(msg, "'W' should be greater than or equal to zero.")
        if (object@W > object@M/2)
            msg <- c(msg, "'W' should be less than or equal to half of 'M'.")
        if (object@A < -object@W)
            msg <- c(msg, "'A' should be greater than or equal to 'minus W'.")
        if (object@A > object@M - 2*object@W)
            msg <- c(msg, "'A' should be less than or equal to 'M minus two W'")
        if (object@boundMin > object@boundMax)
          msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
        msg
    }
)

logicletGml2 <- function(
    parameters,
    T = 262144,
    M = 4.5,
    W = 0.5,
    A = 0,
    transformationId = "defaultLogicletGml2Transform",
    boundMin = -Inf,
    boundMax = Inf)
    new("logicletGml2", parameters = parameters,
        T = T, M = M, W = W, A = A, transformationId = transformationId, boundMin = boundMin, boundMax = boundMax)


## ===================================================================================
## Hyperlog transformation parametrized according to Gating-ML 2.0
## -----------------------------------------------------------------------------------
## Inputs T, M, W, A of type numeric and parameter of type transformation or character
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## -----------------------------------------------------------------------------------
setClass(
    "hyperlogtGml2",
    contains = "singleParameterTransform",
    representation = representation(T = "numeric", M = "numeric", W = "numeric", A = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        parameters = unitytransform(),
        T = 262144,
        M = 4.5,
        W = 0.5,
        A = 0,
        boundMin = -Inf,
        boundMax = Inf),
    validity = function(object)
    {
        msg <- NULL
        if (length(object@parameters) != 1)
            msg <- c(msg, "Logicle transformation is defined for one parameter.")
        if (object@T <= 0)
            msg <- c(msg, "'T' should be greater than zero.")
        if (object@M <= 0)
            msg <- c(msg, "'M' should be greater than zero.")
        if (object@W <= 0)
            msg <- c(msg, "'W' should be greater than zero.")
        if (object@W > object@M/2)
            msg <- c(msg, "'W' should be less than or equal to half of 'M'.")
        if (object@A < -object@W)
            msg <- c(msg, "'A' should be greater than or equal to 'minus W'.")
        if (object@A > object@M - 2*object@W)
            msg <- c(msg, "'A' should be less than or equal to 'M minus two W'")
        if (object@boundMin > object@boundMax)
          msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
        msg
    }
)

hyperlogtGml2 <- function(
    parameters,
    T = 262144,
    M = 4.5,
    W = 0.5,
    A = 0,
    transformationId = "defaultHyperlogtGml2Transform",
    boundMin = -Inf,
    boundMax = Inf)
    new("hyperlogtGml2", parameters = parameters,
        T = T, M = M, W = W, A = A, transformationId = transformationId, boundMin = boundMin, boundMax = boundMax)

## ================================================================================
## Linear transformation parametrized according to Gating-ML 2.0
## --------------------------------------------------------------------------------
## Inputs T, A of type numeric and parameter of type transformation or character
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## --------------------------------------------------------------------------------
setClass(
    "lintGml2",
    contains = "singleParameterTransform",
    representation = representation(T = "numeric", A = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        parameters = unitytransform(),
        T = 262144,
        A = 0,
        boundMin = -Inf,
        boundMax = Inf),
    validity = function(object)
    {
        msg <- NULL
        if (length(object@parameters) != 1)
            msg <- c(msg, "Linear transformation is defined for one parameter.")
        if (object@T <= 0)
            msg <- c(msg, "'T' should be greater than zero.")
        if (object@A < 0)
            msg <- c(msg, "'A' should be greater than or equal to zero.")
        if (object@A > object@T)
            msg <- c(msg, "'A' should be less than or equal to 'T'.")
        if (object@boundMin > object@boundMax)
          msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
        msg
    }
)

lintGml2 <- function(
    parameters,
    T = 262144,
    A = 0,
    transformationId = "defaultLintGml2Transform",
    boundMin = -Inf,
    boundMax = Inf)
    new("lintGml2", parameters = parameters,
        T = T, A = A, transformationId = transformationId, boundMin = boundMin, boundMax = boundMax)


## ================================================================================
## Log transformation parametrized according to Gating-ML 2.0
## --------------------------------------------------------------------------------
## Inputs T, M of type numeric and parameter of type transformation or character
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## --------------------------------------------------------------------------------
setClass(
    "logtGml2",
    contains = "singleParameterTransform",
    representation = representation(T = "numeric", M = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        parameters = unitytransform(),
        T = 262144,
        M = 4.5,
        boundMin = -Inf,
        boundMax = Inf),
    validity = function(object)
    {
        msg <- NULL
        if (length(object@parameters) != 1)
            msg <- c(msg, "Log transformation is defined for one parameter.")
        if (object@T <= 0)
            msg <- c(msg, "'T' should be greater than zero.")
        if (object@M <= 0)
            msg <- c(msg, "'M' should be greater than zero.")
        if (object@boundMin > object@boundMax)
          msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
        msg
    }
)

logtGml2 <- function(
    parameters,
    T = 262144,
    M = 4.5,
    transformationId = "defaultLogGml2Transform",
    boundMin = -Inf,
    boundMax = Inf)
    new("logtGml2", parameters = parameters,
        T = T, M = M, transformationId = transformationId, boundMin = boundMin, boundMax = boundMax)



## ========================================================================================
## Ratio transformation parametrized according to Gating-ML 2.0
## ----------------------------------------------------------------------------------------
## Inputs A, B and C of type numeric and two parameters of type character or transformation
##
## October 2014: additional boundMin and boundMax attributes to all Gating-ML 2.0
## transforms; if the result of the transform is outside of that range then set it
## to the appropriate boundMin/boundMax.
## ----------------------------------------------------------------------------------------
setClass("ratiotGml2",
    contains="transform",
    representation(
        numerator = "transformation", denominator = "transformation",
        pA = "numeric", pB = "numeric", pC = "numeric", boundMin = "numeric", boundMax = "numeric"),
    prototype = prototype(
        numerator=unitytransform(),
        denominator=unitytransform(),
        pA = 1,
        pB = 0,
        pC = 0,
        boundMin = -Inf,
        boundMax = Inf),
    validity = function(object)
    {
      msg <- NULL
      if (object@boundMin > object@boundMax)
        msg <- c(msg, "'boundMin' should be less than or equal to 'boundMax'")
      msg
    }
)

ratiotGml2 <- function(
    numerator = unitytransform(),
    denominator = unitytransform(),
	pA = 1,
	pB = 0,
	pC = 0,
	transformationId = "defaultRatioTransform",
	boundMin = -Inf,
	boundMax = Inf)
{
    if(!is(numerator, "transform")){
        checkClass(numerator, "character", 1)
        numerator <- unitytransform(numerator)
    }
    if(!is(denominator, "transform")){
        checkClass(denominator, "character", 1)
        denominator <- unitytransform(denominator)
    }
    new("ratiotGml2", numerator = numerator, denominator = denominator,
        pA = pA, pB = pB, pC = pC, transformationId = transformationId, 
        boundMin = boundMin, boundMax = boundMax)
}


## ===========================================================================
##  Hyperlog Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## --------------------------------------------------------------------------- 
setClass("hyperlog", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Hyperlog transform is defined for one parameter")
         if(object@a<=0)
             msg <- c(msg, "'a' should be greater than zero")
         if(object@b<=0)
             msg <- c(msg, "'b' should be greater than zero")
         msg
     })

hyperlog <- function(parameters="NULL", a=1, b=1,
                     transformationId="defaultHyperlogTransform")
    new("hyperlog", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  EH Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## --------------------------------------------------------------------------- 
setClass("EHtrans", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "EH transform is defined for one parameter")
         if(object@a<=0)
             msg<-c(msg, "'a' should be greater than zero")
         if(object@b<=0)
             msg<-c( msg, "'b' should be greater than zero")
         msg
     })

EHtrans <- function(parameters, a=1, b=1,
                    transformationId="defaultEHtransTransform")
    new("EHtrans", parameters=parameters, a=a, b=b,
        transformationId=transformationId)
      
          

## ===========================================================================
##  Splitscale Transformation 
## ---------------------------------------------------------------------------
setClass("splitscale", 		
         contains="singleParameterTransform",
         representation=representation(r="numeric",
                                       maxValue="numeric",
                                       transitionChannel="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             r=1,
                             maxValue=1,
                             transitionChannel=4),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg<-c(msg, "Split scale transform is defined for one parameter")
         if(object@r<=0)
             msg <- c(msg, "'r' should be a greater than zero")
         if(object@maxValue<=0)
             msg <- c(msg, "maxValue should be a greater than zero")
         if(object@transitionChannel<0)
             msg <- c(msg, "transitionChannel should be a non negative")
         msg
     })

splitscale <- function(parameters="NULL", r=1, maxValue=1, transitionChannel=4,
                       transformationId="defaultSplitscaleTransform")
    new("splitscale",
        parameters=parameters, r=r, maxValue=maxValue,
        transitionChannel=transitionChannel,
        transformationId=transformationId)



## ===========================================================================
##  Inverse Splitscale Transformation 
## ---------------------------------------------------------------------------
setClass("invsplitscale", 		
         contains="singleParameterTransform",
         representation=representation(r="numeric",
                                       maxValue="numeric",
                                       transitionChannel="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             r=1,
                             maxValue=1,
                             transitionChannel=4),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Split scale transform is defined for one parameter")
         if(object@r<=0)
             msg <- c(msg, "'r' should be a greater than zero")
         if(object@maxValue<=0)
             msg <- c(msg, "'maxValue' should be a greater than zero")
         if(object@transitionChannel<0)
             msg <- c(msg, "'transitionChannel' should be a non negative")
         msg
     })
       
invsplitscale <- function(parameters, r=1, maxValue=1,
                          transitionChannel=4,
                          transformationId="defaultInvsplitscaleTransforms")
    new("invsplitscale",
        parameters=parameters, r=r, maxValue=maxValue,
        transitionChannel=transitionChannel,
        transformationId=transformationId)
      


## ===========================================================================
## Transformation reference
## ---------------------------------------------------------------------------
## Reference to a transformation defined previously
## ---------------------------------------------------------------------------
setClass("transformReference",
         contains="transform",
         representation(searchEnv="environment"))

transformReference <- function(referenceId="defaultTransformReference",
                               searchEnv)
    new("transformReference",
        transformationId=referenceId, searchEnv=searchEnv)
    
