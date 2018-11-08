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
#'     with \code{$}.
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
#'   \item{decompensate}{Reverse the application of a compensation matrix (or a
#'     \code{\linkS4class{compensation}} object) on a \code{flowFrame}
#'     object. This returns a decompensated \code{flowFrame}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   decompensate(flowFrame, matrix)}
#'     \code{   decompensate(flowFrame, data.frame)}
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export 
setClass("transform",
         representation=representation(transformationId="character",
                                       .Data="function"),
         prototype=prototype(transformationId=""))

#' Class "parameters"
#' 
#' A representation of flow parameters that allows for referencing.
#' 
#' 
#' @name parameters-class
#' @aliases parameters-class
#' @docType class
#' @section Objects from the Class: Objects will be created internally whenever
#' necessary and this should not be of any concern to the user.
#' 
#' @slot .Data A list of the individual parameters.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{list}"}, from data part.
#' Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
#' 
#' @author Nishant Gopalakrishnan
#' @keywords classes
#'
#' @export
setClass("parameters", contains="list")

#' Class "transformation"
#' 
#' A virtual class to abstract transformations.
#' 
#' 
#' @name transformation-class
#' @aliases transformation-class transformation
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @section Extends:
#' Class \code{"\linkS4class{characterOrTransformation}"}, directly.
#' @author N. Gopalakrishnan
#' @keywords classes
setClassUnion("transformation", "transform")


#' Class "characterOrTransformation"
#' 
#' A simple union class of \code{character} and \code{\linkS4class{transformation}}.
#' Objects will be created internally whenever necessary and the user should
#' not need to explicitly interact with this class.
#' 
#' @name characterOrTransformation-class
#' @aliases characterOrTransformation-class characterOrTransformation
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @keywords classes
#' @examples
#' 
#' showClass("characterOrTransformation")
#' 
setClassUnion("characterOrTransformation", c("character","transformation"))

#' Class "characterOrParameters"
#' 
#' A simple union class of \code{character} and \code{\linkS4class{parameters}}.
#' Objects will be created internally whenever necessary and the user should
#' not need to explicitly interact with this class.
#' 
#' @name characterOrParameters-class
#' @aliases characterOrParameters-class characterOrParameters
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @keywords classes
#' @examples
#' 
#' showClass("characterOrParameters")
#' 
setClassUnion("characterOrParameters", c("character","parameters"))

#' Class "singleParameterTransform"
#' 
#' A transformation that operates on a single parameter
#' 
#' 
#' @name singleParameterTransform-class
#' @aliases singleParameterTransform-class
#' initialize,singleParameterTransform-method
#' parameters,singleParameterTransform-method
#' @docType class
#' @section Objects from the Class:
#' 
#' Objects can be created by calls of the form
#' \code{new("singleParameterTransform", ...)}.
#' 
#' @slot .Data Object of class \code{"function"}. The transformation.
#' @slot parameters Object of class \code{"transformation"}. The 
#' parameter to transform. Can be a derived parameter from another 
#' transformation.
#' @slot transformationId Object of class \code{"character"}. An 
#' identifier for the object.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author F Hahne
#' @keywords classes
#' @examples
#' 
#' showClass("singleParameterTransform")
#' 
setClass("singleParameterTransform",
         representation=representation(parameters="transformation"),
         contains="transform")

#' Class "nullParameter"
#' 
#' A class used internally for coercing transforms to characters for a return
#' value when a coercion cannot be performed. The user should never need to
#' interact with this class.
#' 
#' @name nullParameter-class
#' @aliases nullParameter-class nullParameter
#' @docType class
#' @section Objects from the Class: Objects will be created internally whenever
#' necessary and this should not be of any concern to the user.
#' @keywords classes
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
#'   more metadata about the results of the filtering operation. See 
#'   the documenation for \code{\link[flowCore:filter-methods]{filter}} 
#'   methods for details.}
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
#'
#' @export
setClass("filter", 
         representation=representation("VIRTUAL",
         filterId="character"),
         prototype=prototype(filterId=""))

#' Class "concreteFilter"
#' 
#' The \code{concreteFilter} serves as a base class for all filters that
#' actually implement a filtering process. At the moment this includes all
#' filters except \code{\linkS4class{filterReference}}, the only non-concrete
#' filter at present.
#' 
#' 
#' @name concreteFilter-class
#' @aliases concreteFilter-class concreteFilter
#' @docType class
#' @section Objects from the Class: Objects of this class should never be
#' created directly. It serves only as a point of inheritance.
#' 
#' @slot filterId The identifier associated with this class.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{filter}"}, directly.
#' 
#' @author B. Ellis
#' @seealso \code{\linkS4class{parameterFilter}}
#' @keywords classes
#'
#' @export
setClass("concreteFilter",
         contains="filter")

                                        # setClass("parameterFilter",
                                        #          representation=representation(parameters="character"),
                                        #          contains="concreteFilter",
                                        #          prototype=prototype(parameters=""))
#' Class "parameterFilter"
#' 
#' A concrete filter that acts on a set of parameters.
#' 
#' 
#' @name parameterFilter-class
#' @aliases parameterFilter-class initialize,parameterFilter-method
#' @docType class
#' @section Objects from the Class: \code{parameterFilter} objects are never
#' created directly. This class serves as an inheritance point for filters that
#' depends on particular parameters.
#' 
#' @slot parameters The names of the parameters employed by this filter.
#' @slot filterId The filter identifier.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{concreteFilter}"}, directly.
#' Class \code{"\linkS4class{filter}"}, by class "concreteFilter", distance 2.
#' 
#' @author B. Ellis
#' @keywords classes
#'
#' @export
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
#' filters(x)
#' 
#' filtersList(x)
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
#' 
#' @export
setClass("filters",
		 contains="list"
		 )
 ## Constructor
#' @export
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
 
#' @export
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
#' @export
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
#' @family Gate classes
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
#' @export
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
#' @export
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
#' @family Gate classes
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{flowSet}}, \code{\link{filter}} for
#' evaluation of \code{quadGates} and \code{\link{split}} for splitting of flow
#' cytometry data sets based on that.
#' @keywords classes methods
#' @examples
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
#' @export
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
#' @export
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
#' @family Gate classes
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
#' @export
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
#' @export
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
#' @family Gate classes
#' @seealso \code{\link{flowFrame}}, \code{\link{filter}}
#' @keywords methods
#'
#' @export
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
#' @export
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
#' @family Gate classes
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
#' @export
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
#' @export
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
#' Class "norm2Filter"
#' 
#' 
#' Class and constructors for a \code{\link{filter}} that fits a bivariate
#' normal distribution to a data set of paired values and selects data points
#' according to their standard deviation from the fitted distribution.
#' 
#' 
#' The filter fits a bivariate normal distribution to the data and selects all
#' events within the Mahalanobis distance multiplied by the \code{scale.factor}
#' argument. The constructor \code{norm2Filter} is a convenience function for
#' object instantiation. Evaluating a \code{curv2Filter} results in an object
#' of class \code{\link{logicalFilterResult}}. Accordingly, \code{norm2Filters}
#' can be used to subset and to split flow cytometry data sets.
#' 
#' @name norm2Filter-class
#' @aliases norm2Filter-class norm2Filter show,norm2Filter-method
#' @docType class
#' @usage
#' norm2Filter(x, y, method="covMcd", scale.factor=1, n=50000,
#' filterId="defaultNorm2Filter")
#' @param x,y Characters giving the names of the measurement parameter on which
#' the filter is supposed to work on. \code{y} can be missing in which case
#' \code{x} is expected to be a character vector of length 2 or a list of
#' characters.
#' @param filterId An optional parameter that sets the \code{filterId} slot of
#' this filter. The object can later be identified by this name.
#' @param scale.factor,n Numerics of length 1, used to set the
#' \code{scale.factor} and n slots of the object.
#' @param method Character in \code{covMcd} or \code{cov.rob}, used to set the
#' \code{method} slot of the object.
#' @return
#' Returns a \code{\link{norm2Filter}} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' 
#' @note
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{norm2Filters}.
#' 
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
#' @slot method One of \code{covMcd} or \code{cov.rob}
#' defining method used for computation of covariance matrix.
#' @slot scale.factor Numeric vector giving factor of standard
#' deviations used for data selection (all points within
#' \code{scalefac} standard deviations are selected).
#' @slot n Object of class \code{"numeric"}, the number of
#' events used to compute the covariance matrix of the bivariate
#' distribution.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter.
#' @slot parameters Object of class \code{"ANY"} describing
#' the parameters used to filter the \code{\link{flowFrame}} or
#' \code{\link{flowSet}}.
#' 
#' @section Objects from the Class:
#'   Objects can be created by calls of the form \code{new("norm2Filter",
#' ...)} or using the constructor \code{norm2Filter}. The constructor
#' is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "norm2Filter")}: The workhorse used to evaluate the filter on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "norm2Filter")}: Print
#'     information about the filter. }
#'   
#' }
#' 
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[MASS]{cov.rob}}, \code{\link[rrcov]{CovMcd}},
#' \code{\link[flowCore:filter-methods]{filter}} for evaluation of
#' \code{norm2Filters} and \code{\link{split}} and \code{\link{Subset}}for
#' splitting and subsetting of flow cytometry data sets based on that.
#' @keywords classes methods
#' @examples
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Create directly. Most likely from a command line
#' norm2Filter("FSC-H", "SSC-H", filterId="myCurv2Filter")
#' 
#' ## To facilitate programmatic construction we also have the following
#' n2f <- norm2Filter(filterId="myNorm2Filter", x=list("FSC-H", "SSC-H"),
#' scale.factor=2)
#' n2f <- norm2Filter(filterId="myNorm2Filter", x=c("FSC-H", "SSC-H"),
#' scale.factor=2)
#' 
#' ## Filtering using norm2Filter
#' fres <- filter(dat, n2f)
#' fres
#' summary(fres)
#' 
#' ## The result of norm2 filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' 
#' @export
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
#' @export
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
#' Class "kmeansFilter"
#' 
#' 
#' A filter that performs one-dimensional k-means (Lloyd-Max) clustering on a
#' single flow parameter.
#' 
#' 
#' The one-dimensional k-means filter is a multiple population filter capable
#' of operating on a single flow parameter. It takes a parameter argument
#' associated with two or more populations and results in the generation of an
#' object of class \code{\link{multipleFilterResult}}.  Populations are
#' considered to be ordered such that the population with the smallest mean
#' intensity will be the first population in the list and the population with
#' the highest mean intensity will be the last population listed.
#' 
#' @name kmeansFilter-class
#' @aliases kmeansFilter kmeansFilter-class length,kmeansFilter-method
#' show,kmeansFilter-method
#' @docType class
#' @usage 
#' kmeansFilter(\dots, filterId="defaultKmeansFilter")
#' @param \dots \code{kmeansFilter} are defined by a single flow parameter and
#' an associated list of \code{k} population names. They can be given as a
#' character vector via a named argument, or as a list with a single named
#' argument. In both cases the name will be used as the flow parameter and the
#' content of the list or of the argument will be used as population names,
#' after coercing to character. For example
#' 
#' \code{kmeansFilter(FSC=c("a", "b", "c"))}
#' 
#' or
#' 
#' \code{kmeansFilter(list(SSC=1:3))}
#' 
#' If the parameter is not fully realized, but instead is the result of a
#' \code{\link[flowCore:transform-class]{transformation}} operation, two
#' arguments need to be passed to the constructor: the first one being the
#' \code{\link[flowCore:transform-class]{transform}} object and the second
#' being a vector of population names which can be coerced to a character. For
#' example
#' 
#' \code{kmeansFilter(tf, c("D", "E"))}
#' 
#' @param filterId An optional parameter that sets the \code{filterId} of the
#' object. The filter can later be identified by this name.
#' @return
#' 
#' Returns a \code{kmeansFilter} object for use in filtering
#' \code{\link[flowCore:flowFrame-class]{flowFrames}} or other flow cytometry
#' objects.
#' @note
#' 
#' See the documentation in the \code{\link[flowViz:flowViz-package]{flowViz}}
#' package for plotting of \code{kmeansFilters}.
#' @section Extends:
#' 
#' Class \code{\linkS4class{parameterFilter}}, directly.
#' 
#' Class \code{\linkS4class{concreteFilter}}, by class \code{parameterFilter},
#' distance 2.
#' 
#' Class \code{\linkS4class{filter}}, by class \code{parameterFilter},
#' distance3.
#' 
#' @slot populations Object of class \code{character}. The
#' names of the \code{k} populations (or clusters) that will be
#' created by the \code{kmeansFilter}. These names will later be used
#' for the respective subpopulations in \code{\link{split}}
#' operations and for the summary of the \code{\link{filterResult}}.
#' @slot parameters Object of class \code{\link{parameters}},
#' defining a single parameter for which the data in the
#' \code{\linkS4class{flowFrame}} is to be clustered. This may also
#' be a \code{\link[flowCore:transform-class]{transformation}} object.
#' @slot filterId Object of class \code{character}, an
#' identifier or name to reference the \code{kmeansFilter} object
#' later on.
#' 
#' @section Objects from the Class:
#' Like all other \code{\linkS4class{filter}} objects in \code{flowCore},
#' \code{kmeansFilter} objects should be instantiated through their
#' constructor \code{kmeansFilter()}. See the \code{Usage} section for
#' details.
#' 
#' @section Methods:
#' 
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "kmeansFilter")}: The workhorse used to evaluate the filter on
#'     data.
#'     
#'     \emph{Usage:}
#'     
#'     This is usually not called directly by the user, but internally by
#'     the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "kmeansFilter")}: Print
#'     information about the filter.
#'     
#'     \emph{Usage:}
#'     
#'     The method is called automatically whenever the object is printed
#'     on the screen. }
#'   
#' }
#' 
#' 
#' @author F. Hahne, B. Ellis, N. LeMeur
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{flowSet}}, \code{\link{filter}} for
#' evaluation of \code{kmeansFilters} and \code{\link{split}} for splitting of
#' flow cytometry data sets based on the result of the filtering operation.
#' @keywords methods classes
#' @examples
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Create the filter
#' kf <- kmeansFilter("FSC-H"=c("Pop1","Pop2","Pop3"), filterId="myKmFilter")
#' 
#' ## Filtering using kmeansFilters
#' fres <- filter(dat, kf)
#' fres
#' summary(fres)
#' names(fres)
#' 
#' ## The result of quadGate filtering are multiple sub-populations
#' ## and we can split our data set accordingly
#' split(dat, fres)
#' 
#' ## We can limit the splitting to one or several sub-populations
#' split(dat, fres, population="Pop1")
#' split(dat, fres, population=list(keep=c("Pop1","Pop2")))
#' 
#' 
#' @export
setClass("kmeansFilter",
         representation=representation(populations="character"),
         prototype=list(filterId="defaultKmeansFilter"),
         contains="parameterFilter")

## Constructor. We allow for the following inputs:
##  ..1 is transform and .2 is some vector that can be coerced to character
##  ..1 is some vector that can be coerced to character
#' @export
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
#' Class "sampleFilter"
#' 
#' 
#' This non-parameter filter selects a number of events from the primary
#' \code{\link{flowFrame}}.
#' 
#' 
#' Selects a number of events without replacement from a \code{flowFrame}.
#' 
#' @name sampleFilter-class
#' @aliases sampleFilter-class sampleFilter show,sampleFilter-method
#' @docType class
#' @usage 
#' sampleFilter(size, filterId="defaultSampleFilter")
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' \code{\link{filter}}. The object can later be identified by this name.
#' @param size The number of events to select.
#' @return
#' 
#' Returns a \code{sampleFilter} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{concreteFilter},
#' distance 2.
#' 
#' @slot size Object of class \code{"numeric"}. Then number of
#' events that are to be selected.
#' @slot filterId A character vector that identifies this
#' \code{filter}.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("sampleFilter",
#' ...)} or using the constructor \code{sampleFilter}. The latter is the
#' recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "sampleFilter")}: The workhorse used to evaluate the gate on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "sampleFilter")}: Print
#'     information about the gate. }
#'   
#' }
#' 
#' @author B. Ellis, F.Hahne
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{filter}} for evaluation of
#' \code{sampleFilters} and \code{\link{split}} and \code{\link{Subset}}for
#' splitting and subsetting of flow cytometry data sets based on that.
#' @keywords methods classes
#' @examples
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' #Create the filter
#' sf <- sampleFilter(filterId="mySampleFilter", size=500)
#' sf
#' 
#' ## Filtering using sampleFilters
#' fres <- filter(dat, sf)
#' fres
#' summary(fres)
#' 
#' ## The result of sample filtering is a logical subset
#' Subset(dat, fres)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' 
#' @export
setClass("sampleFilter",
         representation=representation(size="numeric"),
         contains="concreteFilter",
         prototype=list(size=10000, filterId="defaultSampleFilter"))

##Constructor: We allow for the following inputs:
##  size is always a numeric of length 1
#' @export
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
#' Class "boundaryFilter"
#' 
#' 
#' Class and constructor for data-driven \code{\link{filter}} objects that
#' discard margin events.
#' 
#' 
#' Flow cytomtery instruments usually operate on a given data range, and the
#' limits of this range are stored as keywords in the FSC files. Depending on
#' the amplification settings and the dynamic range of the measured signal,
#' values can occur that are outside of the measurement range, and most
#' instruments will simply pile those values at the minimum or maximum range
#' limit. The \code{boundaryFilter} removes these values, either for a single
#' parameter, or for a combination of parameters. Note that it is often
#' desirable to treat boundary events on a per-parameter basis, since their
#' values might be uninformative for one particular channel, but still be
#' useful in all of the other channels.
#' 
#' The constructor \code{boundaryFilter} is a convenience function for object
#' instantiation. Evaluating a \code{boundaryFilter} results in a single
#' sub-populations, an hence in an object of class \code{\link{filterResult}}.
#' 
#' @name boundaryFilter-class
#' @aliases boundaryFilter-class boundaryFilter show,boundaryFilter-method
#' @docType class
#' @usage 
#' boundaryFilter(x, tolerance=.Machine$double.eps, side=c("both", "lower",
#' "upper"), filterId="defaultBoundaryFilter")
#' @param x Character giving the name(s) of the measurement parameter(s) on
#' which the filter is supposed to work. Note that all events on the margins of
#' ay of the channels provided by \code{x} will be discarded, which is often
#' not desired. Such events may not convey much information in the particular
#' channel on which their value falls on the margin, however they may well be
#' informative in other channels.
#' @param tolerance Numeric vector, used to set the \code{tolerance} slot of
#' the object. Can be set separately for each element in \code{x}. R's
#' recycling rules apply.
#' @param side Character vector, used to set the \code{side} slot of the
#' object.  Can be set separately for each element in \code{x}. R's recycling
#' rules apply.
#' @param filterId An optional parameter that sets the \code{filterId} slot of
#' this filter. The object can later be identified by this name.
#' @return
#' 
#' Returns a \code{boundaryFilter} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
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
#' @slot tolerance Object of class \code{"numeric"}. The
#' machine tolerance used to decide whether an event is on the
#' measurement boundary. Essentially, this is done by evaluating
#' \code{x>minRange+tolerance & x<maxRange-tolerance}.
#' @slot side Object of class \code{"character"}. The margin
#' on which to evaluate the filter. Either \code{upper} for the
#' upper margin or \code{lower} for the lower margin or \code{both}
#' for both margins.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("boundaryFilter",
#' ...)} or using the constructor \code{boundaryFilter}.  Using the
#' constructor is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "boundaryFilter")}: The workhorse used to evaluate the filter on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "boundaryFilter")}: Print
#'     information about the filter. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{flowSet}},
#' \code{\link[flowCore:filter-methods]{filter}} for evaluation of
#' \code{boundaryFilters} and \code{\link{Subset}} for subsetting of flow
#' cytometry data sets based on that.
#' @keywords classes methods
#' @examples
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' ## Create directly. Most likely from a command line
#' boundaryFilter("FSC-H", filterId="myBoundaryFilter")
#' 
#' ## To facilitate programmatic construction we also have the following
#' bf <- boundaryFilter(filterId="myBoundaryFilter", x=c("FSC-H"))
#' 
#' ## Filtering using boundaryFilter
#' fres <- filter(dat, bf)
#' fres
#' summary(fres)
#' 
#' ## We can subset the data with the result from the filtering operation.
#' Subset(dat, fres)
#' 
#' ## A boundaryFilter on the lower margins of several channels
#' bf2 <- boundaryFilter(x=c("FSC-H", "SSC-H"), side="lower")
#' 
#' 
#' @export
setClass("boundaryFilter",
         representation=representation(tolerance="numeric", side="character"),
         contains="parameterFilter",
         prototype=list(tolerance=.Machine$double.eps, filterId="defaultBoundaryFilter",
         side="both"))

##Constructor: We allow for the following inputs:
##  tolerance is always a numeric of length 1
#' @export
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
#' Class "expressionFilter"
#' 
#' 
#' A \code{\link{filter}} holding an expression that can be evaluated to a
#' logical vector or a vector of factors.
#' 
#' 
#' The expression is evaluated in the environment of the flow cytometry values,
#' hence the parameters of a \code{\link{flowFrame}} can be accessed through
#' regular R symbols. The convenience function \code{char2ExpressionFilter}
#' exists to programmatically construct expressions.
#' 
#' @name expressionFilter-class
#' @aliases expressionFilter-class expressionFilter
#' show,expressionFilter-method char2ExpressionFilter
#' @docType class
#' @usage 
#' expressionFilter(expr, ..., filterId="defaultExpressionFilter")
#' char2ExpressionFilter(expr, ..., filterId="defaultExpressionFilter")
#' @param filterId An optional parameter that sets the \code{filterId} of this
#' \code{\link{filter}}. The object can later be identified by this name.
#' @param expr A valid R expression or a character vector that can be parsed
#' into an expression.
#' @param \dots Additional arguments that are passed to the evaluation
#' environment of the expression.
#' @return
#' 
#' Returns a \code{expressionFilter} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{concreteFilter}"}, directly.
#' 
#' Class \code{"\linkS4class{filter}"}, by class \code{concreteFilter},
#' distance 2.
#' 
#' @slot expr The expression that will be evaluated in the
#' context of the flow cytometry values.
#' @slot args An environment providing additional parameters.
#' @slot deparse A character scalar of the deparsed expression.
#' @slot filterId The identifier of the filter.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form
#' \code{new("expressionFilter", ...)}, using the
#' \code{\link{expressionFilter}} constructor or, programmatically, from a
#' character string using the \code{char2ExpressionFilter} function.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "expressionFilter")}: The workhorse used to evaluate the gate on
#'     data. This is usually not called directly by the user, but
#'     internally by calls to the \code{\link{filter}} methods. }
#'   
#'   \item{show}{\code{signature(object = "expressionFilter")}: Print
#'     information about the gate. }
#'   
#' }
#' 
#' @author F. Hahne, B. Ellis
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link{filter}} for evaluation of
#' \code{sampleFilters} and \code{\link{split}} and \code{\link{Subset}}for
#' splitting and subsetting of flow cytometry data sets based on that.
#' @keywords methods classes classes
#' @examples
#' 
#' ## Loading example data
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' 
#' #Create the filter
#' ef <- expressionFilter(`FSC-H` > 200, filterId="myExpressionFilter")
#' ef
#' 
#' ## Filtering using sampeFilters
#' fres <- filter(dat, ef)
#' fres
#' summary(fres)
#' 
#' ## The result of sample filtering is a logical subset
#' newDat <- Subset(dat, fres)
#' all(exprs(newDat)[,"FSC-H"] > 200)
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' split(dat, fres)
#' 
#' ## Programmatically construct an expression
#' dat <- dat[,-8]
#' r <- range(dat)
#' cn <- paste("`", colnames(dat), "`", sep="")
#' exp <- paste(cn, ">", r[1,], "&", cn, "<", r[2,], collapse=" & ")
#' ef2 <- char2ExpressionFilter(exp, filterId="myExpressionFilter")
#' ef2
#' fres2 <- filter(dat, ef2)
#' fres2
#' summary(fres2)
#' 
#' 
#' @export
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
#' @export
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
#' @export
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
#' Class "timeFilter"
#' 
#' 
#' Define a \code{\link{filter}} that removes stretches of unusual data
#' distribution within a single parameter over time. This can be used to
#' correct for problems during data acquisition like air bubbles or clods.
#' 
#' 
#' Clods and disturbances in the laminar flow of a FACS instrument can cause
#' temporal aberrations in the data acquisition that lead to artifactual
#' values. \code{timeFilters} try to identify such stretches of disturbance by
#' computing local variance and location estimates and to remove them from the
#' data.
#' 
#' @name timeFilter-class
#' @aliases timeFilter-class timeFilter timeFilter-class show,timeFilter-method
#' @docType class
#' @usage 
#' timeFilter(..., bandwidth=0.75, binSize, timeParameter,
#' filterId="defaultTimeFilter")
#' @param \dots The names of the parameters on which the filter is supposed to
#' work on. Names can either be given as individual arguments, or as a list or
#' a character vector.
#' @param filterId An optional parameter that sets the \code{filterId} slot of
#' this gate. The object can later be identified by this name.
#' @param bandwidth,binSize Numerics used to set the \code{bandwidth} and
#' \code{binSize} slots of the object.
#' @param timeParameter Character used to set the \code{timeParameter} slot of
#' the object.
#' @return
#' 
#' Returns a \link{timeFilter} object for use in filtering
#' \code{\link{flowFrame}}s or other flow cytometry objects.
#' @note
#' 
#' See the documentation of \code{\link[flowViz:timeLinePlot]{timeLinePlot}} in
#' the \code{\link[flowViz:flowViz-package]{flowViz}} package for details on
#' visualizing temporal problems in flow cytometry data.
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
#' @slot bandwidth Object of class \code{"numeric"}. The
#' sensitivity of the filter, i.e., the amount of local variance of
#' the signal we want to allow.
#' @slot binSize Object of class \code{"numeric"}. The size
#' of the bins used for the local variance and location
#' estimation. If \code{NULL}, a reasonable default is used when
#' evaluating the filter.
#' @slot timeParameter Object of class \code{"character"},
#' used to define the time domain parameter. If \code{NULL}, the
#' filter tries to guess the time domain from the  \code{flowFrame}.
#' @slot parameters Object of class \code{"character"},
#' describing the parameters used to filter the \code{flowFrame}.
#' @slot filterId Object of class \code{"character"},
#' referencing the filter.
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("timeFilter",
#' ...)} or using the constructor \code{timeFilter}. Using the
#' constructor is the recommended way.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{\%in\%}{\code{signature(x = "flowFrame", table =
#'                                   "timeFilter")}: The workhorse used to evaluate the filter on
#'     data. This is usually not called directly by the user. }
#'   
#'   \item{show}{\code{signature(object = "timeFilter")}: Print
#'     information about the filter. }
#'   
#' } 
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\link{flowFrame}}, \code{\link[flowCore:filter-class]{filter}} for
#' evaluation of \code{timeFilters} and \code{\link{split}} and
#' \code{\link{Subset}}for splitting and subsetting of flow cytometry data sets
#' based on that.
#' @keywords classes methods
#' @examples
#' 
#' ## Loading example data
#' data(GvHD)
#' dat <- GvHD[1:10]
#' 
#' ## create the filter
#' tf <- timeFilter("SSC-H", bandwidth=1, filterId="myTimeFilter")
#' tf
#' 
#' ## Visualize problems
#' \dontrun{
#' library(flowViz)
#' timeLinePlot(dat, "SSC-H")
#' }
#' 
#' ## Filtering using timeFilters
#' fres <- filter(dat, tf)
#' fres[[1]]
#' summary(fres[[1]])
#' summary(fres[[7]])
#' 
#' ## The result of rectangle filtering is a logical subset
#' cleanDat <- Subset(dat, fres)
#' 
#' ## Visualizing after cleaning up
#' \dontrun{
#' timeLinePlot(cleanDat, "SSC-H")
#' }
#' 
#' ## We can also split, in which case we get those events in and those
#' ## not in the gate as separate populations
#' allDat <- split(dat[[7]], fres[[7]])
#' 
#' par(mfcol=c(1,3))
#' plot(exprs(dat[[7]])[, "SSC-H"], pch=".")
#' plot(exprs(cleanDat[[7]])[, "SSC-H"], pch=".")
#' plot(exprs(allDat[[2]])[, "SSC-H"], pch=".")
#' 
#' @export 
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
#' @export
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
#' Class filterSet(deprecated)
#' 
#' A container for a collection of related filters.
#' 
#' \code{filterSet} objects are intended to provide a convenient grouping
#' mechanism for a particular gating strategy. To accomplish this, much
#' like the \code{flowSet} object, the \code{filterSet} object introduces
#' reference semantics through the use of an \code{environment}, allowing
#' users to change an upstream filter via the usual assignment mechanism
#' and have that change reflected in all dependent filters. We do this by
#' actually creating two filters for each filter in the \code{filterSet}.
#' The first is the actual concrete filter, which is assigned to a
#' variable of the form \code{.name} where \code{name} is the original
#' filter name. A second \code{filterReference} filter is the created
#' with the original name to point to the internal name. The allows us to
#' evaluate a \code{formula} in the environment without creating a copy
#' of the original filter.
#' 
#' 
#' @name filterSet-class
#' @aliases filterSet-class filterSet
#' filterReference,filterSet,character-method identifier,filterSet-method
#' identifier,filterSet,character-method names,filterSet-method
#' show,filterSet-method sort,filterSet-method %subset%,filterSet,filter-method
#' [,filterSet,character-method [[,filterSet,character-method
#' [[<-,filterSet,NULL,ANY,filter-method [[<-,filterSet,NULL,ANY,formula-method
#' [[<-,filterSet,character,ANY,formula-method
#' [[<-,filterSet,character,ANY,filter-method [[<-,filterSet,missing,ANY,filter
#' [[<-,filterSet,ANY,ANY,filterReference-method
#' [[<-,filterSet,missing,ANY,filter-method
#' @docType class
#' @section Objects from the Class:
#' 
#' There are several ways to create a \code{filterSet} object. There is the
#' \code{\link{filterSet}} constructor, which cre  ates an empty \code{filterSet}
#' object (see the details section for more information). \code{filterSet}
#' objects can also be coerced to and from \code{list} objects using the
#' \code{as} function.
#' 
#' @slot env The environment that actually holds the filters.
#' @slot name A more descriptive name of the set.
#' 
#' @section Methods:
#' \describe{
#'   \item{names}{An unsorted list of the names of the filters contained
#'     within the set.}
#'   
#'   \item{sort}{Returns a topological sort of the names of the filters contained within the 
#'     set. Primarily used by internal functions (such as \code{\link{filter}}), this method is 
#'     also useful for planning gating strategy layouts and the like.}
#'   
#'   \item{filterReference}{Retrieves references to a filter inside a
#'     filterSet}
#'   
#'   \item{[}{Returns the filter reference used inside the filter. See
#'     Details.}
#'   
#'   \item{[[}{Retrieves the actual filters from a filterSet. Note that
#'     composed filters can still contain references.}
#'   
#'   \item{[[<-}{Put a filter into a filterSet. As a convenience, assigning
#'     to the ``""'' or \code{NULL} name will use the filter's name for
#'     assignment. Composed filters can be added easily using formulas
#'     rather than attempting to construct filters the long way. The
#'     formula interface is also lazy, allowing you to add filters in any
#'     order.}
#' }
#'
#'
#' @author B. Ellis
#' @seealso
#' 
#' \code{\link{filterSet}}
#' @keywords classes
#'
#' @export
setClass("filterSet",
         representation=representation(env="environment",
         name="character"),
         prototype=prototype(env=new.env(hash=TRUE, parent=emptyenv()),
         name="Filter Set"))

## constructor
#' @export
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
#' Class filterReference
#' 
#' A reference to another filter inside a reference. Users should generally not
#' be aware that they are using this class, but it is used heavily by
#' \code{"\linkS4class{filterSet}"} classes.
#' 
#' 
#' @name filterReference-class
#' @aliases filterReference-class filterReference
#' filterReference,environment,character-method summary,filterReference-method
#' length,filterReference-method show,filterReference-method
#' eval,filterReference,missing-method
#' @docType class
#' @section Objects from the Class: Objects are generally not created by users
#' so there is no constructor function.
#' 
#' @slot name The R name of the referenced filter.
#' @slot env The environment where the filter must live.
#' @slot filterId The filterId, not really used since you always resolve.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' 
#' @author B. Ellis
#' @seealso \code{"\linkS4class{filterSet}"}
#' @keywords classes
#'
#' @export
setClass("filterReference",
         representation=representation(name="character",
         env="environment"),
         contains="filter")

## Constructor from an environment
#' @export
setMethod("filterReference",
          signature("environment", "character"),
          function(from, name) {
              new("filterReference", name=name, env=from)
          })

## Constructor from another filterSet
#' @export
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
#' Class "setOperationFilter"
#' 
#' This is a Superclass for the unionFilter, intersectFilter, complementFilter
#' and subsetFilter classes, which all consist of two or more component filters
#' and are constructed using set operators (\code{&}, \code{|}, \code{!}, and
#' \code{\%&\%} or \code{\%subset\%} respectively).
#' 
#' 
#' @name setOperationFilter-class
#' @aliases setOperationFilter-class setOperationFilter
#' @docType class
#' 
#' @slot filters Object of class \code{"list"}, containing the component filters.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' @author B. Ellis
#' @family setOperationFilter classes
#' @seealso \code{\link[flowCore:filter-methods]{filter}}
#' @keywords classes
#'
#' @export
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
#' Class unionFilter
#' 
#' This class represents the union of two filters, which is itself a filter
#' that can be incorporated in to further set operations. \code{unionFilter}s
#' are constructed using the binary set operator \code{"|"} with operands
#' consisting of a single \code{filter} or list of \code{filters}.
#' 
#' @name unionFilter-class
#' @aliases unionFilter-class unionFilter show,unionFilter-method
#' @docType class
#' 
#' @slot filters Object of class \code{"list"}, containing the component filters.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' @author B. Ellis
#' @family setOperationFilter classes
#' @seealso \code{\link[flowCore:filter-methods]{filter}}, \code{\linkS4class{setOperationFilter}}
#' @keywords classes
#' 
#' @export 
setClass("unionFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
#' @export
setMethod("|",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("unionFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "or", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
#' @export
setMethod("|",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "|", e2=e2))
#' @export
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
#' Class intersectFilter
#' 
#' This class represents the intersection of two filters, which is itself a filter
#' that can be incorporated in to further set operations. \code{intersectFilter}s
#' are constructed using the binary set operator \code{"&"} with operands consisting
#' of a single filter or list of filters.
#' 
#' @name intersectFilter-class
#' @aliases intersectFilter-class intersectFilter show,intersectFilter-method
#' @docType class
#' 
#' @slot filters Object of class \code{"list"}, containing the component filters.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' @author B. Ellis
#' @family setOperationFilter classes
#' @seealso \code{\link[flowCore:filter-methods]{filter}}, \code{\linkS4class{setOperationFilter}}
#' @keywords classes
#' 
#' @export
setClass("intersectFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
#' @export
setMethod("&",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("intersectFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "and", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
#' @export
setMethod("&",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "&", e2=e2))
#' @export
setMethod("&",
          signature=signature(e1="filter",
          e2="list"),
          definition=function(e1, e2) lapply(e2, "&", e1=e1))



## ===========================================================================
## complementFilter 
## ---------------------------------------------------------------------------
## The complement of a filters, i.e, the logical ! operation.
## ---------------------------------------------------------------------------
#' Class complementFilter
#' 
#' This class represents the logical complement of a single filter, which is 
#' itself a filter that can be incorporated in to further set operations. 
#' \code{complementFilter}s are constructed using the prefix unary set operator 
#' \code{"!"} with a single filter operand.
#' 
#' @name complementFilter-class
#' @aliases complementFilter-class complementFilter show,complementFilter-method
#' @docType class
#' 
#' @slot filters Object of class \code{"list"}, containing the component filters.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' @author B. Ellis
#' @family setOperationFilter classes
#' @seealso \code{\link[flowCore:filter-methods]{filter}}, \code{\linkS4class{setOperationFilter}}
#' @keywords classes
#' 
#' @export
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
#' @export
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
## intersection filter, the only difference is in data-driven filters.
## ---------------------------------------------------------------------------
#' Class subsetFilter
#' 
#' This class represents the action of applying a filter on the subset of
#' data resulting from another filter. This is itself a filter that can be 
#' incorporated in to further set operations. This is similar to an
#' intersectFilter, with behavior only differing if the component filters
#' are data-driven.
#' 
#' \code{subsetFilter}s are constructed using the equivalent binary set operators 
#' \code{"\%&\%"} or \code{"\%subset\%"}. The operator is not symmetric, as the
#' filter on the right-hand side will take the subset of the filter on the
#' left-hand side as input. Left-hand side operands can be a filter, list of
#' filters, or filterSet, while the right-hand side operand must be a single
#' filter.
#' 
#' @name subsetFilter-class
#' @aliases subsetFilter-class subsetFilter show,subsetFilter-method
#' summary,subsetFilter-method
#' @docType class
#' 
#' @slot filters Object of class \code{"list"}, containing the component filters.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' @author B. Ellis
#' @family setOperationFilter classes
#' @seealso \code{\link[flowCore:filter-methods]{filter}}, \code{\linkS4class{setOperationFilter}}
#' @keywords classes
#' 
#' @export
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

#' Take the intersection of two filters
#' 
#' 
#' There are two notions of intersection in \code{flowCore}. First, there is
#' the usual intersection boolean operator \code{&} that has been overridden to
#' allow the intersection of two filters or of a filter and a list for
#' convenience. There is also the \code{\%&\%} or \code{\%subset\%} operator that
#' takes an intersection, but with subset semantics rather than simple
#' intersection semantics. In other words, when taking a subset, calculations
#' from \code{\link[flowCore:filterSummary-class]{summary}} and other methods
#' are taken with respect to the right hand filter. This primarily affects
#' calculations, which are ordinarily calculated with respect to the entire
#' population as well as data-driven gating procedures which will operate only
#' on elements contained by the right hand filter.  This becomes especially
#' important when using filters such as
#' \code{\link[flowCore:norm2Filter-class]{norm2Filter}}
#' 
#' 
#' @name filter-and-methods
#' @aliases intersectFilter-method subsetFilter-method %&% %&%-methods
#' %&%,ANY-method %&%,filter,filter-method %subset%,ANY-method %subset%
#' &,filter,filter-method &,filter,list-method &,list,filter-method
#' %subset%,filter,filter-method %subset%,list,filter-method
#' coerce,intersectFilter,call-method
#' @docType methods
#' 
#' @param e1,e2 \code{\linkS4class{filter}} objects or lists of filter objects
#' 
#' @usage 
#' e1 \%&\% e2
#' e1 \%subset\% e2
#' 
#' @author B. Ellis
#' @keywords methods
## constructor from two filters. %&% is an alias for %subset%
#' @export
setMethod("%subset%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("subsetFilter",
              filters=list(e1, e2), filterId=paste(identifier(e1),"in",
                                    identifier(e2)))
      })
#' @export
setMethod("%&%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2) e1 %subset% e2)

## constructor from a list of filters and a filter
#' @export
setMethod("%subset%",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "%subset%", e2=e2))

## constructor from a filterSet and a filter
#' @export
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
#' Class "filterResult"
#' 
#' Container to store the result of applying a \code{filter} on a
#' \code{flowFrame} object
#' 
#' 
#' @name filterResult-class
#' @aliases filterResult-class filterResult ==,filterResult,flowFrame-method
#' show,filterResult-method [[,filterResult,ANY-method
#' @docType class
#' 
#' @slot frameId Object of class \code{"character"}
#' referencing the \code{flowFrame} object filtered. Used for sanity checking.
#' @slot filterDetails Object of class \code{"list"}
#' describing the filter applied.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filter}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   \item{==}{test equality}
#' }
#' 
#' @author B. Ellis, N. LeMeur
#' @seealso \code{\link[flowCore:filter-methods]{filter}},
#' \code{"\linkS4class{logicalFilterResult}"},
#' \code{"\linkS4class{multipleFilterResult}"},
#' \code{"\linkS4class{randomFilterResult}"}
#' @keywords classes
#' @examples
#' 
#' showClass("filterResult")
#' 
#' @export
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
## Slot subSet is a logical vector indicating the population membership of the
## data in the gated flowFrame.
## ---------------------------------------------------------------------------
#' Class "logicalFilterResult"
#' 
#' Container to store the result of applying a \code{filter} on a
#' \code{flowFrame} object
#' 
#' 
#' @name logicalFilterResult-class
#' @aliases logicalFilterResult-class logicalFilterResult
#' summary,logicalFilterResult-method names,logicalFilterResult-method
#' length,logicalFilterResult-method [[,logicalFilterResult,ANY-method
#' @docType class
#' 
#' @slot subSet Object of class \code{"numeric"}, which is a logical
#' vector indicating the population membership of the data in the gated
#' flowFrame.
#' @slot frameId Object of class \code{"character"}  referencing the 
#' \code{flowFrame} object filtered. Used for sanity checking.
#' @slot filterDetails Object of class \code{"list"} describing the filter 
#' applied.
#' @slot filterId Object of class \code{"character"} referencing the filter 
#' applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filterResult}"}, directly.
#' Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
#' 
#' @author B. Ellis
#' @seealso \code{\link[flowCore:filter-methods]{filter}}
#' @keywords classes
#' @examples
#' 
#' showClass("logicalFilterResult")
#' 
#' @export
setClass("logicalFilterResult",
         representation=representation(subSet="logical"),
         contains="filterResult")



## ===========================================================================
## multipleFilterResult
## ---------------------------------------------------------------------------
## Results from a filtering operation that produces multiple populations.
## Slot subSet is a factor vector indicating the population membership of the
## data in the gated flowFrame. Factor names are used as population names.
## ---------------------------------------------------------------------------
#' Class "multipleFilterResult"
#' 
#' Container to store the result of applying \code{filter} on set of
#' \code{flowFrame} objects
#' 
#' 
#' @name multipleFilterResult-class
#' @aliases multipleFilterResult-class multipleFilterResult
#' length,multipleFilterResult-method names,multipleFilterResult-method
#' names<-,multipleFilterResult-method names<-,multipleFilterResult,ANY-method
#' [[,multipleFilterResult-method [[,multipleFilterResult,ANY-method
#' [,multipleFilterResult,ANY-method summary,multipleFilterResult-method
#' show,multipleFilterResult-method
#' @docType class
#' 
#' @slot subSet Object of class \code{"factor"} indicating the population
#' membership of the data in the gated flowFrame.
#' @slot frameId Object of class \code{"character"}
#' referencing the \code{flowFrame} object filtered. Used for
#' sanity checking.
#' @slot filterDetails Object of class \code{"list"}
#' describing the filter applied.
#' @slot filterId Object of class \code{"character"}
#' referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filterResult}"}, directly.
#' Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{[, [[}{subsetting. If \code{x} is \code{multipleFilterResult},
#'     then \code{x[[i]]} a \code{FilterResult}  object. The semantics is
#'     similar to the behavior of the subsetting operators for lists.}
#'   \item{length}{number of \code{FilterResult} objects in the set.}
#'   \item{names}{names of the  \code{FilterResult} objects in the set.}
#'   \item{summary}{summary \code{FilterResult} objects in the set.}
#' }
#' 
#' @author B. Ellis
#' @seealso \code{\link[flowCore:filterResult-class]{filterResult}}
#' @keywords classes
#' @examples
#' 
#' showClass("multipleFilterResult")
#' 
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
#' Class "manyFilterResult"
#' 
#' The result of a several related, but possibly overlapping filter results.
#' The usual creator of this object will usually be a \code{\link{filter}}
#' operation of \code{\link{filterSet}} object on a \code{\link{flowFrame}}
#' object.
#' 
#' 
#' @name manyFilterResult-class
#' @aliases manyFilterResult-class length,manyFilterResult-method
#' names,manyFilterResult-method [[,manyFilterResult-method
#' [[,manyFilterResult,ANY-method summary,manyFilterResult-method
#' show,manyFilterResult-method as.data.frame.manyFilterResult manyFilterResult
#' parameters,manyFilterResult-method
#' @docType class
#' 
#' @slot subSet Object of class \code{"matrix"}.
#' @slot frameId Object of class \code{"character"} referencing the 
#' \code{flowFrame} object filtered. Used for sanity checking.
#' @slot filterDetails Object of class \code{"list"} describing the
#' filter applied.
#' @slot filterId Object of class \code{"character"} referencing the
#' filter applied.
#' @slot dependency Any dependencies between the filters. Currently
#' not used.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filterResult}"}, directly.
#' Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{[, [[}{subsetting. If \code{x} is \code{manyFilterResult},
#'     then \code{x[[i]]} a \code{filterResult}  object. The semantics is
#'     similar to the behavior of the subsetting operators for lists.}
#'   \item{length}{number of \code{filterResult} objects in the set.}
#'   \item{names}{names of the  \code{filterResult} objects in the set.}
#'   \item{summary}{summary \code{filterResult} objects in the set.}
#' }
#' 
#' @author B. Ellis
#' @seealso \code{\link[flowCore:filterResult-class]{filterResult}}
#' @keywords classes
#' @examples
#' 
#' showClass("manyFilterResult")
#' 
#' @export
setClass("manyFilterResult",
         representation=representation(subSet="matrix",
         dependency="ANY"),
         contains="filterResult")

##constructor
#' @export
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
#' Class "randomFilterResult"
#' 
#' Container to store the result of applying a \code{filter} on a
#' \code{flowFrame} object, with the population membership considered to be
#' stochastic rather than absolute. Currently not utilized.
#' 
#' @name randomFilterResult-class
#' @aliases randomFilterResult-class randomFilterResult
#' @docType class
#' 
#' @slot subSet Object of class \code{"numeric"}, which is a logical vector 
#' indicating the population membership of the data in the gated flowFrame.
#' @slot frameId Object of class \code{"character"} referencing the
#' \code{flowFrame} object filtered. Used for sanity checking.
#' @slot filterDetails Object of class \code{"list"} describing the filter applied.
#' @slot filterId Object of class \code{"character"} referencing the filter applied.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{filterResult}"}, directly.
#' Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
#' 
#' @author B. Ellis
#' @seealso \code{\link[flowCore:filter-methods]{filter}}
#' @keywords classes
#'
#' @export
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
#' Class "filterResultList"
#' 
#' Container to store the result of applying a \code{filter} on a
#' \code{flowSet} object
#' 
#' 
#' @name filterResultList-class
#' @aliases filterResultList-class filterResultList
#' [,filterResultList,ANY-method [[,filterResultList,ANY-method
#' names,filterResultList-method parameters,filterResultList-method
#' show,filterResultList-method split,flowSet,filterResultList-method
#' summary,filterResultList-method
#' @docType class
#' @section Objects from the Class:
#' 
#' Objects are created by applying a \code{\link{filter}} on a
#' \code{\link{flowSet}}. The user doesn't have to deal with manual object
#' instantiation.
#' 
#' @slot .Data Object of class \code{"list"}. The class
#' directly extends \code{list}, and this slot holds the list data.
#' @slot frameId Object of class \code{"character"} The IDs of
#' the \code{\link[flowCore:flowFrame-class]{flowFrames}} in the filtered
#' \code{\link{flowSet}}.
#' @slot filterDetails Object of class \code{"list"}. Since
#' \code{filterResultList} inherits from \code{\link{filterResult}},
#' this slot has to be set. It contains only the input filter.
#' @slot filterId Object of class \code{"character"}. The
#' identifier for the object.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{list}"}, from data part.
#' Class \code{"\linkS4class{filterResult}"}, directly.
#' Class \code{"\linkS4class{concreteFilter}"}, by class "filterResult", distance 2.
#' Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 3.
#' 
#' @section Methods:
#' \describe{
#'   \item{[}{\code{signature(x = "filterResultList", i = "ANY")}: Subset
#'     to \code{filterResultList}. }
#'   \item{[[}{\code{signature(x = "filterResultList", i = "ANY")}: Subset
#'     to individual \code{\link{filterResult}}. }
#'   
#'   \item{names}{\code{signature(x = "filterResultList")}: Accessor to
#'     the frameId slot. }
#'   
#'   \item{parameters}{\code{signature(object = "filterResultList")}:
#'       Return parameters on which data has been filtered. }
#'   
#'   \item{show}{\code{signature(object = "filterResultList")}: Print
#'     details about the object. }
#'   
#'   \item{split}{\code{signature(x = "flowSet", f =
#'                                  "filterResultList")}: Split a \code{\link{flowSet}} based on the
#'     results in the \code{filterResultlIst}. See \code{\link{split}}
#'     for details. }
#'   
#'   \item{summary}{\code{signature(object = "filterResultList")}:
#'       Summarize the filtering operation. This creates a
#'     \code{\link[flowCore:filterSummaryList-class]{filterSummaryList}}
#'     object. } 
#' }
#' 
#' @author Florian Hahne
#' @seealso \code{\linkS4class{filter}}, \code{\linkS4class{filterResult}},
#' \code{\linkS4class{logicalFilterResult}},
#' \code{\linkS4class{multipleFilterResult}},
#' \code{\linkS4class{randomFilterResult}}
#' @keywords classes
#' @examples
#' 
#' library(flowStats)
#' ## Loading example data and creating a curv1Filter
#' data(GvHD)
#' dat <- GvHD[1:3]
#' c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)
#' 
#' ## applying the filter
#' fres <- filter(dat, c1f)
#' fres
#' 
#' ## subsetting the list
#' fres[[1]]
#' fres[1:2]
#' 
#' ## details about the object
#' parameters(fres)
#' names(fres)
#' summary(fres)
#' 
#' ## splitting based on the filterResults
#' split(dat, fres)
#' 
#' @export
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
#' Class "filterList"
#' 
#' Container for a list of \code{\link[flowCore:filter-methods]{filter}}
#' objects. The class mainly exists for method dispatch.
#' 
#' 
#' @name filterList-class
#' @aliases filterList-class filterList show,filterList-method
#' identifier,filterList-method identifier<-,filterList,character-method
#' @docType class
#' @usage filterList(x, filterId=identifier(x[[1]]))
#' @param x A list of \code{\link{filter}} objects.
#' @param filterId The global identifier of the filter list. As default, we
#' take the filterId of the first \code{filter} object in \code{x}.
#' @return
#' 
#' A \code{filterList} object for the constructor.
#' @section Objects from the Class: Objects are created from regular lists
#' using the constructor \code{filterList}.
#' 
#' @slot .Data Object of class \code{"list"}. The class
#' directly extends \code{list}, and this slot holds the list data.
#' @slot filterId Object of class \code{"character"}. The
#' identifier for the object.
#' 
#' @section Extends:
#' 
#' Class \code{"\linkS4class{list}"}, from data part.
#' 
#' @section Methods:
#' 
#' \describe{
#'  \item{show}{\code{signature(object = "filterList")}: Print
#'  details about the object. }
#'  
#'  \item{identifier, identifier<-}{\code{signature(object =
#'  "filterList")}: Accessor and replacement method for the object's
#'  filterId slot. }
#'  }
#' 
#' @author Florian Hahne
#' @seealso \code{\link[flowCore:filter-methods]{filter}},
#' @keywords classes
#' @examples
#' 
#' f1 <- rectangleGate(FSC=c(100,200), filterId="testFilter")
#' f2 <- rectangleGate(FSC=c(200,400))
#' fl <- filterList(list(a=f1, b=f2))
#' fl
#' identifier(fl)
#' 
#'
#' @export
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
#' @export
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
## In the case of multipleFilterResults, the individual slots(except 'count')
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
#' Class "filterSummary"
#' 
#' Class and methods to handle the summary information of a gating operation.
#' 
#' 
#' Calling \code{summary} on a \code{\link{filterResult}} object prints summary
#' information on the screen, but also creates objects of class
#' \code{filterSummary} for computational access.
#' 
#' @name filterSummary-class
#' @aliases filterSummary-class filterSummary summary,filterResult-method
#' [[,filterSummary,numeric-method [[,filterSummary,character-method
#' $,filterSummary-method coerce,filterSummary,data.frame-method
#' length,filterSummary-method names,filterSummary-method
#' print,filterSummary-method show,filterSummary-method toTable
#' toTable,filterSummary-method
#' @docType class
#' @usage
#' \S4method{summary}{filterResult}(object, \dots)
#' @param object An object inheriting from class \code{\link{filterResult}}
#' which is to be summarized.
#' @param \dots Further arguments that are passed to the generic.
#' @return
#' 
#' An object of class \code{filterSummary} for the \code{summary} constructor,
#' a named list for the subsetting operators. The \code{$} operator returns a
#' named vector of the respective value, where each named element corresponds
#' to one sub-population.
#' @section Objects from the Class:
#' 
#' Objects are created by calling \code{summary} on a \code{link{filterResult}}
#' object. The user doesn't have to deal with manual object instantiation.
#' 
#' @slot name Object of class \code{"character"} The name(s) of
#' the populations created in the filtering operation. For a
#' \code{\link{logicalFilterResult}} this is just a single value; the
#' name of the \code{link{filter}}.
#' @slot true Object of class \code{"numeric"}. The number of
#' events within the population(s).
#' @slot count Object of class \code{"numeric"}. The total
#' number of events in the gated \code{\link{flowFrame}}.
#' @slot p Object of class \code{"numeric"} The percentage of
#' cells in the population(s).
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{[[}{\code{signature(x = "filterSummary", i = "numeric")}:
#'       Subset the \code{filterSummary} to a single population. This only
#'     makes sense for
#'     \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResults}}.
#'     The output is a list of summary statistics. }
#'   
#'   \item{[[}{\code{signature(x = "filterSummary", i = "character")}:
#'       see above }
#'   
#'   \item{$}{\code{signature(x = "filterSummary", name = "ANY")}: A
#'     list-like accessor to the slots and more. Valid values are
#'     \code{n} and \code{count} (those are identical), \code{true} and
#'     \code{in} (identical), \code{false} and \code{out} (identical),
#'     \code{name}, \code{p} and \code{q} (\code{1-p}).  }
#'   
#'   \item{coerce}{\code{signature(from = "filterSummary", to =
#'                                   "data.frame")}: Coerce object to \code{data.frame}. }
#'   
#'   \item{length}{\code{signature(x = "filterSummary")}: The number of
#'     populations in the \code{fitlerSummary}. }
#'   
#'   \item{names}{\code{signature(x = "filterSummary")}: The names of the
#'     populations in the \code{filterSummary}. }
#'   
#'   \item{print}{\code{signature(x = "filterSummary")}: Print details
#'     about the object. }
#'   
#'   \item{show}{\code{signature(object = "filterSummary")}: Print
#'     details about the object.}
#'   
#'   \item{toTable}{\code{signature(x = "filterSummary")}: Coerce object
#'     to \code{data.frame}. }
#' }
#' 
#' @author Florian Hahne, Byron Ellis
#' @seealso
#' 
#' \code{\linkS4class{filterResult}}, \code{\linkS4class{logicalFilterResult}},
#' \code{\linkS4class{multipleFilterResult}}, \code{\linkS4class{flowFrame}}
#' \code{\linkS4class{filterSummaryList}}
#' @keywords classes
#' @examples
#' 
#' library(flowStats)
#' 
#' ## Loading example data, creating and applying a curv1Filter
#' dat <- read.FCS(system.file("extdata","0877408774.B08",
#' package="flowCore"))
#' c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)
#' fres <- filter(dat, c1f)
#' 
#' ## creating and showing the summary
#' summary(fres)
#' s <- summary(fres)
#' 
#' ## subsetting
#' s[[1]]
#' s[["peak 2"]]
#' 
#' ##accessing details
#' s$true
#' s$n
#' toTable(s)
#' 
#' 
#' @export
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
#' Class "filterSummaryList"
#' 
#' 
#' Class and methods to handle summary statistics for from filtering operations
#' on whole \code{\link[flowCore:flowSet-class]{flowSets}}.
#' 
#' 
#' Calling \code{summary} on a \code{\link{filterResultList}} object prints summary
#' information on the screen, but also creates objects of class
#' \code{filterSummaryList} for computational access.
#' 
#' @name filterSummaryList-class
#' @aliases filterSummaryList-class filterSummaryList
#' toTable,filterSummaryList-method
#' @docType class
#' @section Usage:
#' summary(object, \dots)
#' @param object An object of class.
#' \code{\link[flowCore:filterResultList-class]{filterResultList}} which is to
#' be summarized.
#' @param \dots Further arguments that are passed to the generic.
#' @return
#' 
#' An object of class \code{filterSummaryList}.
#' @section Objects from the Class:
#' 
#' Objects are created by calling \code{summary} on a
#' \code{link{filterResultList}} object. The user doesn't have to deal with
#' manual object instantiation.
#' 
#' @slot .Data Object of class \code{"list"}. The class
#' directly extends \code{list}, and this slot holds the list data.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{list}"}, from \code{.Data} part.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{toTable}{\code{signature(x = "filterSummaryList")}: Coerce
#'     object to \code{data.frame}. Additional factors are added to
#'     indicate list items in the original object. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{filterResult}}, \code{\linkS4class{filterResultList}},
#' \code{\linkS4class{logicalFilterResult}},
#' \code{\linkS4class{multipleFilterResult}}, \code{\linkS4class{flowFrame}}
#' \code{\linkS4class{filterSummary}}
#' @keywords classes
#' @examples
#' 
#' library(flowStats)
#' 
#' ## Loading example data, creating and applying a curv1Filter
#' data(GvHD)
#' dat <- GvHD[1:3]
#' c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)
#' fres <- filter(dat, c1f)
#' 
#' ## creating and showing the summary
#' summary(fres)
#' s <- summary(fres)
#' 
#' ## subsetting
#' s[[1]]
#' 
#' ##accessing details
#' toTable(s)
#' 
#' 
#' @export
setClass("filterSummaryList",
         contains="list")



## ===========================================================================
## transform functions
## ---------------------------------------------------------------------------
## Constructors for the different varieties of transforms. All of these
## create objects of the basic class 'transform', unless stated otherwise.
## ---------------------------------------------------------------------------
#' Create the definition of a linear transformation function to be applied on a
#' data set
#' 
#' Create the definition of the linear Transformation that will be applied on
#' some parameter via the \code{transform} method.  The definition of this
#' function is currently x <- a*x+b
#' 
#' @usage linearTransform(transformationId="defaultLinearTransform", a = 1, b = 0)
#' @param transformationId character string to identify the transformation
#' @param a double that corresponds to the multiplicative factor in the
#' equation
#' @param b double that corresponds to the additive factor in the equation
#' @return Returns an object of class \code{transform}.
#' @author N. LeMeur
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   linearTrans <- linearTransform(transformationId="Linear-transformation", a=2, b=0)
#'   dataTransform <- transform(samp, transformList('FSC-H' ,linearTrans))
#' 
#' 
#' @export
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
#' Create the definition of a quadratic transformation function to be applied
#' on a data set
#' 
#' Create the definition of the quadratic Transformation that will be applied
#' on some parameter via the \code{transform} method.  The definition of this
#' function is currently x <- a*x\^2 + b*x + c
#' 
#' @usage quadraticTransform(transformationId="defaultQuadraticTransform", a = 1, b = 1, c = 0)
#' @param transformationId character string to identify the transformation
#' @param a double that corresponds to the quadratic coefficient in the
#' equation
#' @param b double that corresponds to the linear coefficient in the equation
#' @param c double that corresponds to the intercept in the equation
#' @return Returns an object of class \code{transform}.
#' @author N. Le Meur
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   quadTrans <- quadraticTransform(transformationId="Quadratic-transformation", a=1, b=1, c=0)
#'   dataTransform <- transform(samp, transformList('FSC-H', quadTrans))
#' 
#' 
#' @export
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
#' Create the definition of a ln transformation function (natural logarthim) to
#' be applied on a data set
#' 
#' Create the definition of the ln Transformation that will be applied on some
#' parameter via the \code{transform} method.  The definition of this function
#' is currently x<-log(x)*(r/d).  The transformation would normally be used to
#' convert to a linear valued parameter to the natural logarithm scale.
#' Typically r and d are both equal to 1.0. Both must be positive.
#' 
#' @usage lnTransform(transformationId="defaultLnTransform", r=1, d=1)
#' @param transformationId character string to identify the transformation
#' @param r positive double that corresponds to a scale factor.
#' @param d positive double that corresponds to a scale factor
#' @return Returns an object of class \code{transform}.
#' @author B. Ellis and N. LeMeur
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#'   data(GvHD)
#'   lnTrans <- lnTransform(transformationId="ln-transformation", r=1, d=1)
#'   ln1 <- transform(GvHD, transformList('FSC-H', lnTrans))
#' 
#' opar = par(mfcol=c(2, 1))
#' plot(density(exprs(GvHD[[1]])[ ,1]), main="Original")
#' plot(density(exprs(ln1[[1]])[ ,1]), main="Ln Transform")
#' 
#' 
#' @export
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
#' Create the definition of a log transformation function (base specified by
#' user) to be applied on a data set
#' 
#' Create the definition of the log Transformation that will be applied on some
#' parameter via the \code{transform} method.  The definition of this function
#' is currently x<-log(x,logbase)*(r/d).  The transformation would normally be
#' used to convert to a linear valued parameter to the natural logarithm scale.
#' Typically r and d are both equal to 1.0. Both must be positive.  logbase =
#' 10 corresponds to base 10 logarithm.
#' 
#' @usage logTransform(transformationId="defaultLogTransform", logbase=10, r=1, d=1)
#' @param transformationId character string to identify the transformation
#' @param logbase positive double that corresponds to the base of the
#' logarithm.
#' @param r positive double that corresponds to a scale factor.
#' @param d positive double that corresponds to a scale factor
#' @return Returns an object of class \code{transform}.
#' @author B. Ellis, N. LeMeur
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   logTrans <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
#'   trans <- transformList('FSC-H', logTrans)
#'   dataTransform <- transform(samp, trans)
#' 
#' @export
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
#' Compute a transform using the 'biexponential' function
#' 
#' The 'biexponential' is an over-parameterized inverse of the hyperbolic sine.
#' The function to be inverted takes the form biexp(x) =
#' a*exp(b*(x-w))-c*exp(-d*(x-w))+f with default parameters selected to
#' correspond to the hyperbolic sine.
#' 
#' @usage
#' biexponentialTransform(transformationId="defaultBiexponentialTransform", 
#'                        a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0, 
#'                        tol = .Machine$double.eps^0.25, maxit = as.integer(5000))
#' @param transformationId A name to assign to the transformation. Used by the
#' transform/filter integration routines.
#' @param a See the function description above. Defaults to 0.5
#' @param b See the function description above. Defaults to 1.0
#' @param c See the function description above. Defaults to 0.5 (the same as
#' \code{a})
#' @param d See the function description above. Defaults to 1 (the same as
#' \code{b})
#' @param f A constant bias for the intercept. Defaults to 0.
#' @param w A constant bias for the 0 point of the data. Defaults to 0.
#' @param tol A tolerance to pass to the inversion routine
#' (\code{\link{uniroot}} usually)
#' @param maxit A maximum number of iterations to use, also passed to
#' \code{\link{uniroot}}
#' @return Returns values giving the inverse of the biexponential within a
#' certain tolerance. This function should be used with care as numerical
#' inversion routines often have problems with the inversion process due to the
#' large range of values that are essentially 0. Do not be surprised if you end
#' up with population splitting about \code{w} and other odd artifacts.
#' @author B. Ellis, N Gopalakrishnan
#' @family Transform functions
#' @seealso \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' # Construct some "flow-like" data which tends to be hetereoscedastic.
#' data(GvHD)
#' biexp  <- biexponentialTransform("myTransform")
#' 
#' after.1 <- transform(GvHD, transformList('FSC-H', biexp))
#' 
#' biexp  <- biexponentialTransform("myTransform",w=10)
#' after.2 <- transform(GvHD, transformList('FSC-H', biexp))
#' 
#' opar = par(mfcol=c(3, 1))
#' plot(density(exprs(GvHD[[1]])[, 1]), main="Original")
#' plot(density(exprs(after.1[[1]])[, 1]), main="Standard Transform")
#' plot(density(exprs(after.2[[1]])[, 1]), main="Shifted Zero Point")
#'
#' @export
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
#' Computes a transform using the 'logicle_transform' function
#' 
#' 
#' Logicle transformation creates a subset of
#' \code{\link{biexponentialTransform}} hyperbolic sine transformation
#' functions that provides several advantages over linear/log transformations
#' for display of flow cytometry data. (The logicleTransform method makes use
#' of the C++ implementation of the logicle transform contributed by Wayne
#' Moore et al.)
#' 
#' 
#' @aliases logicleTransform estimateLogicle
#' @usage 
#' logicleTransform(transformationId="defaultLogicleTransform", w = 0.5, t = 262144,
#'                  m = 4.5, a = 0)
#'                  estimateLogicle(x, channels,...) 
#' @param transformationId A name to assign to the transformation. Used by the
#' transform/filter routines.
#' @param w w is the linearization width in asymptotic decades. w should be > 0
#' and determines the slope of transformation at zero.  w can be estimated
#' using the equation w=(m-log10(t/abs(r)))/2, where r is the most negative
#' value to be included in the display
#' @param t Top of the scale data value, e.g, 10000 for common 4 decade data or
#' 262144 for a 18 bit data range. t should be greater than zero
#' @param m m is the full width of the transformed display in asymptotic
#' decades. m should be greater than zero
#' @param a Additional negative range to be included in the display in
#' asymptotic decades. Positive values of the argument brings additional
#' negative input values into the transformed display viewing area. Default
#' value is zero corresponding to a Standard logicle function.
#' @param x Input flow frame for which the logicle transformations are to be
#' estimated.
#' @param channels channels or markers for which the logicle transformation is
#' to be estimated.
#' @param ... other arguments:
#' 
#' q: a numeric type specifying quantile value, default is 0.05
#' @author Wayne Moore, N Gopalakrishnan
#' @family Transform functions
#' @seealso \code{\link[flowCore]{inverseLogicleTransform}},
#' \code{\link[flowCore]{estimateLogicle} }
#' @references Parks D.R., Roederer M., Moore W.A.(2006) A new "logicle"
#' display method avoids deceptive effects of logarithmic scaling for low
#' signals and compensated data. CytometryA, 96(6):541-51.
#' @keywords methods
#' @examples
#' 
#' data(GvHD)
#' samp <- GvHD[[1]] 
#' ## User defined logicle function
#' lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)
#' trans <- transformList(c("FL1-H", "FL2-H"), lgcl)
#' after <- transform(samp, trans)
#' invLgcl <- inverseLogicleTransform(trans = lgcl)
#' trans <- transformList(c("FL1-H", "FL2-H"), invLgcl)
#' before <- transform (after,  trans)
#' 
#' ## Automatically estimate the logicle transformation based on the data
#' lgcl <- estimateLogicle(samp, channels = c("FL1-H", "FL2-H", "FL3-H", "FL2-A", "FL4-H"))
#' ## transform  parameters using the estimated logicle transformation
#' after <- transform(samp, lgcl)
#' 
#' 
#' @export
logicleTransform <- function(transformationId="defaultLogicleTransform", 
        w = 0.5, t = 262144, m = 4.5, a = 0) {

    k <- new("transform", .Data=function(x) 
            x <- logicle_transform(as.double(x), as.double(t),as.double(w), as.double(m), as.double(a), FALSE)
            )            
    k@transformationId <- transformationId
    k
}

### Inverse logicle transformation constructor
#' Computes the inverse of the transform defined by the 'logicleTransform'
#' function or the transformList generated by 'estimateLogicle' function
#' 
#' inverseLogicleTransform can be use to compute the inverse of the Logicle
#' transformation. The parameters w, t, m, a for calculating the inverse are
#' obtained from the 'trans' input passed to the 'inverseLogicleTransform'
#' function. (The inverseLogicleTransform method makes use of the C++
#' implementation of the inverse logicle transform contributed by Wayne Moore
#' et al.)
#' 
#' @usage inverseLogicleTransform(trans,transformationId,...)
#' @param trans An object of class 'transform' created using the
#' 'logicleTransform' function or class 'transformList' created by
#' 'estimateLogicle'.  The parameters w, t, m, a for calculating the inverse
#' are obtained from the 'trans' input passed to the 'inverseLogicleTransform'
#' function.
#' @param transformationId A name to assigned to the inverse transformation.
#' Used by the transform routines.
#' @param ...  not used.
#' @author Wayne Moore, N. Gopalakrishnan
#' @family Transform functions
#' @seealso \code{\link[flowCore]{logicleTransform}}
#' @references Parks D.R., Roederer M., Moore W.A.(2006) A new "logicle"
#' display method avoids deceptive effects of logarithmic scaling for low
#' signals and compensated data. CytometryA, 96(6):541-51.
#' @keywords methods
#' @examples
#' 
#' data(GvHD)
#' samp <- GvHD[[1]] 
#' 
#' #########inverse the transform object###############
#' logicle  <- logicleTransform(t = 10000, w = 0.5, m = 4.5 , a =0 ,"logicle")
#' ## transform FL1-H parameter using logicle transformation
#' after <- transform(samp, transformList('FL1-H', logicle))
#' 
#' ## Inverse transform the logicle transformed data to retrieve the original data
#' invLogicle <- inverseLogicleTransform(trans = logicle)
#' before <- transform (after, transformList('FL1-H', invLogicle))
#' 
#' #########inverse the transformList object###############
#' translist <- estimateLogicle(samp, c("FL1-H", "FL2-H"))
#' after <- transform(samp, translist)
#' ## Inverse 
#' invLogicle <- inverseLogicleTransform(translist)
#' before <- transform (after, invLogicle)
#' 
#' @export
inverseLogicleTransform <- function(trans, transformationId, ...)UseMethod("inverseLogicleTransform")
#' @export
inverseLogicleTransform.default <- function(trans, transformationId, ...) {
  
    stop("trans has to be an object of class \"transform\"
            created using the \"logicleTransform\" function\n
         or a 'transformList' created by 'estimateLogicle'\n")
}
#' @export
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
#' @export
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

#' @export
estimateLogicle <- function(x, channels, ...)UseMethod("estimateLogicle")
#' @export
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
#' Create the definition of a truncate transformation function to be applied on
#' a data set
#' 
#' Create the definition of the truncate Transformation that will be applied on
#' some parameter via the \code{transform} method.  The definition of this
#' function is currently x[x<a] <- a.  Hence, all values less than a are
#' replaced by a. The typical use would be to replace all values less than 1 by
#' 1.
#' 
#' @usage truncateTransform(transformationId="defaultTruncateTransform", a=1)
#' @param transformationId character string to identify the transformation
#' @param a double that corresponds to the value at which to truncate
#' @return Returns an object of class \code{transform}.
#' @author P. Haaland
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   truncateTrans <- truncateTransform(transformationId="Truncate-transformation", a=5)
#'   dataTransform <- transform(samp,transformList('FSC-H', truncateTrans))
#' 
#' 
#' @export
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
#' Create the definition of a scale transformation function to be applied on a
#' data set
#' 
#' Create the definition of the scale Transformation that will be applied on
#' some parameter via the \code{transform} method.  The definition of this
#' function is currently x = (x-a)/(b-a).  The transformation would normally be
#' used to convert to a 0-1 scale. In this case, b would be the maximum
#' possible value and a would be the minimum possible value.
#' 
#' @usage scaleTransform(transformationId="defaultScaleTransform", a, b)
#' @param transformationId character string to identify the transformation
#' @param a double that corresponds to the value that will be transformed to 0
#' @param b double that corresponds to the value that will be transformed to 1
#' @return Returns an object of class \code{transform}.
#' @author P. Haaland
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   scaleTrans <- scaleTransform(transformationId="Truncate-transformation", a=1, b=10^4)
#'   dataTransform <- transform(samp, transformList('FSC-H', scaleTrans))
#' 
#' @export
scaleTransform <- function(transformationId="defaultScaleTransform",
                           a=1, b=10^4)
{
    t <- new("transform", .Data=function(x) (x-a)/(b-a))
    t@transformationId <- transformationId
    t
}

## Split-scale transformation constructor
#' Compute the split-scale transformation describe by FL. Battye
#' 
#' The split scale transformation described by Francis L. Battye [B15] (Figure
#' 13) consists of a logarithmic scale at high values and a linear scale at low
#' values with a fixed transition point chosen so that the slope (first
#' derivative) of the transform is continuous at that point. The scale extends
#' to the negative of the transition value that is reached at the bottom of the
#' display.
#' 
#' @usage 
#' splitScaleTransform(transformationId="defaultSplitscaleTransform",
#'                     maxValue=1023, transitionChannel=64, r=192)
#' @param transformationId A name to assign to the transformation. Used by the
#' transform/filter integration routines.
#' @param maxValue Maximum value the transformation is applied to, e.g., 1023
#' @param transitionChannel Where to split the linear versus the logarithmic
#' transformation, e.g., 64
#' @param r Range of the logarithm part of the display, ie. it may be expressed
#' as the maxChannel - transitionChannel considering the maxChannel as the
#' maximum value to be obtained after the transformation.
#' @return Returns values giving the inverse of the biexponential within a
#' certain tolerance. This function should be used with care as numerical
#' inversion routines often have problems with the inversion process due to the
#' large range of values that are essentially 0. Do not be surprised if you end
#' up with population splitting about \code{w} and other odd artifacts.
#' @author N. LeMeur
#' @family Transform functions
#' @seealso \code{\link{transform}}
#' @references Battye F.L. A Mathematically Simple Alternative to the
#' Logarithmic Transform for Flow Cytometric Fluorescence Data Displays.
#' http://www.wehi.edu.au/cytometry/Abstracts/AFCG05B.html.
#' @keywords methods
#' @examples
#' 
#' data(GvHD)
#' ssTransform  <- splitScaleTransform("mySplitTransform")
#' after.1 <- transform(GvHD, transformList('FSC-H', ssTransform))
#' 
#' opar = par(mfcol=c(2, 1))
#' plot(density(exprs(GvHD[[1]])[, 1]), main="Original")
#' plot(density(exprs(after.1[[1]])[, 1]), main="Split-scale Transform")
#' 
#' @export
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
#' Create the definition of an arcsinh transformation function (base specified
#' by user) to be applied on a data set
#' 
#' Create the definition of the arcsinh Transformation that will be applied on
#' some parameter via the \code{transform} method.  The definition of this
#' function is currently x<-asinh(a+b*x)+c).  The transformation would normally
#' be used to convert to a linear valued parameter to the natural logarithm
#' scale. By default a and b are both equal to 1 and c to 0.
#' 
#' @usage
#' arcsinhTransform(transformationId="defaultArcsinhTransform", a=1, b=1, c=0)
#' @param transformationId character string to identify the transformation
#' @param a positive double that corresponds to a shift about 0.
#' @param b positive double that corresponds to a scale factor.
#' @param c positive double
#' @return Returns an object of class \code{transform}.
#' @author B. Ellis
#' @family Transform functions
#' @seealso \code{\link{transform-class}}, \code{\link{transform}},
#' \code{asinh}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#'   asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
#'   translist <- transformList('FSC-H', asinhTrans) 
#'   dataTransform <- transform(samp, translist)
#' 
#' @export
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
#' Class "parameterTransform"
#' 
#' Link a transformation to particular flow parameters
#' 
#' 
#' @name parameterTransform-class
#' @aliases parameterTransform-class parameterTransform
#' @docType class
#' 
#' @slot .Data Object of class \code{"function"}, the
#' transformation function.
#' @slot parameters Object of class \code{"character"} The
#' parameters the transformation is applied to.
#' @slot transformationId Object of class
#' \code{"character"}. The identifier for the object.
#' 
#' @section Objects from the Class:
#' 
#' Objects are created by using the \code{\%on\%} operator and are usually not
#' directly instantiated by the user.
#' @section Extends:
#' 
#' Class \code{"\linkS4class{transform}"}, directly.
#' Class \code{"\linkS4class{function}"}, by class "transform", distance 2.
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{\%on\%}{\code{signature(e1 = "filter", e2 =
#'                                   "parameterTransform")}: Apply the transformation. }
#'   \item{\%on\%}{\code{signature(e1 = "parameterTransform", e2 =
#'                                   "flowFrame")}: see above }
#'   \item{parameters}{\code{signature(object = "parameterTransform")}:
#'       Accessor to the parameters slot }
#' }
#' 
#' @author Byron Ellis
#' @keywords classes
#'
#' @export
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
#' A class for mapping transforms between parameters
#' 
#' 
#' This class provides a mapping between parameters and transformed parameters
#' via a function.
#' 
#' 
#' @name transformMap-class
#' @aliases transformMap-class transformMap show,transformMap-method
#' @docType class
#' 
#' @slot output Name of the transformed parameter.
#' @slot input Name of the parameter to transform.
#' @slot f Function used to accomplish the transform.
#' 
#' @section Objects from the Class:
#' 
#' Objects of this type are not usually created by the user, except perhaps in
#' special circumstances. They are generally automatically created by the
#' inline \code{\link[flowCore:transform-class]{transform}} process during the
#' creation of a \code{\link{transformFilter}}, or by a call to the
#' \code{\link{transformList}} constructor.
#' 
#' @section Methods:
#' \describe{
#'   \item{show}{\code{signature(object = "transformList")}: Print details
#'     about the object. }
#' }
#' 
#' @author B. Ellis, F. Hahne
#' @seealso
#' 
#' \code{\link{transform}}, \code{\link{transformList}}
#' @keywords classes
#' @examples
#' 
#' new("transformMap", input="FSC-H", output="FSC-H", f=log)
#' 
#' 
#' @export 
setClass("transformMap",
         representation=representation(output="character",
         input="character",
         f="function"))



## ===========================================================================
## transformList
## ---------------------------------------------------------------------------
## A list of transformMaps
## ---------------------------------------------------------------------------
#' Class "transformList"
#' 
#' A list of transformMaps to be applied to a list of parameters.
#' 
#' 
#' @name transformList-class
#' @aliases transformList-class transformList colnames,transformList-method
#' c,transformList-method identifier,transformList-method
#' identifier<-,transformList,character-method
#' @docType class
#' @usage transformList(from, tfun, to=from, transformationId =
#' "defaultTransformation")
#' 
#' @param from,to Characters giving the names of the measurement parameter on
#' which to transform on and into which the result is supposed to be stored. If
#' both are equal, the existing parameters will be overwritten.
#' @param tfun A list if functions or a character vector of the names of the
#' functions used to transform the data. R's recycling rules apply, so a single
#' function can be given to be used on all parameters.
#' @param transformationId The identifier for the object.
#' 
#' @slot transforms Object of class \code{"list"}, where each
#' list item is of class \code{\link{transformMap}}.
#' @slot transformationId Object of class \code{"character"},
#' the identifier for the object.
#' 
#' @section Objects from the Class:
#' 
#' Objects can be created by calls of the form \code{new("transformList",
#' ...)}, by calling the \code{\link{transform}} method with key-value pair
#' arguments of the form \code{key} equals character and \code{value} equals
#' function, or by using the constructor \code{transformList}. See below for
#' details
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{colnames}{\code{signature(x = "transformList")}: This returns
#'     the names of the parameters that are to be transformed. }
#'   
#'   \item{c}{\code{signature(x = "transformList")}: Concatenate
#'     \code{transformList}s or regular lists and \code{transformLists}. }
#'   
#'   \item{\%on\%}{\code{signature(e1 = "transformList", e2 =
#'                                   "flowFrame")}: Perform a transformation using the
#'     \code{transformList} on a \code{\link{flowFrame}} or
#'     \code{\link{flowSet}}. }
#' }
#' 
#' @author B. Ellis, F. Hahne
#' @seealso \code{\link{transform}}, \code{\link{transformMap}}
#' @keywords classes
#' @examples
#' 
#' tl <- transformList(c("FSC-H", "SSC-H"), list(log, asinh))
#' colnames(tl)
#' c(tl, transformList("FL1-H", "linearTransform"))
#' data(GvHD)
#' transform(GvHD[[1]], tl)
#' 
#' 
#' @export 
setClass("transformList",
         representation=representation(transforms="list",
                                       transformationId="character"),
         prototype=prototype(transformationId="defaultTransformation"),
         validity=function(object)
         if(all(sapply(object@transforms, is, "transformMap"))) TRUE else
         stop("All list items of a 'transformList' must be of class ",
              "'transformMap.'", call.=FALSE))

## constructor
#' @export
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
#' 
#' A class for encapsulating a filter to be performed on transformed parameters
#' 
#' 
#' The \code{transformFilter} class is a mechanism for including one or more
#' variable transformations into the filtering process. Using a special case of
#' \code{\link[flowCore:transform-class]{transform}} we can introduce
#' transformations inline with the filtering process eliminating the need to
#' process \code{\link[flowCore:flowFrame-class]{flowFrame}} objects before
#' applying a filter.
#' 
#' 
#' @name transformFilter-class
#' @aliases transformFilter-class transformFilter show,transformFilter-method
#' @docType class
#' 
#' @slot transforms A list of transforms to perform on the
#' target \code{\link[flowCore:flowFrame-class]{flowFrame}}
#' @slot filter The filter to be applied to the transformed
#' frame
#' @slot filterId The name of the filter (chosen
#' automatically)
#' 
#' @section Objects from the Class:
#' 
#' Objects of this type are not generally created ``by hand''. They are a side
#' effect of the use of the \code{\link[flowCore:filter-on-methods]{\%on\%}}
#' method with a \code{\link[flowCore:filter-methods]{filter}} object on the
#' left hand side and a
#' \code{\link[flowCore:transformList-class]{transformList}} on the right hand
#' side.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{filter}"}, directly.
#' 
#' @author B. Ellis
#' @seealso
#' 
#' \code{"\linkS4class{filter}"}, \code{"\linkS4class{transform}"},
#' \code{\link[flowCore:transform-class]{transform}}
#' @keywords classes
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata", "0877408774.B08", package="flowCore"))
#' 
#' ## Gate this object after log transforming the forward and side
#' ## scatter variables
#' filter(samp, norm2Filter("FSC-H", "SSC-H", scale.factor=2)
#'        %on% transform("FSC-H"=log,"SSC-H"=log))
#' 
#' 
#' @export 
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
#' Class "compensation"
#' 
#' 
#' Class and methods to compensate for spillover between channels by applying a
#' spillover matrix to a \code{flowSet} or a \code{flowFrame} assuming a simple
#' linear combination of values.
#' 
#' 
#' The essential premise of compensation is that some fluorochromes may
#' register signals in detectors that do not correspond to their primary
#' detector (usually a photomultiplier tube). To compensate for this fact, some
#' sort of standard is used to obtain the background signal (no dye) and the
#' amount of signal on secondary channels for each fluorochrome relative to the
#' signal on their primary channel.
#' 
#' To calculate the spillover percentage we use either the mean or the median
#' (more often the latter) of the secondary signal minus the background signal
#' for each dye to obtain \code{n} by \code{n} matrix, \code{S}, of so-called
#' spillover values, expressed as a percentage of the primary channel. The
#' observed values are then considered to be a linear combination of the true
#' fluorescence and the spillover from each other channel so we can obtain the
#' true values by simply multiplying by the inverse of the spillover matrix.
#' 
#' The spillover matrix can be obtained through several means. Some flow
#' cytometers provide a spillover matrix calculated during acquisition,
#' possibly by the operator, that is made available in the metadata of the
#' flowFrame.  While there is a theoretical standard keyword \code{$SPILL} it
#' can also be found in the \code{SPILLOVER} or \code{SPILL} keyword depending
#' on the cytometry. More commonly the spillover matrix is calculated using a
#' series of compensation cells or beads collected before the experiment. If
#' you have set of FCS files with one file per fluorochrome as well as an
#' unstained FCS file you can use the
#' \code{\link[flowCore:spillover-flowSet]{spillover}} method for
#' \code{\link[flowCore:flowSet-class]{flowSets}} to automatically calculate a
#' spillover matrix.
#' 
#' The \code{compensation} class is essentially a wrapper around a
#' \code{matrix} that allows for transformed parameters and method dispatch.
#' 
#' @name compensation-class
#' @aliases compensation-class compensation identifier,compensation-method
#' parameters,compensation-method identifier<-,compensation,character-method
#' show,compensation-method compensate
#' @docType class
#' @usage
#' compensation(\dots, spillover, compensationId="defaultCompensation")
#' 
#' compensate(x, spillover, \dots)
#' @param spillover The spillover or compensation matrix.
#' @param compensationId The identifier for the compensation object.
#' @param x An object of class \code{\linkS4class{flowFrame}} or
#' \code{\linkS4class{flowSet}}.
#' @param \dots Further arguments.
#' 
#' The constructor is designed to be useful in both programmatic and
#' interactive settings, and \dots{} serves as a container for possible
#' arguments. The following combinations of values are allowed:
#' 
#' Elements in \dots{} are \code{character} scalars of parameter names or
#' \code{\linkS4class{transform}} objects and the colnames in \code{spillover}
#' match to these parameter names.
#' 
#' The first element in \dots{} is a \code{character} vector of parameter names
#' or a list of \code{character} scalars or \code{\linkS4class{transform}}
#' objects and the colnames in \code{spillover} match to these parameter names.
#' 
#' Argument \code{spillover} is missing and the first element in \dots{} is a
#' \code{matrix}, in which case it is assumed to be the spillover matrix.
#' 
#' \dots{} is missing, in which case all parameter names are taken from the
#' colnames of \code{spillover}.
#' 
#' @return
#' 
#' A \code{compensation} object for the constructor.
#' 
#' A \code{\linkS4class{flowFrame}} or \code{\linkS4class{flowSet}} for the
#' \code{compensate} methods.
#' @section Objects from the Class:
#' 
#' Objects should be created using the constructor \code{compensation()}. See
#' the \code{Usage} and \code{Arguments} sections for details.
#' 
#' @slot spillover Object of class \code{matrix}; the
#' spillover matrix.
#' @slot compensationId Object of class \code{character}. An
#' identifier for the object.
#' @slot parameters Object of class \code{parameters}. The
#' flow parameters for which the compensation is defined. This can
#' also be objects of class \code{\linkS4class{transform}}, in which
#' case the compensation is performed on the compensated parameters.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{compensate}{\code{signature(x = "flowFrame", spillover =
#'                                       "compensation")}: Apply the compensation defined in a
#'     \code{compensation} object on a \code{\linkS4class{flowFrame}}.
#'     This returns a compensated \code{flowFrame}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   compensate(flowFrame, compensation)}
#'     
#'   }
#'   
#'   \item{compensate}{\code{signature(x = "flowFrame", spillover =
#'                                       "matrix")}: Apply a compensation matrix to a
#'     \code{\linkS4class{flowFrame}}.  This returns a compensated
#'     \code{flowFrame}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   compensate(flowFrame, matrix)}
#'     
#'   }
#'   
#'   \item{compensate}{\code{signature(x = "flowFrame", spillover =
#'                                       "data.frame")}:Try to coerce the \code{data.frame} to a
#'     \code{matrix} and apply that to a
#'     \code{\linkS4class{flowFrame}}.  This returns a compensated
#'     \code{flowFrame}.
#'     
#'     \emph{Usage:}
#'     
#'     \code{   compensate(flowFrame, data.frame)}
#'     
#'   }
#'   
#'   \item{identifier, identifier<-}{\code{signature(object =
#'                                                     "compensation")}: Accessor and replacement methods for the
#'     \code{compensationId} slot.
#'     
#'     \emph{Usage:}
#'     
#'     
#'     \code{   identifier(compensation)}
#'     
#'     \code{   identifier(compensation) <- value}
#'     
#'   }
#'   
#'   
#'   \item{parameters}{\code{signature(object =
#'                                       "compensation")}: Get the parameter names of the
#'     \code{compensation} object. This method also tries to resolve
#'     all \code{\link[flowCore:transform-class]{transforms}} and
#'     \code{\link[flowCore:transformReference-class]{transformReferences}}
#'     before returning the parameters as character vectors. Unresolvable
#'     references return \code{NA}.
#'     
#'     \emph{Usage:}
#'     
#'     
#'     \code{   parameters(compensation)}
#'     
#'     
#'   }
#'   
#'   
#'   \item{show}{\code{signature(object = "compensation")}: Print details
#'     about the object.
#'     
#'     \emph{Usage:}
#'     
#'     This method is automatically called when the object is printed on
#'     the screen.
#'     
#'   }  
#' }
#' 
#' @author F.Hahne, B. Ellis, N. Le Meur
#' @seealso
#' 
#' \code{\link[flowCore:spillover-flowSet]{spillover}}
#' @keywords methods classes
#' @examples
#' 
#' ## Read sample data and a sample spillover matrix
#' samp   <- read.flowSet(path=system.file("extdata", "compdata", "data",
#'           package="flowCore")) 
#' cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
#' comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
#' comp.mat
#' 
#' ## compensate using the spillover matrix directly
#' summary(samp)
#' samp <- compensate(samp, comp.mat)
#' summary(samp)
#' 
#' ## create a compensation object and compensate using that
#' comp <- compensation(comp.mat)
#' compensate(samp, comp)
#' 
#' ## demo the sample-specific compensation
#' ## create a list of comps (each element could be a 
#' ## different compensation tailored for the specific sample)
#' comps <- sapply(sampleNames(samp), function(sn)comp, simplify = FALSE)
#' # the names of comps must be matched to sample names of the flowset
#' compensate(samp, comps)
#' 
#' @export
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
#' @export
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
#' Class "compensatedParameter"
#' 
#' 
#' Emission spectral overlap can be corrected by subtracting the amount of
#' spectral overlap from the total detected signals. This compensation process
#' can be described by using spillover matrices.
#' 
#' The compensatedParameter class allows for compensation of specific parameters
#' the user is interested in by creating compensatedParameter objects and
#' evaluating them. This allows for use of compensatedParameter in gate
#' definitions.
#' 
#' 
#' @name compensatedParameter-class
#' @aliases compensatedParameter-class compensatedParameter
#' eval,compensatedParameter,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument. The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class:
#' 
#' Objects can be created by calls to the constructor of the form
#' \code{compensatedParameter(parameters,spillRefId,transformationId,searchEnv)}.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot parameters Object of class \code{"character"} -- the flow
#' parameters to be  compensated.
#' @slot spillRefId Object of class \code{"character"} -- the name of the
#' compensation object (The compensation object contains the spillover Matrix).
#' @slot searchEnv Object of class \code{"environment"} -environment in
#' which the compensation object is defined.
#' @slot transformationId Object of class \code{"character"} -- a unique Id to
#' reference the compensatedParameter object.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author Gopalakrishnan N,F.Hahne
#' @seealso compensation
#' @keywords classes
#' @examples
#' 
#' samp   <- read.flowSet(path=system.file("extdata", "compdata", "data", package="flowCore"))
#' cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
#' comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
#' comp.mat
#' 
#' ## create a compensation object 
#' comp <- compensation(comp.mat,compensationId="comp1")
#' ## create a compensated parameter object 
#' cPar1<-compensatedParameter(c("FL1-H","FL3-H"),"comp",searchEnv=.GlobalEnv)
#' compOut<-eval(cPar1)(exprs(samp[[1]]))
#' 
#' 
#' @export
setClass("compensatedParameter",
          contains=c("transform"),
          representation=representation(parameters="character",spillRefId="character",
                                        searchEnv="environment"
                                       )
        )

## Constructor
#' @export
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

#' Class "fcReference" and its subclasses
#' 
#' Classes and methods to provide reference-based semantics for flow cytometry
#' workflows.
#' 
#' 
#' These classes provide references to objects within an R environment and
#' allow for method dispatch based on the nature of the referenced object. The
#' parent \code{fcReference} class is used for references to all R objects,
#' unless there exists a more specific subclass. \code{fcTreeReference},
#' \code{fcViewReference}, and \code{fcActionReference} are used to reference
#' to \code{\link[graph:graphNEL-class]{graphNEL}}, \code{\link{view}}, and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' \code{fcDataReference} should be used for \code{\link{flowFrame}} or
#' \code{\link{flowSet}} objects, whereas \code{fcFilterResultReference},
#' \code{fcFilterReference}, \code{fcTransformReference},
#' \code{fcCompensateReference}, and \code{fcNormalizationReference} link to
#' \code{\link{filterResult}}, \code{\link{filter}},
#' \code{\linkS4class{transform}} and \code{\linkS4class{compensation}}
#' objects. \code{fsStructureReference} only exists to jointly dispatch on
#' certain subgroups of references.
#' 
#' @name fcReference-class
#' @aliases fcReference-class fcReference fcStructureReference-class
#' fcTreeReference-class fcTreeReference fcAliasReference-class
#' fcAliasReference fcDataReference-class fcDataReference
#' fcActionReference-class fcActionReference fcViewReference-class
#' fcViewReference fcFilterResultReference-class fcFilterResultReference
#' fcFilterReference-class fcFilterReference fcCompensateReference-class
#' fcCompensateReference fcTransformReference-class fcTransformReference
#' fcNormalizationReference-class fcNormalizationReference
#' fcNullReference-class fcSubsettingReference fcSubsettingReference-class
#' assign fcNullReference get,fcReference,missing,missing-method
#' get,fcNullReference,missing,missing-method identifier,fcReference-method
#' isNull isNull,fcReference-method Rm Rm,fcReference,missing,character-method
#' Rm,fcReference,workFlow,character-method
#' Rm,fcNullReference,missing,character-method show,fcReference-method
#' show,fcNullReference-method fcJournalReference-class fcJournalReference
#' @docType class
#' @usage 
#' fcReference(ID=paste("genericRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcTreeReference(ID=paste("treeRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcAliasReference(ID=paste("aliasRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcDataReference(ID=paste("dataRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcActionReference(ID=paste("actionRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcViewReference(ID=paste("viewRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcFilterResultReference(ID=paste("fresRef", guid(), sep="_"),
#'                         env=new.env(parent=emptyenv()))
#' 
#' fcFilterReference(ID=paste("filterRef", guid(), sep="_"), env=new.env(parent=emptyenv()))
#' 
#' fcCompensateReference(ID=paste("compRef",
#'                                guid(), sep="_"),
#'                       env=new.env(parent=emptyenv()))
#' 
#' fcNormalizationReference(ID=paste("normRef",
#'                                   guid(), sep="_"),
#'                          env=new.env(parent=emptyenv()))
#' 
#' fcSubsettingReference(ID=paste("subRef",
#'                                guid(), sep="_"),
#'                       env=new.env(parent=emptyenv()))
#' 
#' fcTransformReference(ID=paste("transRef",
#'                               guid(), sep="_"),
#'                      env=new.env(parent=emptyenv()))
#' 
#' fcNullReference(...)
#' 
#' 
#' assign(x, value, pos = -1, envir = as.environment(pos), inherits = FALSE, 
#'        immediate = TRUE)
#' 
#' \S4method{get}{fcReference,missing,missing}(x, pos = -1, envir = as.environment(pos), mode = "any", 
#'                                             inherits = TRUE)
#' 
#' isNull(f)
#' 
#' Rm(symbol, envir, subSymbol, ...)
#' 
#' @param x,f,symbol An object of class or inheriting from class
#' \code{fcReference}.
#' @param ID The reference identifier.
#' @param value An arbitrary R object which is supposed to be assigned to the
#' environment in the \code{\link{workFlow}} object and to which a reference is
#' returned.
#' @param env An environment, usually within a \code{\link{workFlow}} object.
#' @param pos,envir Objects of class \code{\link{workFlow}}.
#' @param inherits,immediate,mode,subSymbol,\dots Further arguments from the
#' generics that are not used in this context.
#' @return
#' 
#' An object of class \code{fcReference} or one of its subclasses for the
#' \code{assign} constructor.
#' 
#' The object referenced to for the \code{get} method.
#' 
#' A character string of the object symbol for the \code{identifier} method.
#' 
#' A logical scalar for the \code{isNull} method.
#' @section Extends:
#' 
#' \code{fcStructureReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcTreeReference}:
#' 
#' Class \code{"\linkS4class{fcStructureReference}"}, directly.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcStructureReference",
#' distance 2.
#' 
#' \code{fcAliasReference}:
#' 
#' Class \code{"\linkS4class{fcStructureReference}"}, directly.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcStructureReference",
#' distance 2.
#' 
#' \code{fcDataReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcActionReference}:
#' 
#' Class \code{"\linkS4class{fcStructureReference}"}, directly.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcStructureReference",
#' distance 2.
#' 
#' \code{fcViewReference}:
#' 
#' Class \code{"\linkS4class{fcStructureReference}"}, directly.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcStructureReference",
#' distance 2.
#' 
#' \code{fcFilterResultReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcFilterReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcCompensateReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcTransformReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcNormalizationReference}:
#' 
#' Class \code{"\linkS4class{fcReference}"}, directly.
#' 
#' \code{fcNullReference}:
#' 
#' Class \code{"\linkS4class{fcDataReference}"}, directly.  Class
#' \code{"\linkS4class{fcActionReference}"}, directly.  Class
#' \code{"\linkS4class{fcViewReference}"}, directly.  Class
#' \code{"\linkS4class{fcFilterResultReference}"}, directly.  Class
#' \code{"\linkS4class{fcFilterReference}"}, directly.  Class
#' \code{"\linkS4class{fcCompensateReference}"}, directly.  Class
#' \code{"\linkS4class{fcTransformReference}"}, directly.  Class
#' \code{"\linkS4class{fcNormalizationReference}"}, directly.  Class
#' \code{"\linkS4class{fcTreeReference}"}, directly.  Class
#' \code{"\linkS4class{fcAliasReference}"}, directly.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcDataReference", distance2.
#' Class \code{"\linkS4class{fcStructureReference}"}, by class
#' "fcActionReference", distance 2.  Class \code{"\linkS4class{fcReference}"},
#' by class "fcActionReference", distance 3.  Class
#' \code{"\linkS4class{fcStructureReference}"}, by class "fcViewReference",
#' distance 2.  Class \code{"\linkS4class{fcReference}"}, by class
#' "fcViewReference", distance3.  Class \code{"\linkS4class{fcReference}"}, by
#' class "fcFilterResultReference", distance 2.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcFilterReference", distance
#' 2.  Class \code{"\linkS4class{fcReference}"}, by class
#' "fcCompensateReference", distance 2.  Class
#' \code{"\linkS4class{fcReference}"}, by class "fcTransformReference",
#' distance 2.  Class \code{"\linkS4class{fcStructureReference}"}, by class
#' "fcTreeReference", distance 2.  Class \code{"\linkS4class{fcReference}"}, by
#' class "fcTreeReference", distance 3.
#' 
#' @section Objects from the Class:
#' Objects should be created via the \code{assign} constructor. Whenever
#' an object is assigned to a \code{\link{workFlow}} using the
#' \code{assign} method, an appropriate instance of class
#' \code{fcReference} or one of its subclasses is returned. In addition,
#' there are the usual constructor functions of same names as the classes
#' that can be used for object instantiation without assignment. Note
#' that this might lead to unresolvable references unless the object
#' referenced to is available in the environment.
#' 
#' @slot ID Object of class \code{"character"} The name of the
#' object in \code{env} referenced to.
#' @slot env Object of class \code{"environment"} An
#' environment that contains the referenced objects. Usually, this
#' will be the environment that's part of a \code{\link{workFlow}}
#' object.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{get}{\code{signature(x = "fcReference", pos = "missing", envir
#'                              = "missing", mode = "missing", inherits = "missing")}: Resolve
#'     the reference, i.e., get the object from the environment. }
#'   
#'   \item{get}{\code{signature(x = "fcNullReference", pos = "missing",
#'                              envir = "missing", mode = "missing", inherits = "missing")}: Resolve
#'     the reference. This always returns \code{NULL}. }
#'   
#'   \item{identifier}{\code{signature(object = "fcReference")}: Return
#'     a character string of the object name. }
#'   
#'   \item{isNull}{\code{signature(f = "fcReference")}: Check whether a
#'     \code{fcReference} is a \code{NULL} reference. Note that this is
#'     different from a unresolvable reference. }
#'   
#'   \item{Rm}{\code{signature(symbol = "fcReference", envir = "missing",
#'                             subSymbol = "character")}: Remove the object referenced to by a
#'     \code{fcReference} from its environment. The argument
#'     \code{subSymbol} will be automatically set by the generic and
#'     should never be provided by the user. }
#'   
#'   \item{Rm}{\code{signature(symbol = "fcReference", envir =
#'                               "workFlow", subSymbol = "character")}: Remove the object referenced to by a
#'     \code{fcReference} from a \code{\link{workFlow}}. The argument
#'     \code{subSymbol} will be automatically set by the generic and
#'     should never be provided by the user. }
#'   
#'   \item{Rm}{\code{signature(symbol = "fcNullReference", envir =
#'                               "missing", subSymbol = "character")}: Essentially, this doesn't do
#'       anything since there is no object referenced to. }
#'     
#'     \item{show}{\code{signature(object = "fcReference")}: Print details
#'       about the object. }
#'     
#'     \item{show}{\code{signature(object = "fcNullReference")}: Print details
#'       about the object. }
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\link{workFlow}}
#' @keywords classes
#' @examples
#' 
#' showClass("fcReference")
#' 
#' @export

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
#' @export
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
#' @export
setClass("fcStructureReference",
         contains=list("VIRTUAL",
         "fcReference"))



## ===========================================================================
## fcTreeReference
## ---------------------------------------------------------------------------
## A reference to a graphNEL object representing the workflow tree
## ---------------------------------------------------------------------------
#' @export
setClass("fcTreeReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("treeRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcJournalReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("journalRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcAliasReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("aliasRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcDataReference",
         contains="fcReference",
         prototype=prototype(ID=paste("dataRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcActionReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("actionRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcViewReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("viewRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcFilterResultReference",
         contains="fcReference",
         prototype=prototype(ID=paste("fresRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcFilterReference",
         contains="fcReference",
         prototype=prototype(ID=paste("filterRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcCompensateReference",
         contains="fcReference",
         prototype=prototype(ID=paste("compRef", guid(), sep="_"))
         )

## constructor
#' @export
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
## objects until transformation is a more useful class. We need this to
## store a trasnforamtion within a transformActionItem without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
#' @export
setClass("fcTransformReference",
         contains="fcReference",
         prototype=prototype(ID=paste("transRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcNormalizationReference",
         contains="fcReference",
         prototype=prototype(ID=paste("normRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
setClass("fcSubsettingReference",
         contains="fcReference",
         prototype=prototype(ID=paste("subRef", guid(), sep="_"))
         )

## constructor
#' @export
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
#' @export
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

#' @export
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
#' Class "normalization"
#' 
#' 
#' Class and methods to normalize a a \code{flowSet} using a potentially
#' complex normalization function.
#' 
#' Data normalization of a \code{flowSet} is a rather fuzzy concept and the
#' class mainly exists for method dispatch in the workflow tools. The idea is
#' to have a rather general function that takes a \code{flowSet} and a list of
#' parameter names as input and applies any kind of normalization to the
#' respective data columns. The output of the function has to be a
#' \code{flowSet} again. Although we don't formally check for it, the
#' dimensions of the input and of the output set should remain the same.
#' Additional arguments may be passed to the normalization function via the
#' \code{arguments} list. Internally we evaluate the function using
#' \code{\link{do.call}} and one should check its documentation for details.
#' 
#' Currently, the most prominent example for a normalization function is
#' warping, as provided by the \code{flowStats} package.
#' 
#' @name normalization-class
#' @aliases normalization-class normalization normalize
#' add,workFlow,normalization-method
#' identifier<-,normalization,character-method identifier,normalization-method
#' normalize,flowSet,normalization-method parameters,normalization-method
#' @docType class
#' @usage
#' normalization(parameters, normalizationId="defaultNormalization",
#'               normFunction, arguments=list())
#'
#' normalize(data, x,...)
#' @param parameters Character vector of parameter names.
#' @param normalizationId The identifier for the normalization object.
#' @param x An object of class \code{\linkS4class{flowSet}}.
#' @param normFunction The normalization function
#' @param arguments The list of additional arguments to \code{normFunction}
#' @param data The \code{flowSet} to normalize.
#' @param \dots other arguments: see
#' \code{\link[flowStats:normalize-methods]{normalize-methods}}for details.
#' 
#' @return
#' 
#' A \code{normalization} object for the constructor.
#' 
#' A \code{\linkS4class{flowSet}} for the \code{normalize} methods.
#' @section Objects from the Class:
#' 
#' Objects should be created using the constructor \code{normalization()}. See
#' the \code{Usage} and \code{Arguments} sections for details.
#' 
#' @slot parameters Object of class \code{"character"}. The
#' flow parameters that are supposed to be normalized by the
#' normalization function.
#' @slot normalizationId Object of class \code{"character"}. An
#' identifier for the object.
#' @slot normFunction Object of class \code{"function"} The
#' normalization function. It has to take two mandatory arguments:
#' \code{x}, the \code{flowSet}, and \code{parameters}, a character
#' of parameter names that are to be normalized by the
#' function. Additional arguments have to be passed in via
#' \code{arguments}.
#' @slot arguments Object of class \code{"list"} A names list
#' of additional arguments. Can be \code{NULL}.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "normalization")}: The constructor for the workFlow. }
#'   
#'   \item{identifier<-}{\code{signature(object = "normalization", value
#'                                       = "character")}: Set method for the identifier slot. }
#'   
#'   \item{identifier}{\code{signature(object = "normalization")}: Get
#'     method for the identifier slot. }
#'   
#'   \item{normalize}{\code{signature(data = "flowSet", x =
#'                                      "normalization")}: Apply a normalization to a \code{\linkS4class{flowSet}}. }
#'   
#'   \item{parameters}{\code{signature(object = "normalization")}: The
#'     more generic constructor. }
#' }
#' @author F. Hahne
#' @keywords methods classes
#'
#' @export
setClass("normalization",
         representation(parameters="character",
                        normalizationId="character",
                        normFunction="function",
                        arguments="list"),
         prototype=prototype(normalizationId="defaultNormalization",
                             normFunction=function(x) x)
         )

## constructor
#' @export
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
#' Class "characterOrNumeric"
#' 
#' A simple union class of \code{character} and \code{numeric}.
#' Objects will be created internally whenever necessary and the user should
#' not need to explicitly interact with this class.
#' 
#' @name characterOrNumeric-class
#' @aliases characterOrNumeric-class characterOrNumeric
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' @keywords classes
#' @examples
#' 
#' showClass("characterOrNumeric")
#' 
setClassUnion("characterOrNumeric", c("character","numeric"))

#' Class "subsetting"
#' 
#' 
#' Class and methods to subset a a \code{flowSet}. This is only needed for
#' method dispatch in the workFlow framework.
#' 
#' The class mainly exists for method dispatch in the workflow tools.
#' 
#' @name subsetting-class
#' @aliases subsetting-class subsetting add,workFlow,subsetting-method
#' add,workFlow,numeric-method add,workFlow,character-method
#' add,workFlow,logical-method identifier<-,subsetting,character-method
#' identifier,subsetting-method show,subsetting-method
#' @docType class
#' @usage subsetting(indices, subsettingId="defaultSubsetting")
#' @param indices Character or numeric vector of sample names.
#' @param subsettingId The identifier for the subsetting object.
#' @return
#' 
#' A \code{subsetting} object.
#' @section Objects from the Class:
#' 
#' Objects should be created using the constructor \code{subsetting()}. See the
#' \code{Usage} and \code{Arguments} sections for details.
#' 
#' @slot subsettingId Object of class \code{"character"}. An
#' identifier for the object.
#' @slot indices Object of class
#' \code{"numericOrCharacter"}. Indices into a \code{flowSet}.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "subsetting")}: The constructor for the workFlow. }
#'   
#'   \item{identifier<-}{\code{signature(object = "subsetting", value
#'                                       = "character")}: Set method for the identifier slot. }
#'   
#'   \item{identifier}{\code{signature(object = "subsetting")}: Get
#'     method for the identifier slot. }
#'   
#'   \item{show}{\code{signature(object = "subsetting")}: Show details
#'     about the object. }
#'   
#' }
#' 
#' 
#' @author F. Hahne
#' @keywords methods classes
#'
#' @export
setClass("subsetting",
         representation(subsettingId="character",
                        indices="characterOrNumeric"),
         prototype=prototype(subsettingId="defaultSubsetting")
         )

## constructor
#' @export
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
#' Class "workFlow"(deprecated)
#' 
#' Class and methods to organize standard flow cytometry data analysis
#' operation in a concise workflow.
#' 
#' 
#' \code{workFlow} objects organize standard flow data analysis operations like
#' gating, compensation and transformation in one single object. The user can
#' interact with a \code{workFlow} object (e.g. adding operations, removing
#' them, summarizing the results) without having to keep track of intermediate
#' objects and names.
#' 
#' The integral part of a \code{workFlow} is an evaluation environment which
#' holds all objects that are created during the analysis. The structure of the
#' whole workflow is a tree, where nodes represent \code{link{view}}s (or
#' results of an operation) and edges represent
#' \code{\link[flowCore:actionItem-class]{actionItems}} (or the operations
#' themselves).
#' 
#' @name workFlow-class
#' @aliases workFlow-class workFlow add add,workFlow,concreteFilter-method
#' add,workFlow,filterList-method add,workFlow,transformList-method
#' add,workFlow,compensation-method
#' assign,ANY,ANY,missing,workFlow,missing,missing-method
#' assign,missing,ANY,workFlow,missing,missing,missing-method
#' assign,missing,ANY,missing,workFlow,missing,missing-method
#' assign,character,ANY,workFlow,missing,missing,missing-method
#' assign,fcReference,ANY,workFlow,missing,missing,missing-method
#' [,workFlow,ANY-method [[,workFlow,ANY-method $,workFlow-method
#' get,character,workFlow,missing-method get,character,missing,workFlow-method
#' ls ls,workFlow,missing,missing,missing,missing-method
#' ls,workFlow,missing,missing,missing,character-method
#' mget,character,workFlow-method names,workFlow-method
#' plot,workFlow,missing-method Rm,character,workFlow,character-method
#' show,workFlow-method views views,workFlow-method actions
#' actions,workFlow-method alias,workFlow-method alias,environment-method
#' nodes,workFlow-method summary,workFlow-method journal
#' journal,workFlow-method journal,environment-method tree tree,workFlow-method
#' tree,environment-method undo
#' @docType class
#' @usage 
#' workFlow(data, name = "default", env = new.env(parent = emptyenv()))
#' 
#' undo(wf, n=1)
#' @param data An object of class \code{\link{flowFrame}} or
#' \code{\link{flowSet}} for which a basic \code{\link{view}} is created.
#' @param name A more human-readable name of the view.
#' @param env Object of class \code{environment}. The evaluation environment
#' used for the \code{\link{workFlow}}.
#' @param wf Object of class \code{workFlow}.
#' @param n The number of operations to undo.
#' @return
#' 
#' A \code{workFlow} object for the constructor
#' 
#' Both \code{applyParentFilter} and \code{undo} are called for their
#' side-effects.
#' @section Objects from the Class:
#' 
#' Objects should be created using the constructor \code{workFlow}, which takes
#' a \code{\link{flowFrame}} or \code{\link{flowSet}} as only mandatory input
#' and creates a basic view for that.
#' 
#' @slot name Object of class \code{"character"}. The name of
#' the workFlow object.
#' @slot tree Object of class \code{"fcTreeReference"}. A
#' reference to the \code{\link[graph:graphNEL-class]{graphNEL}}
#' objects representing the view structure of the workflow.
#' @slot alias Object of class \code{"fcAliasReference"}. A
#' reference to the alias table.
#' @slot journal Object of class \code{"fcJournalReference"}. A
#' reference to the journal.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment for the workflow in which all objects will be
#' stored.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "concreteFilter")}: Create a new \code{\link{gateActionItem}} and
#'     \code{\link{gateView}} from a \code{\link{filter}} and assign those to
#'     the workflow. }
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "filterList")}: Create a new \code{\link{gateActionItem}} and
#'     \code{\link{gateView}} from a \code{\link{filterList}}
#'     and assign those to the workflow. }
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "transformList")}: Create a new \code{\link{transformActionItem}} and
#'     \code{\link{transformView}} from a \code{\link{transform}}
#'     and assign those to the workflow. }
#'   
#'   \item{add}{\code{signature(wf = "workFlow", action =
#'                                "compensation")}: Create a new \code{\link{compensateActionItem}} and
#'     \code{\link{compensateView}} from a \code{\link{compensation}}
#'     and assign those to the workflow.}
#'   
#'   \item{assign}{\code{signature(x = "ANY", value = "ANY", pos =
#'                                   "missing", envir = "workFlow", inherits = "missing", immediate =
#'                                   "missing")}: Assign an object to the environment in the
#'     \code{workFlow} object and return a \code{\link{fcReference}} to
#'     it. The symbol for the object is created as a unique identifier. }
#'   
#'   \item{assign}{\code{signature(x = "missing", value = "ANY", pos =
#'                                   "workFlow", envir = "missing", inherits = "missing", immediate =
#'                                   "missing")}: see above }
#'   
#'   \item{assign}{\code{signature(x = "missing", value = "ANY", pos =
#'                                   "missing", envir = "workFlow", inherits = "missing", immediate =
#'                                   "missing")}: same as above, but provide custom symbol for the
#'     assignment.}
#'   
#'   \item{assign}{\code{signature(x = "character", value = "ANY", pos =
#'                                   "workFlow", envir = "missing", inherits = "missing", immediate =
#'                                   "missing")}: see above}
#'   
#'   \item{assign}{\code{signature(x = "fcReference", value = "ANY", pos
#'                                 = "workFlow", envir = "missing", inherits = "missing", immediate
#'                                 = "missing")}: same as above, but assign object using an
#'     existing \code{\link{fcReference}}. Note that assigning
#'     \code{NULL} essentially removes the original object. }
#'   
#'   \item{[}{\code{signature(x = "workFlow", i = "ANY")}: Cast a useful
#'     error message. }
#'   
#'   \item{[[}{\code{signature(x = "workFlow", i = "ANY")}: Treat the
#'     \code{workFlow} object as a regular environment. Essentially, this
#'     is equivalent to \code{get(x, i)}.  }
#'   
#'   \item{$}{\code{signature(x = "workFlow", name = "character")}:
#'       Allow for list-like access. Note that completion is only
#'     available for \code{\link[flowCore:view-class]{views}} since all other
#'     objects in the environment are considered to be internal. }
#'   
#'   \item{get}{\code{signature(x = "character", pos = "workFlow", envir
#'                              = "missing", mode = "missing", inherits = "missing")}: Get an
#'     object identified by symbol \code{x} from the environment in the
#'     \code{workFlow}. }
#'   
#'   \item{get}{\code{signature(x = "character", pos = "missing", envir =
#'                                "workFlow", mode = "missing", inherits = "missing")}: see above }
#'   
#'   \item{ls}{\code{signature(name = "workFlow", pos = "missing", envir
#'                             = "missing", all.names = "missing", pattern = "missing")}: List
#'     the content of the environment in the \code{workFlow}. }
#'   
#'   \item{ls}{\code{signature(name = "workFlow", pos = "missing", envir
#'                             = "missing", all.names = "missing", pattern = "character")}: see
#'     above}
#'   
#'   \item{mget}{\code{signature(x = "character", envir = "workFlow",
#'                               mode = "missing", ifnotfound = "missing", inherits =
#'                                 "missing")}: Get multiple objects identified by the symbols in
#'     \code{x} from the environment in the \code{workFlow}. }
#'   
#'   \item{names}{\code{signature(x = "workFlow")}: List the identifiers
#'     for all \code{\link[flowCore:view-class]{views}} and
#'     \code{\link[flowCore:actionItem-class]{actionItems}} in the
#'     \code{workFlow}. }
#'   
#'   \item{plot}{\code{signature(x = "workFlow", y = "missing")}: Plot
#'     the structure of the \code{workFlow} tree. }
#'   
#'   \item{Rm}{\code{signature(symbol = "character", envir = "workFlow",
#'                             subSymbol = "character")}: Remove the object identified by the
#'     symbol \code{symbol} from the \code{workFlow}. }
#'   
#'   \item{undo}{\code{signature(wf = "workFlow", n = "numeric")}: Undo
#'     the last \code{n} operations on the  \code{workFlow}. }
#'   
#'   \item{show}{\code{signature(object = "workFlow")}: Print details
#'     about the object. }
#'   
#'   \item{summary}{\code{signature(object = "workFlow")}: Summarize a
#'     view in the \code{workFlow}. }
#'   
#'   \item{nodes}{\code{signature(object = "workFlow")}: Return a named vector
#'     of node ids where the names are the human readable names stored in
#'     the alias table. }
#'   
#'   \item{actions}{\code{signature(x = "workFlow")}: List the names of
#'     the \code{\link[flowCore:actionItem-class]{actionItems}} in the
#'     \code{workFlow}. }
#'   
#'   \item{views}{\code{signature(x = "workFlow")}: List the names of
#'     only the \code{\link[flowCore:view-class]{views}} in the
#'     \code{workFlow}. }
#'   
#'   \item{alias}{\code{signature(object = "workFlow")}: Return the alias
#'     table for the \code{workFlow}.}
#'   
#'   \item{alias}{\code{signature(object = "environment")}: Return the
#'     alias table from a generic environment. The method tries to find
#'     'fcAliasRef' among the object symbols in the environment. }
#'   
#'   \item{journal}{\code{signature(object = "workFlow")}: Return the
#'     journal for the \code{workFlow}.}
#'   
#'   \item{journal}{\code{signature(object = "environment")}: Return the
#'     journal from a generic environment. The method tries to find
#'     'fcJournalRef' among the object symbols in the environment. }
#'   
#'   \item{tree}{\code{signature(object = "workFlow")}: Return the
#'     tree of the \code{workFlow}.}
#'   
#'   \item{journal}{\code{signature(object = "environment")}: Return the
#'     tree from a generic environment. The method tries to find
#'     'fcTreeRef' among the object symbols in the environment. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{"\linkS4class{view}"}, \code{"\linkS4class{actionItem}"}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export 
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
#' @export
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
#' @export
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
#' @export
setMethod("assign",
          signature=signature(x="missing",
          value="ANY",
          pos="missing",
          envir="workFlow",
          inherits="missing",
          immediate="missing"),
          definition=function(value, envir) assign(value=value, pos=envir))

## Assign to a particular symbol (potentially overwriting existing ones)
#' @export
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
#' @export
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
#' @export
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
#' Class "actionItem"
#' 
#' Class and method to capture standard operations in a flow cytometry
#' workflow.
#' 
#' 
#' \code{actionItems} provide a means to bind standard operations on flow
#' cytometry data in a workflow. Usually, the user doesn't have to create these
#' objects, instead they will be automatically created when applying one of the
#' standard operations (gating, transformation, compensation) to a
#' \code{\link{workFlow}} object. Each \code{actionItem} creates one or several
#' new \code{\link[flowCore:view-class]{views}}, which again can be the basis
#' to apply further operations. One can conceptualize \code{actionItems} being
#' the edges in the workflow tree connecting
#' \code{\link[flowCore:view-class]{views}}, which are the nodes of the tree.
#' There are more specific subclasses for the three possible types of
#' operation: \code{\link{gateActionItem}} for gating operations,
#' \code{\link{transformActionItem}} for transformations, and
#' \code{\link{compensateActionItem}} for compensation operations. See their
#' documentation for details.
#' 
#' @name actionItem-class
#' @aliases actionItem-class actionItem identifier,actionItem-method
#' names,actionItem-method alias,actionItem-method parent,actionItem-method
#' Rm,actionItem,workFlow,character-method
#' @docType class
#' @section Objects from the Class: A virtual Class: No objects may be created
#' from it.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot parentView Object of class \code{"fcViewReference"}. A
#' reference to the parent \code{\link{view}} the \code{actionItem}
#' is applied on.
#' @slot alias Object of class \code{"fcAliasReference"}. A
#' reference to the alias table.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{identifier}{\code{signature(object = "actionItem")}: Accessor
#'     for the \code{ID} slot. }
#'   
#'   \item{names}{\code{signature(x = "actionItem")}: Accessor
#'     for the \code{name} slot.}
#'   
#'   \item{parent}{\code{signature(object = "actionItem")}: Accessor for
#'     the \code{parentView} slot. Note that the reference is resolved,
#'     i.e., the \code{\link{view}} object is returned. }
#'   
#'   \item{alias}{\code{signature(object = "actionItem")}: Get the alias table
#'     from a \code{actionItem}. }
#'   
#'   \item{Rm}{\code{signature(symbol = "actionItem", envir = "workFlow",
#'                             subSymbol = "character")}:  Remove a \code{actionItem} from a
#'     \code{\link{workFlow}}. This method is recursive and will also
#'     remove all dependent \code{\link[flowCore:view-class]{views}} and
#'     \code{actionItems}. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{gateActionItem}},
#' \code{\linkS4class{transformActionItem}},
#' \code{\linkS4class{compensateActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' Class "gateActionItem"
#' 
#' 
#' Class and method to capture gating operations in a flow cytometry workflow.
#' 
#' 
#' \code{gateActionItems} provide a means to bind gating operations in a
#' workflow. Each \code{gateActionItem} represents a single
#' \code{\link{filter}}.
#' 
#' @name gateActionItem-class
#' @aliases gateActionItem-class gateActionItem gate gate,gateActionItem-method
#' print,gateActionItem-method Rm,gateActionItem,workFlow,character-method
#' show,gateActionItem-method summary,gateActionItem-method
#' @docType class
#' @usage
#' gateActionItem(ID = paste("gateActionRef", guid(), sep = "_"), name =
#'                paste("action", identifier(get(gate)), sep = "_"), parentView, gate,
#'                filterResult, workflow)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param parentView,gate,filterResult References to the parent
#' \code{\link{view}}, \code{\link{filter}}, and \code{\link{filterResult}}
#' objects, respectively.
#' @return
#' 
#' A reference to the \code{gateActionItem} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{gateActionItem} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{gateActionItem} from a \code{\link{filter}} object and directly
#' assigns it to a \code{\link{workFlow}}. Alternatively, one can use the
#' \code{gateActionItem} constructor function for more programmatic access.
#' 
#' @slot gate Object of class \code{"fcFilterReference"}. A
#' reference to the \code{\link{filter}} that is used for the gating
#' operation.
#' @slot filterResult Object of class
#' \code{"fcFilterResultReference"}. A reference to the
#' \code{\link{filterResult}} produced by the gating operation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name
#' @slot parentView Object of class \code{"fcViewReference"}. A
#' reference to the parent \code{\link{view}} the \code{gateActionItem}
#' is applied on.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{gate}{\code{signature(object = "gateActionItem")}: Accessor to
#'     the \code{gate} slot. Note that this resolved the reference, i.e.,
#'     the \code{\link{filter}} object is returned. }
#'   
#'   \item{print}{\code{signature(x = "gateActionItem")}: Print details
#'     about the object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "gateActionItem", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a \code{gateActionItem}
#'     from a \code{\link{workFlow}}. This method is recursive and will
#'     also remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItem}}s.}
#'   
#'   \item{show}{\code{signature(object = "gateActionItem")}: Print details
#'     about the object. }
#'   
#'   \item{summary}{\code{signature(object = "gateActionItem")}:
#'       Summarize the gating operation and return the appropriate
#'     \code{\link{filterSummary}} object. }
#'   
#' }
#' 
#' @section Extends:
#' Class \code{"\linkS4class{actionItem}"}, directly.
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{actionItem}},
#' \code{\linkS4class{transformActionItem}},
#' \code{\linkS4class{compensateActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' @export
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
#' Class "transformActionItem"
#' 
#' 
#' Class and method to capture transformation operations in a flow cytometry
#' workflow.
#' 
#' 
#' \code{transformActionItems} provide a means to bind transformation
#' operations in a workflow. Each \code{transformActionItem} represents a
#' single \code{\link{transform}}.
#' 
#' @name transformActionItem-class
#' @aliases transformActionItem-class transformActionItem
#' print,transformActionItem-method
#' Rm,transformActionItem,workFlow,character-method
#' show,transformActionItem-method
#' @docType class
#' @usage
#' transformActionItem(ID = paste("transActionRef", guid(), sep = "_"),
#'                     name=paste("action", identifier(get(transform)), sep = "_"), 
#'                     parentView, transform, workflow)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param parentView,transform References to the parent \code{\link{view}} and
#' \code{\link{transform}} objects, respectively.
#' @return
#' 
#' A reference to the \code{transformActionItem} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{transformActionItem} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{transformActionItem} from a \code{\link{transform}} object and
#' directly assigns it to a \code{\link{workFlow}}. Alternatively, one can use
#' the \code{transformActionItem} constructor function for more programmatic
#' access.
#' 
#' @slot transform Object of class
#' \code{"fcTransformReference"}. A reference to the
#' \code{\link{transform}} object that is used for the
#' transformation operation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot parentView Object of class
#' \code{"fcViewReference"}. A reference to the parent
#' \code{\link{view}} the \code{transformActionItem} is applied on.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{actionItem}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   \item{print}{\code{signature(x = "transformActionItem")}: Print
#'     details about the object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "transformActionItem", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{transformActionItem} from a \code{\link{workFlow}}. This
#'     method is recursive and will also remove all dependent
#'     \code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}.}
#'   
#'   \item{show}{\code{signature(object = "transformActionItem")}: Print
#'     details about the object. }
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{actionItem}},
#' \code{\linkS4class{gateActionItem}},
#' \code{\linkS4class{compensateActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export 
setClass("transformActionItem",
         contains="actionItem",
         representation=representation(transform="fcTransformReference"))

## The constructor creates the transformActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
#' @export
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
#' Class "compensateActionItem"
#' 
#' 
#' Class and method to capture compensation operations in a flow cytometry
#' workflow.
#' 
#' 
#' \code{compensateActionItems} provide a means to bind compensation operations
#' in a workflow. Each \code{compensateActionItem} represents a single
#' \code{\link{compensation}}.
#' 
#' @name compensateActionItem-class
#' @aliases compensateActionItem-class compensateActionItem
#' print,compensateActionItem-method
#' Rm,compensateActionItem,workFlow,character-method
#' show,compensateActionItem-method
#' @docType class
#' @usage
#' compensateActionItem(ID = paste("compActionRef", guid(), sep = "_"),
#'                      name = paste("action", identifier(get(compensate)), sep = "_"),
#'                      parentView, compensate, workflow)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param parentView,compensate References to the parent \code{\link{view}} and
#' \code{\link{compensation}} objects, respectively.
#' @return
#' 
#' A reference to the \code{compensateActionItem} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{compensateActionItem} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{compensateActionItem} from a \code{\link{compensation}} object and
#' directly assigns it to a \code{\link{workFlow}}. Alternatively, one can use
#' the \code{compensateActionItem} constructor function for more programmatic
#' access.
#' 
#' @slot compensate Object of class
#' \code{"fcCompensateReference"}. A reference to the
#' \code{\link{compensation}} object that is used for the
#' compensation operation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot parentView Object of class
#' \code{"fcViewReference"}. A reference to the parent
#' \code{\link{view}} the \code{compensateActionItem} is applied on.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{actionItem}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{print}{\code{signature(x = "compensateActionItem")}: Print
#'     details about the object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "compensateActionItem", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{compensateActionItem} from a \code{\link{workFlow}}. This
#'     method is recursive and will also remove all dependent
#'     \code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
#'   
#'   \item{show}{\code{signature(object = "compensateActionItem")}: Print
#'     details about the object. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{actionItem}},
#' \code{\linkS4class{gateActionItem}},
#' \code{\linkS4class{transformActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
setClass("compensateActionItem",
         contains="actionItem",
         representation=representation(compensate="fcCompensateReference"))

## The constructor creates the compensateActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
#' @export
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
#' Class "normalizeActionItem"
#' 
#' 
#' Class and method to capture normalization operations in a flow cytometry
#' workflow.
#' 
#' 
#' \code{normalizeActionItems} provide a means to bind normalization operations
#' like warping in a workflow. Each \code{normalizeActionItem} represents a
#' single \code{\link{normalization}}.
#' 
#' @name normalizeActionItem-class
#' @aliases normalizeActionItem-class normalizeActionItem
#' print,normalizeActionItem-method
#' Rm,normalizeActionItem,workFlow,character-method
#' show,normalizeActionItem-method
#' @docType class
#' @usage
#' normalizeActionItem(ID = paste("normActionRef", guid(), sep = "_"),
#'                     name = paste("action", identifier(get(normalization)), sep = "_"),
#'                     parentView, normalization, workflow)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param parentView,normalization References to the parent \code{\link{view}}
#' and \code{\link{normalization}} objects, respectively.
#' @return
#' 
#' A reference to the \code{normalizeActionItem} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{normalizeActionItem} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{normalizeActionItem} from a \code{\link{normalization}} object and
#' directly assigns it to a \code{\link{workFlow}}. Alternatively, one can use
#' the \code{normalizeActionItem} constructor function for more programmatic
#' access.
#' 
#' @slot normalization Object of class
#' \code{"fcNormalizationReference"}. A reference to the
#' \code{\link{normalization}} object that is used for the
#' compensation operation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot parentView Object of class
#' \code{"fcViewReference"}. A reference to the parent
#' \code{\link{view}} the \code{normalizeActionItem} is applied on.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{actionItem}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{print}{\code{signature(x = "normalizeActionItem")}: Print
#'     details about the object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "normalizeActionItem", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{normalizeActionItem} from a \code{\link{workFlow}}. This
#'     method is recursive and will also remove all dependent
#'     \code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
#'   
#'   \item{show}{\code{signature(object = "normalizeActionItem")}: Print
#'     details about the object. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{actionItem}},
#' \code{\linkS4class{gateActionItem}},
#' \code{\linkS4class{transformActionItem}},
#' \code{\linkS4class{compensateActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
setClass("normalizeActionItem",
         contains="actionItem",
         representation=representation(normalization="fcNormalizationReference"))

## The constructor creates the normalizeActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
#' @export
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
#' Class "subsettingActionItem"
#' 
#' 
#' Class and method to capture subsetting operations in a flow cytometry
#' workflow.
#' 
#' 
#' \code{subsettingActionItems} provide a means to bind subsetting operations
#' in a workflow. Each \code{subsettingActionItem} represents a single
#' \code{\link{subsetting}}.
#' 
#' @name subsettingActionItem-class
#' @aliases subsettingActionItem-class subsettingActionItem
#' print,subsettingActionItem-method
#' Rm,subsettingActionItem,workFlow,character-method
#' show,subsettingActionItem-method
#' @docType class
#' @usage 
#' subsettingActionItem(ID = paste("subActionRef", guid(), sep = "_"),
#'                      name = paste("action", identifier(get(subsetting)), sep = "_"),
#'                      parentView, subsetting, workflow)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param parentView,subsetting References to the parent \code{\link{view}} and
#' \code{\link{subsetting}} objects, respectively.
#' @return
#' 
#' A reference to the \code{subsettingActionItem} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{subsettingActionItem} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{subsettingActionItem} from a \code{\link{normalization}} object and
#' directly assigns it to a \code{\link{workFlow}}. Alternatively, one can use
#' the \code{subsettingActionItem} constructor function for more programmatic
#' access.
#' 
#' @slot subsetting Object of class
#' \code{"fcSubsettingReference"}. A reference to the
#' \code{\link{subsetting}} object that is used for the
#' operation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the \code{actionItem}.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot parentView Object of class
#' \code{"fcViewReference"}. A reference to the parent
#' \code{\link{view}} the \code{subsettingActionItem} is applied on.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{actionItem}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{print}{\code{signature(x = "subsettingActionItem")}: Print
#'     details about the object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "subsettingActionItem", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{subsettingActionItem} from a \code{\link{workFlow}}. This
#'     method is recursive and will also remove all dependent
#'     \code{views} and \code{\link[flowCore:actionItem-class]{actionItems}}. }
#'   
#'   \item{show}{\code{signature(object = "subsettingActionItem")}: Print
#'     details about the object. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{actionItem}},
#' \code{\linkS4class{gateActionItem}},
#' \code{\linkS4class{transformActionItem}},
#' \code{\linkS4class{compensateActionItem}}, \code{\linkS4class{view}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
setClass("subsettingActionItem",
         contains="actionItem",
         representation=representation(subsetting="fcSubsettingReference"))

## The constructor creates the subsettingActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
#' @export
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
#' Class "view"
#' 
#' Class and method to capture the results of standard operations (called
#' "views" here) in a flow cytometry workflow.
#' 
#' 
#' \code{Views} provide a means to bind the results of standard operations on
#' flow cytometry data in a workflow. Each view can be considered the outcome
#' of one operation. There are more specific subclasses for the three possible
#' types of operation: \code{\link{gateView}} for gating operations,
#' \code{\link{transformView}} for transformations, and
#' \code{\link{compensateView}} for compensation operations. See their
#' documentation for details.
#' 
#' @name view-class
#' @aliases view-class view action,view-method action alias alias,view-method
#' Data,view-method Data names,view-method parent parent,view-method
#' parent,NULL-method print,view-method Rm,view,workFlow,character-method
#' show,view-method identifier,view-method identifier,NULL-method
#' xyplot,formula,view-method xyplot,view,missing-method
#' @docType class
#' @usage
#' view(workflow, ID=paste("viewRef", guid(), sep="_"),
#'      name="default", data, action)
#' 
#' parent(object)
#' 
#' Data(object)
#' 
#' action(object)
#' 
#' alias(object, \dots)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param object An object of class \code{view} or one of its subclasses.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action References to the data and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @param \dots Further arguments that get passed to the generic.
#' @return
#' 
#' A reference to the view that is created inside the \code{\link{workFlow}}
#' environment as a side effect of calling the constructor.
#' 
#' The parent view (i.e., the view based on which the current view was created)
#' for the parent method.
#' @section Objects from the Class:
#' 
#' Objects should be created using the constructor \code{view}, which also
#' assigns the view to a \code{\link{workFlow}} object.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the
#' \code{\link[flowCore:actionItem-class]{actionItem}} that generated
#' the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' @slot alias Object of class \code{"fcAliasReference"}. A
#' reference to the alias table.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view. See
#' \code{\link{gateView}} for details on copying and subsetting of the raw
#' data in the context of gating.
#' 
#' @section Methods:
#' \describe{
#'   \item{action}{\code{signature(object = "view")}: Accessor for the
#'     \code{action} slot. Note that this returns the actual
#'     \code{\link[flowCore:actionItem-class]{actionItem}} object, i.e.,
#'     the reference gets resolved. } 
#'   
#'   \item{Data}{\code{signature(object = "view")}:  Accessor for the
#'     \code{data} slot. Note that this returns the actual data  object,
#'     i.e., the reference gets resolved. }
#'   
#'   \item{names}{\code{signature(x = "view")}: Accessor to the
#'     \code{name} slot. }
#'   
#'   \item{alias}{\code{signature(object = "view")}: Get the alias table
#'     from a \code{view}. }
#'   
#'   \item{parent}{\code{signature(object = "view")}: The parent view,
#'     i.e., the view based on which the current view was created. }
#'   
#'   \item{print}{\code{signature(x = "view")}: Print details about the
#'     object. }
#'   
#'   \item{Rm}{\code{signature(symbol = "view", envir = "workFlow",
#'                             subSymbol = "character")}: Remove a \code{view} from a
#'     \code{\link{workFlow}}. This method is recursive and will also
#'     remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItems}}. }
#'   
#'   \item{show}{\code{signature(object = "view")}: Print details about the
#'     object. }
#'   
#'   \item{xyplot}{\code{signature(x = "formula", data = "view")}: Plot
#'     the data underlying the view. }
#'   
#'   \item{xyplot}{\code{signature(x = "view", data = "missing")}: Plot
#'     the data underlying the view. }
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{gateView}},
#' \code{\linkS4class{transformView}}, \code{\linkS4class{compensateView}},
#' \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export 
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
#' @export
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
#' Class "gateView"
#' 
#' 
#' Class and method to capture the result of gating operations in a flow
#' cytometry workflow.
#' 
#' 
#' \code{gateViews} provide a means to bind the results of gating operations in
#' a workflow. Each \code{gateView} represents one of the populations that
#' arise from the gating.
#' \code{\link[flowCore:logicalFilterResult-class]{logicalFilterResults}}
#' create two \code{gateViews} (events in the gate and events not in the gate),
#' \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResults}} one
#' \code{view} for each population. See the documentation of the parent class
#' \code{\link{view}} for more details.
#' 
#' @name gateView-class
#' @aliases gateView-class gateView Rm,gateView,workFlow,character-method
#' summary,gateView-method xyplot,formula,gateView-method
#' @docType class
#' @usage 
#' gateView(workflow, ID=paste("gateViewRef", guid(), sep="_"),
#'          name="default", action, data, indices, 
#'          filterResult, frEntry)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action,filterResult References to the data,
#' \code{\link{filterResult}}, and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @param indices A logical vector of indices in the parent data.
#' @param frEntry A character vector indicating the name of the population in
#' the \code{\link{filterResult}}.
#' @return
#' 
#' A reference to the \code{gateView} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{gateView} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{gateView} from a \code{\link{filter}} object and directly assigns it
#' to a \code{\link{workFlow}}. Alternatively, one can use the \code{gateView}
#' constructor function for more programmatic access.
#' 
#' @slot indices Object of class \code{"logical"}. The indices
#' in the parent data for events that are within the filter.
#' @slot filterResult Object of class
#' \code{"fcFilterResultReference"}. A reference to the outcome of
#' the filtering operation.
#' @slot frEntry Object of class \code{"character"} The
#' population in the \code{\link{filterResult}} that corresponds to the
#' current view. See details for further explanation.
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the
#' \code{\link[flowCore:actionItem-class]{actionItem}} that generated
#' the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the
#' \code{\link[flowCore:workFlow-class]{workFlow}}.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view. Subsets of
#' the data are only generated when a a further action is invoked on
#' a particular \code{gateView}. Summary statistics about the view
#' can be acquired through the usual process of summarizing
#' \code{\link[flowCore:filterResult-class]{filterResults}}.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{view}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{Rm}{\code{signature(symbol = "gateView", envir = "workFlow",
#'                             subSymbol = "character")}: Remove a \code{gateView} from a
#'     \code{\link{workFlow}}. This method is recursive and will also
#'     remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItem}}s.  }
#'   
#'   \item{summary}{\code{signature(x = "formula", data = "gateView")}:
#'       Summarize the gating operation. }
#'   
#'   \item{xyplot}{\code{signature(x = "formula", data = "gateView")}:
#'       Plot the data of the \code{gateView} along with the gate. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{view}},
#' \code{\linkS4class{transformView}}, \code{\linkS4class{compensateView}},
#' \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' Class "transformView"
#' 
#' 
#' Class and method to capture the result of transformation operations in a
#' flow cytometry workflow.
#' 
#' 
#' @name transformView-class
#' @aliases transformView-class transformView
#' Rm,transformView,workFlow,character-method
#' @docType class
#' @usage 
#' transformView(workflow, ID=paste("transViewRef", guid(), sep="_"),
#'               name="default", action, data)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action References to the data and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @return
#' 
#' A reference to the \code{transformView} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{transformView} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{transformView} from a \code{\link{transform}} object and directly
#' assigns it to a \code{\link{workFlow}}. Alternatively, one can use the
#' \code{transformView} constructor function for more programmatic access.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the \code{\link[flowCore:actionItem-class]{actionItem}} 
#' that generated the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{view}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{Rm}{\code{signature(symbol = "transformView", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{transformView} from a \code{\link{workFlow}}. This method is
#'     recursive and will also remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItems}}. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{view}},
#' \code{\linkS4class{gateView}}, \code{\linkS4class{compensateView}},
#' \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export 
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
#' @export
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
#' Class "compensateView"
#' 
#' 
#' Class and method to capture the result of compensation operations in a flow
#' cytometry workflow.
#' 
#' 
#' @name compensateView-class
#' @aliases compensateView-class compensateView
#' Rm,compensateView,workFlow,character-method
#' @docType class
#' @usage
#' compensateView(workflow, ID=paste("compViewRef", guid(), sep="_"),
#'                name="default", action, data)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action References to the data and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @return
#' 
#' A reference to the \code{compensateView} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{compensateView} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{compensateView} from a \code{\link{compensation}} object and directly
#' assigns it to a \code{\link{workFlow}}. Alternatively, one can use the
#' \code{compensateView} constructor function for more programmatic access.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the \code{\link[flowCore:actionItem-class]{actionItem}}
#' that generated the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{view}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{Rm}{\code{signature(symbol = "compensateView", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{compensateView} from a \code{\link{workFlow}}. This method is
#'     recursive and will also remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItem}}s. }
#'   
#' }
#' 
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{view}},
#' \code{\linkS4class{gateView}}, \code{\linkS4class{transformView}},
#' \code{\linkS4class{normalizeView}}, \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' @export
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
#' Class "normalizeView"
#' 
#' 
#' Class and method to capture the result of normalization operations in a flow
#' cytometry workflow.
#' 
#' 
#' @name normalizeView-class
#' @aliases normalizeView-class normalizeView
#' Rm,normalizeView,workFlow,character-method
#' @docType class
#' @usage 
#' normalizeView(workflow, ID=paste("normViewRef", guid(), sep="_"),
#'               name="default", action, data)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action References to the data and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @return
#' 
#' A reference to the \code{normalizeView} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{normalizeView} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{normalizeView} from a \code{\link{normalization}} object and directly
#' assigns it to a \code{\link{workFlow}}. Alternatively, one can use the
#' \code{normalizeView} constructor function for more programmatic access.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the \code{\link[flowCore:actionItem-class]{actionItem}} 
#' that generated the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{view}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{Rm}{\code{signature(symbol = "normalizeView", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{normalizeView} from a \code{\link{workFlow}}. This method is
#'     recursive and will also remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItem}}s. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{view}},
#' \code{\linkS4class{gateView}}, \code{\linkS4class{transformView}},
#' \code{\linkS4class{compensateView}}, \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' @export
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
#' @export
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
#' Class "subsettingView"
#' 
#' 
#' Class and method to capture the result of subsetting operations in a flow
#' cytometry workflow.
#' 
#' 
#' @name subsettingView-class
#' @aliases subsettingView-class subsettingView
#' Rm,subsettingView,workFlow,character-method
#' @docType class
#' @usage 
#' subsettingView(workflow, ID=paste("subViewRef", guid(), sep="_"),
#'                name="default", action, data)
#' @param workflow An object of class \code{\link{workFlow}} for which a view
#' is to be created.
#' @param ID A unique identifier of the view, most likely created by using the
#' internal \code{guid} function.
#' @param name A more human-readable name of the view.
#' @param data,action References to the data and
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects, respectively.
#' @return
#' 
#' A reference to the \code{subsettingView} that is created inside the
#' \code{\link{workFlow}} environment as a side effect of calling the
#' \code{add} method.
#' 
#' A \code{subsettingView} object for the constructor.
#' @section Objects from the Class:
#' 
#' Objects should be created using the \code{add} method, which creates a
#' \code{subsettingView} from a \code{\link{subsetting}} object and directly
#' assigns it to a \code{\link{workFlow}}. Alternatively, one can use the
#' \code{subsettingView} constructor function for more programmatic access.
#' 
#' @slot ID Object of class \code{"character"}. A unique
#' identifier for the view.
#' @slot name Object of class \code{"character"}. A more
#' human-readable name.
#' @slot action Object of class \code{"fcActionReference"}. A
#' reference to the
#' \code{\link[flowCore:actionItem-class]{actionItem}} that generated
#' the view.
#' @slot env Object of class \code{"environment"}. The
#' evaluation environment in the \code{\link{workFlow}}.
#' @slot data Object of class \code{"fcDataReference"} A
#' reference to the data that is associated to the view.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{view}"}, directly.
#' 
#' @section Methods:
#' \describe{
#'   
#'   \item{Rm}{\code{signature(symbol = "subsettingView", envir =
#'                               "workFlow", subSymbol = "character")}: Remove a
#'     \code{subsettingView} from a \code{\link{workFlow}}. This method is
#'     recursive and will also remove all dependent \code{views} and
#'     \code{\link[flowCore:actionItem-class]{actionItem}}s. }
#'   
#' }
#' 
#' @author Florian Hahne
#' @seealso
#' 
#' \code{\linkS4class{workFlow}}, \code{\linkS4class{view}},
#' \code{\linkS4class{gateView}}, \code{\linkS4class{transformView}},
#' \code{\linkS4class{compensateView}}, \code{\linkS4class{actionItem}}
#' @keywords classes
#' @examples
#' 
#' showClass("view")
#' 
#' @export
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
#' @export
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
#' @export
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


#' @export
setMethod("add",
          signature=signature(wf="workFlow", action="character"),
          definition=function(wf, action, ...)
      {
          action <- subsetting(action)
          add(wf, action,...)
      })

#' @export
setMethod("add",
          signature=signature(wf="workFlow", action="numeric"),
          definition=function(wf, action, ...)
      {
          action <- subsetting(action)
          add(wf, action,...)
      })

#' @export
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
## objects which can be evaluated to retrieve the corresponding columns from the
## data frame
## ---------------------------------------------------------------------------
#' Class "unitytransform"
#' 
#' Unity transform class transforms parameters names provided as characters
#' into unity transform objects which can be evaluated to retrieve the
#' corresponding columns from the data frame
#' 
#' 
#' @name unitytransform-class
#' @aliases unitytransform-class unitytransform show,unitytransform-method
#' eval,unitytransform,missing-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{unitytransform(parameters,transformationId)}.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot parameters Object of class \code{"character"} -- the flow
#' parameters to be transformed.
#' @slot transformationId Object of class \code{"character"} -- a unique Id to
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso dg1polynomial, ratio
#' @family mathematical transform classes
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   un1<-unitytransform(c("FSC-H","SSC-H"),transformationId="un1")
#'   transOut<-eval(un1)(exprs(dat))
#' 
#' @export 
setClass("unitytransform",
	 contains="transform",
	 representation=representation(parameters="character"))

#' @export
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
#' Class "dg1polynomial"
#' 
#' dg1polynomial allows for scaling,linear combination and translation within a
#' single transformation defined by the function
#' \deqn{ f(parameter_1,...,parameter_n,a_1,...,a_n,b) = b + \Sigma_{i=1}^n
#' a_i*parameter_i }
#' 
#' 
#' @name dg1polynomial-class
#' @aliases dg1polynomial-class dg1polynomial eval,dg1polynomial,missing-method
#' initialize,dg1polynomial-method parameters<-,dg1polynomial,character-method
#' parameters<-,dg1polynomial,parameters-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column.(See example below)
#' @section Objects from the Class: Objects can be created by using the
#' constructor \code{dg1polynomial(parameter,a,b,transformationId)}.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot parameters Object of class \code{"parameters"} --the flow parameters
#' that are to be transformed.
#' @slot a Object of class \code{"numeric"} -- coefficients of length equal
#' to the number of flow parameters.
#' @slot b Object of class \code{"numeric"} -- coefficient of length 1 that
#' performs the translation.
#' @slot transformationId Object of class \code{"character"} unique ID to
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso ratio,quadratic,squareroot
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   dg1<-dg1polynomial(c("FSC-H","SSC-H"),a=c(1,2),b=1,transformationId="dg1")
#'   transOut<-eval(dg1)(exprs(dat))
#' 
#' @export
setClass("dg1polynomial", 		
         contains="transform",
         representation=representation(parameters="parameters",
                                       a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=new("parameters"),
                             a=1,
                             b=1))

#' @export
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
#' Class "ratio"
#' 
#' ratio transform calculates the ratio of two parameters defined by the
#' function \deqn{f(parameter_1,parameter_2)=\frac{parameter_1}{parameter_2}}
#' 
#' 
#' @name ratio-class
#' @aliases ratio-class ratio eval,ratio,missing-method initialize,ratio-method
#' @docType class
#' @note The ratio transformation object can be evaluated using the eval method
#' by passing the data frame as an argument.The transformed parameters are
#' returned as matrix with one column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{ratio(parameter1,parameter2,transformationId) }.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot numerator Object of class \code{"transformation"} -- flow parameter
#' to be transformed
#' @slot denominator Object of class \code{"transformation"} -- flow parameter
#' to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso  dg1polynomial,quadratic,squareroot
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   rat1<-ratio("FSC-H","SSC-H",transformationId="rat1")
#'   transOut<-eval(rat1)(exprs(dat))
#' 
#' @export
setClass("ratio",
         contains="transform",
         representation(numerator="transformation",
                        denominator="transformation"),
	 prototype=prototype(numerator=unitytransform(),
                             denominator=unitytransform()))

#' @export
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
#' Class "quadratic"
#' 
#' Quadratic transform class which represents a transformation defined by the 
#' function \deqn{f(parameter,a)=a*parameter^2}
#' 
#' 
#' @name quadratic-class
#' @aliases quadratic quadratic-class quadratic eval,quadratic,missing-method
#' @docType class
#' @note The quadratic transformation object can be evaluated using the eval
#' method by passing the data frame as an argument.The transformed parameters
#' are returned as a column vector. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{quadratic(parameters,a,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- non-zero multiplicative 
#' constant.
#' @slot parameters Object of class \code{"transformation"} -- flow 
#' parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique 
#' ID to reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", 
#' distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform",
#' distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform",
#' distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso dg1polynomial,ratio,squareroot
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   quad1<-quadratic(parameters="FSC-H",a=2,transformationId="quad1")
#'   transOut<-eval(quad1)(exprs(dat))
#' 
#' @export
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
             msg <- c(msg, "'a' should be non-zero")
         msg
     })

#' @export
quadratic <- function(parameters="NULL", a=1,
                      transformationId="defaultQuadraticTransform")
    new("quadratic",parameters=parameters,a=a,
        transformationId=transformationId)

          

## ===========================================================================
## Squareroot transformation
## ---------------------------------------------------------------------------
#' Class "squareroot"
#' 
#' Square root transform class, which represents a transformation defined by the 
#' function \deqn{f(parameter,a)= \sqrt{ |{\frac{parameter}{a}|}}}
#' 
#' 
#' @name squareroot-class
#' @aliases squareroot-class squareroot squareroot eval,squareroot,missing-method
#' @docType class
#' @note The squareroot transformation object can be evaluated using the eval
#' method by passing the data frame as an argument.The transformed parameters
#' are returned as a column vector. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{squareroot(parameters,a,transformationId)}
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @slot .Data Object of class \code{"function"}
#' @slot a Object of class \code{"numeric"} -- non-zero multiplicative 
#' constant
#' @slot parameters Object of class \code{"transformation"} -- flow 
#' parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique 
#' ID to reference the transformation.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso dg1polynomial, ratio, quadratic
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   sqrt1<-squareroot(parameters="FSC-H",a=2,transformationId="sqrt1")
#'   transOut<-eval(sqrt1)(exprs(dat))
#' 
#' @export
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
             msg <- c(msg, "Coefficien> t should be non-zero")
         msg
     })

#' @export
squareroot <- function(parameters, a=1,
                       transformationId="defaultSquarerootTransform")
    new("squareroot", parameters=parameters, a=a,
        transformationId=transformationId)



## ===========================================================================
##  Logarithmic Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
#' Class "logarithm"
#' 
#' Logartithmic transform class, which represents a transformation defined by
#' the function
#' 
#' \deqn{f(parameter,a,b)= ln(a*prarameter)*b ~~~~a*parameter>0} \deqn{0
#' ~~~~a*parameter<=0}
#' 
#' 
#' @name logarithm-class
#' @aliases logarithm-class logarithm eval,logarithm,missing-method
#' @docType class
#' @note The logarithm transformation object can be evaluated using the eval
#' method by passing the data frame as an argument.The transformed parameters
#' are returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{logarithm(parameters,a,b,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}
#' @slot a Object of class \code{"numeric"} -- non-zero multiplicative constant.
#' @slot b Object of class \code{"numeric"} -- non-zero multiplicative constant.
#' @slot parameters Object of class \code{"transformation"} -- flow parameters to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso exponential, quadratic
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'  dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   lg1<-logarithm(parameters="FSC-H",a=2,b=1,transformationId="lg1")
#'   transOut<-eval(lg1)(exprs(dat))
#' 
#' @export
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
             msg <- c(msg, "'a' should be a non-zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non-zero number")
         msg
     })

#' @export
logarithm <- function(parameters, a=1, b=1,
                      transformationId="defaultLogarithmTransform")
    new("logarithm", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Exponential Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
#' Class "exponential"
#' 
#' Exponential transform class, which represents a transformation given by the 
#' function \deqn{f(parameter,a,b)=e^{parameter/b}*\frac{1}{a}}
#' 
#' 
#' @name exponential-class
#' @aliases exponential-class exponential eval,exponential,missing-method
#' @docType class
#' @note The exponential transformation object can be evaluated using the eval
#' method by passing the data frame as an argument.The transformed parameters
#' are returned as a matrix with a single column
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor\code{exponential(parameters,a,b)}.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- non-zero constant.
#' @slot b Object of class \code{"numeric"}- non-zero constant.
#' @slot parameters Object of class \code{"transformation"} -- flow 
#' parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- 
#' unique ID to reference the transformation
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#'
#' @author Gopalakrishnan N, F.Hahne
#' @seealso logarithm
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'  dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   exp1<-exponential(parameters="FSC-H",a=1,b=37,transformationId="exp1")
#'   transOut<-eval(exp1)(exprs(dat))
#' 
#' @export
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
             msg<-c(msg,"'a' should be a non-zero number")
         if(object@b==0)
             msg <- c(msg,"'b' should be a non-zero number")
         msg  
     })

#' @export
exponential <- function(parameters, a=1, b=1,
                        transformationId="defaultExponentialTransformation")
    new("exponential", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Inverse hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
#' Class "asinht"
#' 
#' Inverse hyperbolic sine transform class, which represents a transformation 
#' defined by the function: 
#' \deqn{f(parameter,a,b)=sinh^{-1}(a*parameter)*b}
#' This definition is such that it can function as an inverse of 
#' \code{\linkS4class{sinht}} using the same definitions of the constants a
#' and b.
#' 
#' @name asinht-class
#' @aliases asinht-class asinht eval,asinht,missing-method
#' @docType class
#' @note The inverse hyperbolic sin transformation object can be evaluated
#' using the eval method by passing the data frame as an argument.The
#' transformed parameters are returned as a matrix with a single column. (See
#' example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{asinht(parameter,a,b,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- non-zero constant.
#' @slot b Object of class \code{"numeric"} -- non-zero constant.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter
#' to be transformed
#' @slot transformationId Object of class \code{"character"} -- unique ID to 
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso sinht
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'  dat <- read.FCS(system.file("extdata","0877408774.B08",  package="flowCore"))
#'   asinh1<-asinht(parameters="FSC-H",a=2,b=1,transformationId="asinH1")
#'   transOut<-eval(asinh1)(exprs(dat))
#' 
#' @export
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
             msg <- c(msg, "'a' should be a non-zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non-zero number")
         msg
     })

#' @export
asinht <- function(parameters="NULL", a=1, b=1,
                   transformationId="defaultAsinhTransform")
    new("asinht", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
#' Class "sinht"
#' 
#' Hyperbolic sin transform class, which represents a transformation 
#' defined by the function: 
#' \deqn{f(parameter,a,b)=sinh(parameter/b)/a} 
#' This definition is such that it can function as an inverse of 
#' \code{\linkS4class{asinht}} using the same definitions of the constants a
#' and b.
#' 
#' @name sinht-class
#' @aliases sinht-class sinht eval,sinht,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column.(See example below)
#' 
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{sinht(parameter,a,b,transformationId)}.
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- non-zero constant.
#' @slot b Object of class \code{"numeric"} -- non-zero constant.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter
#' to be transformed
#' @slot transformationId Object of class \code{"character"} -- unique ID to 
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso asinht
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'  dat <- read.FCS(system.file("extdata","0877408774.B08",  package="flowCore"))
#'  sinh1<-sinht(parameters="FSC-H",a=1,b=2000,transformationId="sinH1")
#'  transOut<-eval(sinh1)(exprs(dat))
#' 
#' @export
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
             msg <- c(msg, "'a' should be a non-zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non-zero number")
         msg
     })

#' @export
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
#' Class asinhtGml2
#' 
#' Inverse hyperbolic sin transformation as parameterized in Gating-ML 2.0. 
#' 
#' asinhtGml2 is defined by the following function: 
#' \deqn{bound(f, boundMin, boundMax) = max(min(f,boundMax),boundMin))} where 
#' \deqn{f(parameter, T, M, A) = (asinh(parameter * sinh(M * ln(10)) / T) +A * ln(10)) / ((M + A) * ln(10))}
#' 
#' This transformation is equivalent to Logicle(T, 0, M, A) (i.e., with W=0).
#' It provides an inverse hyperbolic sine transformation that maps a data value
#' onto the interval [0,1] such that: 
#' \itemize{ 
#' \item The top of scale value (i.e., T ) is mapped to 1.  
#' \item Large data values are mapped to locations similar to an 
#' (M + A)-decade logarithmic scale.  
#' \item A decades of negative data are brought on scale.
#' }
#' 
#' In addition, if a boundary is defined by the boundMin and/or boundMax
#' parameters, then the result of this transformation is restricted to the
#' [boundMin,boundMax] interval. Specifically, should the result of the f
#' function be less than boundMin, then let the result of this transformation
#' be boundMin. Analogically, should the result of the f function be more than
#' boundMax, then let the result of this transformation be boundMax. The
#' boundMin parameter shall not be greater than the boundMax parameter.
#' 
#' 
#' @name asinhtGml2-class
#' @aliases asinhtGml2-class asinhtGml2 eval,asinhtGml2,missing-method
#' @docType class
#' @note The inverse hyperbolic sin transformation object can be evaluated
#' using the eval method by passing the data frame as an argument. The
#' transformed parameters are returned as a matrix with a single column. (See
#' example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{asinhtGml2(parameter, T, M, A, transformationId, boundMin, boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot T Object of class \code{numeric} -- positive constant (top of scale value).
#' @slot M Object of class \code{numeric} -- positive constant (desired number of decades).
#' @slot A Object of class \code{numeric} -- non-negative constant that is less than or equal 
#' to M (desired number of additional negative decades).
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{\linkS4class{singleParameterTransform}}, directly.
#' 
#' Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
#' 
#' Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
#' 
#' Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
#' 
#' @author Spidlen, J.
#' @seealso \code{\link{asinht}}, \code{\link{transform-class}},
#' \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' @keywords classes
#' @examples
#' 
#' myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myASinH1 <- asinhtGml2(parameters = "FSC-H", T = 1000, M = 4.5, 
#'     A = 0, transformationId="myASinH1")
#' transOut <- eval(myASinH1)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class logicletGml2
#' 
#' Logicle transformation as published by Moore and Parks.
#' 
#' logicletGml2 is defined by the
#' following function: 
#' \deqn{bound(logicle, boundMin, boundMax) = max(min(logicle,boundMax),boundMin))} 
#' where \deqn{logicle(x, T, W, M, A) = root(B(y, T, W, M, A) - x)} and \eqn{B} 
#' is a modified biexponential function: 
#' \deqn{B(y, T, W, M, A) = ae^{by} - ce^{-dy} - f} where 
#' \itemize{
#' \item x is the value that is being transformed (an FCS dimension value).
#' Typically, x is less than or equal to T, although the transformation
#' function is also defined for x greater than T.
#' \item y is the result of the transformation.
#' \item T is greater than zero and represents the top of
#' scale value.
#' \item M is greater than zero and represents the number of
#' decades that the true logarithmic scale approached at the high end of the
#' Logicle scale would cover in the plot range.
#' \item W is non-negative and not greater than half of M and represents the 
#' number of such decades in the approximately linear region. The choice of 
#' \eqn{W = M/2} specifies a scale that is essentially linear over the whole 
#' range except for a small region of large data values. For situations in which 
#' values of W approaching \eqn{M/2} might be chosen, ordinary linear display scales 
#' will usually be more appropriate. The choice of \eqn{W = 0} gives essentially the 
#' hyperbolic sine function.
#' \item A is the number of additional decades of negative data
#' values to be included. A shall be greater than or equal to \eqn{-W}, and
#' less than or equal to \eqn{M - 2W}
#' \item root is a standard root finding
#' algorithm (e.g., Newton's method) that finds y such as \eqn{B(y, T, W, M, A)
#' = x}.
#' } 
#' and \eqn{a}, \eqn{b}, \eqn{c}, \eqn{d} and \eqn{f} are defined by
#' means of \eqn{T}, \eqn{W}, \eqn{M}, \eqn{A}, \eqn{w}, \eqn{x0}, \eqn{x1},
#' \eqn{x2}, \eqn{ca} and \eqn{fa} as: 
#' \deqn{w = W/(M+A)} \deqn{x2 = A/(M+A)}
#' \deqn{x1 = x2 + w} 
#' \deqn{x0 = x2 + 2*w} 
#' \deqn{b = (M + A)*ln(10)} and
#' \eqn{d} is a constant so that \deqn{2*(ln(d) - ln(b)) + w*(d + b) = 0} given
#' \eqn{b} and \eqn{w}, and 
#' \deqn{ca = e^{x0*(b+d)}} 
#' \deqn{fa = e^{b*x1} - (ca/(e^{d*x1}))} 
#' \deqn{a = T / (e^b - fa - (ca/e^d)) } \deqn{c = ca * a}
#' \deqn{f = fa * a}
#' 
#' The Logicle scale is the inverse of a modified biexponential function. It
#' provides a Logicle display that maps scale values onto the \eqn{[0,1]}
#' interval such that the data value \eqn{T} is mapped to 1, large data values
#' are mapped to locations similar to an (M + A)-decade logarithmic scale, and
#' A decades of negative data are brought on scale. For implementation
#' purposes, it is recommended to follow guidance in Moore and Parks
#' publication.
#' 
#' In addition, if a boundary is defined by the boundMin and/or boundMax
#' parameters, then the result of this transformation is restricted to the
#' [boundMin,boundMax] interval. Specifically, should the result of the logicle
#' function be less than boundMin, then let the result of this transformation
#' be boundMin. Analogically, should the result of the logicle function be more
#' than boundMax, then let the result of this transformation be boundMax. The
#' boundMin parameter shall not be greater than the boundMax parameter.
#' 
#' 
#' @name logicletGml2-class
#' @aliases logicletGml2-class logicletGml2 eval,logicletGml2,missing-method
#' @docType class
#' @note Please note that \code{logicletGml2} and
#' \code{\link{logicleTransform}} are similar transformations; however, the
#' Gating-ML 2.0 compliant \code{logicletGml2} brings "reasonable" data values
#' to the scale of \eqn{[0,1]} while the \code{\link{logicleTransform}} scales
#' these values to \eqn{[0,M]}.
#' 
#' The logicle transformation object can be evaluated using the eval method by
#' passing the data frame as an argument. The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{logicletGml2(parameter, T, M, W, A, transformationId, boundMin,
#' boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot T Object of class \code{numeric} -- positive constant (top of scale value).
#' @slot M Object of class \code{numeric} -- positive constant (desired number of decades).
#' @slot W Object of class \code{numeric} -- non-negative constant that is not greater than half of M
#' (the number of such decades in the approximately linear region).
#' @slot A Object of class \code{numeric} -- a constant that is greater than or equal to -W, and also
#' less than or equal to M-2W. (A represents the number of additional decades of negative data values to 
#' be included.)
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{\linkS4class{singleParameterTransform}}, directly.
#' 
#' Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
#' 
#' Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
#' 
#' Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
#' 
#' @author Spidlen, J., Moore, W.
#' @seealso \code{\link{logicleTransform}}, \code{\link{transform-class}},
#' \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' 
#' Moore, WA and Parks, DR. Update for the logicle data scale including
#' operational code implementations. Cytometry A., 2012:81A(4):273-277.
#' 
#' Parks, DR and Roederer, M and Moore, WA. A new "Logicle" display method
#' avoids deceptive effects of logarithmic scaling for low signals and
#' compensated data. Cytometry A., 2006:69(6):541-551.
#' @keywords classes
#' @examples
#' 
#' myDataIn  <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myLogicle <- logicletGml2(parameters = "FSC-H", T = 1023, M = 4.5, 
#'     W = 0.5, A = 0, transformationId="myLogicle")
#' transOut  <- eval(myLogicle)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class hyperlogtGml2
#' 
#' Hyperlog transformation parameterized according to Gating-ML 2.0.
#' 
#' hyperlogtGml2 is defined by the following function: 
#' \deqn{bound(hyperlog, boundMin, boundMax) = max(min(hyperlog,boundMax),boundMin))} 
#' where \deqn{hyperlog(x, T, W, M, A) = root(EH(y, T, W, M, A) - x)} and 
#' \eqn{EH} is defined as: 
#' \deqn{EH(y, T, W, M, A) = ae^{by} + cy - f} where 
#' \itemize{ 
#' \item x is the value that is being
#' transformed (an FCS dimension value). Typically, x is less than or equal to
#' T, although the transformation function is also defined for x greater than
#' T.
#' \item y is the result of the transformation.
#' \item T is greater than zero and represents the top of scale value.
#' \item M is greater than zero and represents the number of decades that the 
#' true logarithmic scale approached at the high end of the Hyperlog scale would 
#' cover in the plot range.
#' \item W is positive and not greater than half of M and represents the number of 
#' such decades in the approximately linear region.
#' \item A is the number of additional decades of negative data values to be included. A
#' shall be greater than or equal to \eqn{-W}, and less than or equal to \eqn{M
#' - 2W}
#' \item root is a standard root finding algorithm (e.g., Newton's
#' method) that finds y such as \eqn{B(y, T, W, M, A) = x}. } and \eqn{a},
#' \eqn{b}, \eqn{c} and \eqn{f} are defined by means of \eqn{T}, \eqn{W},
#' \eqn{M}, \eqn{A}, \eqn{w}, \eqn{x0}, \eqn{x1}, \eqn{x2}, \eqn{e0}, \eqn{ca}
#' and \eqn{fa} as: 
#' \deqn{w = W/(M+A)} 
#' \deqn{x2 = A/(M+A)} 
#' \deqn{x1 = x2 + w}
#' \deqn{x0 = x2 + 2*w} 
#' \deqn{b = (M + A)*ln(10)} 
#' \deqn{e0 = e^{b*x0}} 
#' \deqn{ca= e0/w} 
#' \deqn{fa = e^{b*x1} + ca*x1} 
#' \deqn{a = T / (e^b + ca - fa)} 
#' \deqn{c = ca * a} 
#' \deqn{f = fa * a}
#' 
#' In addition, if a boundary is defined by the boundMin and/or boundMax
#' parameters, then the result of this transformation is restricted to the
#' [boundMin,boundMax] interval. Specifically, should the result of the
#' hyperlog function be less than boundMin, then let the result of this
#' transformation be boundMin. Analogically, should the result of the hyperlog
#' function be more than boundMax, then let the result of this transformation
#' be boundMax. The boundMin parameter shall not be greater than the boundMax
#' parameter.
#' 
#' 
#' @name hyperlogtGml2-class
#' @aliases hyperlogtGml2-class hyperlogtGml2 eval,hyperlogtGml2,missing-method
#' @docType class
#' @note That \code{hyperlogtGml2} transformation brings "reasonable" data
#' values to the scale of \eqn{[0,1]}.  The transformation is somewhat similar
#' to \code{\link{logicletGml2}}. (See Gating-ML 2.0 for detailed comparison)
#' 
#' The hyperlog transformation object can be evaluated using the eval method by
#' passing the data frame as an argument. The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{hyperlogtGml2(parameter, T, M, W, A, transformationId, boundMin,
#' boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot T Object of class \code{numeric} -- positive constant (top of scale value).
#' @slot M Object of class \code{numeric} -- positive constant (desired number of decades).
#' @slot W Object of class \code{numeric} -- positive constant that is not greater than half of M
#' (the number of such decades in the approximately linear region)
#' @slot A Object of class \code{numeric} -- a constant that is greater than or equal to -W, and also
#' less than or equal to M-2W. (A represents the number of additional decades of negative data values to 
#' be included.)
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{\linkS4class{singleParameterTransform}}, directly.
#' 
#' Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
#' 
#' Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
#' 
#' Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
#' 
#' @author Spidlen, J., Moore, W.
#' @seealso \code{\link{hyperlog}}, \code{\link{logicleTransform}},
#' \code{\link{transform-class}}, \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' @keywords classes
#' @examples
#' 
#' myDataIn  <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myHyperLg <- hyperlogtGml2(parameters = "FSC-H", T = 1023, M = 4.5, 
#'     W = 0.5, A = 0, transformationId="myHyperLg")
#' transOut  <- eval(myHyperLg)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class lintGml2
#' 
#' Linear transformation as parameterized in Gating-ML 2.0.
#' 
#' lintGml2 is defined by the following function: 
#' \deqn{bound(f, boundMin, boundMax) = max(min(f,boundMax),boundMin))} where 
#' \deqn{f(parameter, T, A) = (parameter + A) / (T + A)}
#' 
#' This transformation provides a linear display that maps scale values from
#' the \eqn{[-A, T]} interval to the \eqn{[0, 1]} interval.  However, it is
#' defined for all \eqn{x in R} including outside of the \eqn{[-A, T]}
#' interval.
#' 
#' In addition, if a boundary is defined by the boundMin and/or boundMax
#' parameters, then the result of this transformation is restricted to the
#' [boundMin,boundMax] interval. Specifically, should the result of the f
#' function be less than boundMin, then let the result of this transformation
#' be boundMin. Analogically, should the result of the f function be more than
#' boundMax, then let the result of this transformation be boundMax. The
#' boundMin parameter shall not be greater than the boundMax parameter.
#' 
#' 
#' @name lintGml2-class
#' @aliases lintGml2-class lintGml2 eval,lintGml2,missing-method
#' @docType class
#' @note The linear transformation object can be evaluated using the eval
#' method by passing the data frame as an argument. The transformed parameters
#' are returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{lintGml2(parameter, T, A, transformationId, boundMin, boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot T Object of class \code{numeric} -- positive constant (top of scale value).
#' @slot A Object of class \code{numeric} -- non-negative constant that is less than or equal
#' to T; it is determining the bottom end of the transformation.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{\linkS4class{singleParameterTransform}}, directly.
#' 
#' Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
#' 
#' Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
#' 
#' Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
#' 
#' @author Spidlen, J.
#' @seealso \code{\link{linearTransform}}, \code{\link{transform-class}},
#' \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' @keywords classes
#' @examples
#' 
#' myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myLinTr1 <- lintGml2(parameters = "FSC-H", T = 1000, A = 0, 
#'     transformationId="myLinTr1")
#' transOut <- eval(myLinTr1)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class logtGml2
#' 
#' Log transformation as parameterized in Gating-ML 2.0.
#' 
#' logtGml2 is defined by the following function: 
#' \deqn{bound(f, boundMin, boundMax) = max(min(f,boundMax),boundMin))} where 
#' \deqn{f(parameter, T, M) = (1/M) * log10(x/T) + 1}
#' 
#' This transformation provides a logarithmic display that maps scale values
#' from the \eqn{(0, T]} interval to the \eqn{(-Inf, 1]} interval such that the
#' data value T is mapped to 1 and M decades of data are mapped into the
#' interval.  Also, the limit for x going to 0 is -Inf.
#' 
#' In addition, if a boundary is defined by the boundMin and/or boundMax
#' parameters, then the result of this transformation is restricted to the
#' [boundMin,boundMax] interval. Specifically, should the result of the f
#' function be less than boundMin, then let the result of this transformation
#' be boundMin. Analogically, should the result of the f function be more than
#' boundMax, then let the result of this transformation be boundMax. The
#' boundMin parameter shall not be greater than the boundMax parameter.
#' 
#' 
#' @name logtGml2-class
#' @aliases logtGml2-class logtGml2 eval,logtGml2,missing-method
#' @docType class
#' @note The log transformation object can be evaluated using the eval method
#' by passing the data frame as an argument. The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{logtGml2(parameter, T, M, transformationId, boundMin, boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot T Object of class \code{numeric} -- positive constant (top of scale value).
#' @slot M Object of class \code{numeric} -- positive constant (number of decades).
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{\linkS4class{singleParameterTransform}}, directly.
#' 
#' Class \code{\linkS4class{transform}}, by class singleParameterTransform, distance 2.
#' 
#' Class \code{\linkS4class{transformation}}, by class singleParameterTransform, distance 3.
#' 
#' Class \code{\linkS4class{characterOrTransformation}}, by class singleParameterTransform, distance 4.
#' 
#' @author Spidlen, J.
#' @seealso \code{\link{logTransform}}, \code{\link{transform-class}},
#' \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' @keywords classes
#' @examples
#' 
#' myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myLogTr1 <- logtGml2(parameters = "FSC-H", T = 1023, M = 4.5, 
#'     transformationId="myLogTr1")
#' transOut <- eval(myLogTr1)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class "ratiotGml2"
#' 
#' Ratio transformation as parameterized in Gating-ML 2.0.
#' 
#' ratiotGml2 is defined by the following function: 
#' \deqn{bound(f, boundMin, boundMax) =
#' max(min(f,boundMax),boundMin))} where 
#' \deqn{f(p1, p2, A, B, C) = A * (p1 - B) / (p2 - C)}
#' 
#' If a boundary is defined by the boundMin and/or boundMax parameters, then
#' the result of this transformation is restricted to the [boundMin,boundMax]
#' interval. Specifically, should the result of the f function be less than
#' boundMin, then let the result of this transformation be boundMin.
#' Analogically, should the result of the f function be more than boundMax,
#' then let the result of this transformation be boundMax. The boundMin
#' parameter shall not be greater than the boundMax parameter.
#' 
#' 
#' @name ratiotGml2-class
#' @aliases ratiotGml2-class ratiotGml2 eval,ratiotGml2,missing-method
#' initialize,ratiotGml2-method parameters,ratiotGml2-method
#' @docType class
#' @note The ratiotGml2 transformation object can be evaluated using the eval
#' method by passing the data frame as an argument. The transformed parameters
#' are returned as matrix with one column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' 
#' \code{ratiotGml2(p1, p2, A, B, C, transformationId, boundMin, boundMax)}
#' 
#' @slot .Data Object of class \code{function}.
#' @slot numerator Object of class \code{"transformation"} -- flow parameter to be 
#' used as numerator in the transformation function.
#' @slot denominator Object of class \code{"transformation"} -- flow parameter to be 
#' used as denominator in the transformation function.
#' @slot pA Object of class \code{numeric} constant A.
#' @slot pB Object of class \code{numeric} constant B.
#' @slot pC Object of class \code{numeric} constant C.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference 
#' the transformation.
#' @slot boundMin Object of class \code{numeric} -- lower bound of the transformation, default -Inf.
#' @slot boundMax Object of class \code{numeric} -- upper bound of the transformation, default Inf.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author Spidlen, J.
#' @seealso \code{\link{ratio}}, \code{\link{transform-class}},
#' \code{\link{transform}}
#' @family mathematical transform classes
#' @references Gating-ML 2.0: International Society for Advancement of
#' Cytometry (ISAC) standard for representing gating descriptions in flow
#' cytometry. \url{http://flowcyt.sourceforge.net/gating/20141009.pdf}
#' @keywords classes
#' @examples
#' 
#' myDataIn <- read.FCS(system.file("extdata", "0877408774.B08", 
#'     package="flowCore"))
#' myRatioT <- ratiotGml2("FSC-H", "SSC-H", pA = 2, pB = 3, 
#'     pC = -10, transformationId = "myRatioT")
#' transOut <- eval(myRatioT)(exprs(myDataIn))
#' 
#' @export
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

#' @export
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
#' Class "hyperlog"
#' 
#' Hyperlog transformation of a parameter is defined by the function
#' \deqn{f(parameter,a,b)=root{EH(y,a,b)-parameter}}
#' where EH is a function defined by \deqn{EH(y,a,b) = 10^{(\frac{y}{a})} +
#' \frac{b*y}{a}-1, y>=0}
#' \deqn{EH(y,a,b)= -10^{(\frac{-y}{a})} + \frac{b*y}{a}+1, y<0}
#' 
#' 
#' @name hyperlog-class
#' @aliases hyperlog-class hyperlog eval,hyperlog,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{hyperlog(parameter,a,b,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- numeric constant
#' treater than zero.
#' @slot b Object of class \code{"numeric"} numeric constant greater than zero.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be 
#' transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to 
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso EHtrans
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'   package="flowCore"))
#'   hlog1<-hyperlog("FSC-H",a=1,b=1,transformationId="hlog1")
#'   transOut<-eval(hlog1)(exprs(dat))
#' 
#' @export
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

#' @export
hyperlog <- function(parameters="NULL", a=1, b=1,
                     transformationId="defaultHyperlogTransform")
    new("hyperlog", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  EH Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
#' Class "EHtrans"
#' 
#' EH transformation of a parameter is defined by the function
#' \deqn{EH(parameter,a,b)= 10^{(\frac{parameter}{a})} +
#' \frac{b*parameter}{a}-1, parameter>=0}
#' \deqn{-10^{(\frac{-parameter}{a})} + \frac{b*parameter}{a}+1, parameter<0}
#' 
#' 
#' @name EHtrans-class
#' @aliases EHtrans-class EHtrans eval,EHtrans,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor \code{EHtrans(parameters,a,b,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot a Object of class \code{"numeric"} -- numeric constant greater than zero.
#' @slot b Object of class \code{"numeric"} -- numeric constant greater than zero.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be 
#' transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso hyperlog
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry V 1.5
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",
#'                   package="flowCore"))
#'   eh1<-EHtrans("FSC-H",a=1250,b=4,transformationId="eh1")
#'   transOut<-eval(eh1)(exprs(dat))
#' 
#' @export
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

#' @export
EHtrans <- function(parameters, a=1, b=1,
                    transformationId="defaultEHtransTransform")
    new("EHtrans", parameters=parameters, a=a, b=b,
        transformationId=transformationId)
      
          

## ===========================================================================
##  Splitscale Transformation 
## ---------------------------------------------------------------------------
#' Class "splitscale"
#' 
#' The split scale transformation class defines a transformation that has a
#' logarithmic scale at high values and a linear scale at low values. The
#' transition points are chosen so that the slope of the transformation is
#' continuous at the transition points.
#' 
#' The split scale transformation is defined by the function
#' 
#' \deqn{f(parameter,r,maxValue,transitionChannel) = a*parameter+ b, parameter<=t}
#' \deqn{(parameter,r,maxValue,transitionChannel) = log_{10}(c*parameter)*\frac{r}{d}, parameter > t } where,
#' \deqn{b=\frac{transitionChannel}{2}}
#' \deqn{d=\frac{2*log_{10}(e)*r}{transitionChannel} + log_{10}(maxValue) }
#' \deqn{t=10^{log_{10}t}} \deqn{a= \frac{transitionChannel}{2*t}}
#' \deqn{log_{10}ct=\frac{(a*t+b)*d}{r}} \deqn{c=10^{log_{10}ct}}
#' 
#' 
#' @name splitscale-class
#' @aliases splitscale-class splitscale eval,splitscale,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' \code{splitscale(parameters,r,maxValue,transitionChannel,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot r Object of class \code{"numeric"} -- a positive value indicating the range of the logarithmic 
#' part of the display.
#' @slot maxValue Object of class \code{"numeric"} -- a positive value indicating the maximum value the transformation
#' is applied to.
#' @slot transitionChannel Object of class \code{"numeric"} -- non negative value that indicates where to 
#' split the linear vs. logarithmic transformation.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N, F.Hahne
#' @seealso invsplitscale
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
#'   sp1<-splitscale("FSC-H",r=768,maxValue=10000,transitionChannel=256)
#'   transOut<-eval(sp1)(exprs(dat))
#' 
#' @export
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

#' @export
splitscale <- function(parameters="NULL", r=1, maxValue=1, transitionChannel=4,
                       transformationId="defaultSplitscaleTransform")
    new("splitscale",
        parameters=parameters, r=r, maxValue=maxValue,
        transitionChannel=transitionChannel,
        transformationId=transformationId)



## ===========================================================================
##  Inverse Splitscale Transformation 
## ---------------------------------------------------------------------------
#' Class "invsplitscale"
#' 
#' As its name suggests, the inverse split scale transformation class represents
#' the inverse transformation of a split scale transformation that has a logarithmic scale at 
#' high values and a linear scale at low values.
#' 
#' The inverse split scale transformation is defined by the function
#' \deqn{f(parameter,r,maxValue,transitionChannel)  \frac{(parameter-b)}{a}, parameter<=t*a + b}
#' \deqn{f(parameter,r,maxValue,transitionChannel) = \frac{10^{parameter*\frac{d}{r}}}{c}, parameter > t*a+b }
#' where 
#' \deqn{b=\frac{transitionChannel}{2}}
#' \deqn{d=\frac{2*log_{10}(e)*r}{transitionChannel} + log_{10}(maxValue) }
#' \deqn{t=10^{log_{10}t}} \deqn{a= \frac{transitionChannel}{2*t}}
#' \deqn{log_{10}ct=\frac{(a*t+b)*d}{r}} \deqn{c=10^{log_{10}ct}}
#' 
#' 
#' @name invsplitscale-class
#' @aliases invsplitscale-class invsplitscale eval,invsplitscale,missing-method
#' @docType class
#' @note The transformation object can be evaluated using the eval method by
#' passing the data frame as an argument.The transformed parameters are
#' returned as a matrix with a single column. (See example below)
#' @section Objects from the Class: Objects can be created by calls to the
#' constructor
#' \code{invsplitscale(parameters,r,maxValue,transitionChannel,transformationId)}
#' 
#' @slot .Data Object of class \code{"function"}.
#' @slot r Object of class \code{"numeric"} -- a positive value indicating
#' the range of the logarithmic part of the dispmlay.
#' @slot maxValue Object of class \code{"numeric"} -- a positive value 
#' indicating the maximum value the transformation is applied to.
#' @slot transitionChannel Object of class \code{"numeric"} -- non negative 
#' value that indicates where to split the linear vs. logarithmic transformation.
#' @slot parameters Object of class \code{"transformation"} -- flow parameter
#' to be transformed.
#' @slot transformationId Object of class \code{"character"} -- unique ID to
#' reference the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{singleParameterTransform}"}, directly.
#' 
#' Class \code{"\linkS4class{transform}"}, by class "singleParameterTransform", distance 2.
#' 
#' Class \code{"\linkS4class{transformation}"}, by class "singleParameterTransform", distance 3.
#' 
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "singleParameterTransform", distance 4.
#' 
#' @author Gopalakrishnan N,F.Hahne
#' @seealso splitscale
#' @family mathematical transform classes
#' @references Gating-ML Candidate Recommendation for Gating Description in
#' Flow Cytometry
#' @keywords classes
#' @examples
#' 
#'   dat <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
#'   sp1<-invsplitscale("FSC-H",r=512,maxValue=2000,transitionChannel=512)
#'   transOut<-eval(sp1)(exprs(dat))
#' 
#' @export
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
       
#' @export
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
#' Class "transformReference"
#' 
#' Class allowing for reference of transforms, for instance as parameters.
#' 
#' 
#' @name transformReference-class
#' @aliases transformReference-class transformReference
#' parameters,transformReference-method eval,transformReference,missing-method
#' @docType class
#' @section Objects from the Class: Objects will be created internally whenever
#' necessary and this should not be of any concern to the user.
#' 
#' @slot .Data The list of references.
#' @slot searchEnv The environment into which the reference points.
#' @slot transformationId The name of the transformation.
#' 
#' @section Extends:
#' Class \code{"\linkS4class{transform}"}, directly.
#' Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
#' Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
#' 
#' @author N. Gopalakrishnan
#' @keywords classes
#'
#' @export 
setClass("transformReference",
         contains="transform",
         representation(searchEnv="environment"))

#' @export
transformReference <- function(referenceId="defaultTransformReference",
                               searchEnv)
    new("transformReference",
        transformationId=referenceId, searchEnv=searchEnv)
    
