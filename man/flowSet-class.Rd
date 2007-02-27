\name{flowSet-class}
\docType{class}
\alias{flowSet}
\alias{flowSet-class}
\alias{[,flowSet-method}
\alias{[[,flowSet-method}
\alias{[[<-,flowSet-method}
\alias{Subset,flowSet,ANY}
\alias{fsApply,flowSet-method}
\alias{show,flowSet-method}
\alias{length,flowSet-method}
\alias{pData,flowSet-method}
\alias{colnames,flowSet-method}
\alias{colnames<-,flowSet-method}
\alias{sampleNames,flowSet-method}
\alias{phenoData,flowSet-method}
\alias{phenoData<-,flowSet,ANY}
\alias{phenoData<-,flowSet,phenoData-method}
\alias{keyword,flowSet,list-method}
\alias{compensate,flowSet,matrix-method}
\alias{spillover,flowSet-method}
\alias{transform,flowSet-method}
\alias{split, flowSet,ANY}
\alias{rbind2, flowSet,missing}
\alias{rbind2, flowSet,flowSet-method}


\title{'flowSet': a class for storing flow cytometry raw data from
quantitative cell-based assays}

\description{This class is a container for a set of
  \code{\link[flowCore]{flowFrame}} objects}


\section{Creating Objects}{
Objects can be created using the function
\code{\link[flowCore]{read.flowSet}} or via\cr
\code{  new('flowSet',}\cr
\code{    frames = ....,   # environment with flowFrames}\cr
\code{    phenoData = .... # object of class AnnotatedDataFrame}\cr
\code{    colnames = ....  # object of class character}\cr
\code{  )}\cr
}

\section{Slots}{\describe{
\item{\code{frames}:}{An \code{\link[base:environment]{environment}}
containing one or more \code{\link[flowCore:flowFrame-class]{flowFrame}} objects.}
\item{\code{phenoData}:}{A
\code{\link[Biobase:AnnotatedDataFrame]{AnnotatedDataFrame}}. Each row
corresponds to one of the flowFrames in the \code{frames} slot.
It is mandatory that the pData has column named \code{name}}
\item{\code{colnames}:}{A \code{character} object with the
(common) column names of all the data matrices in the flowFrames.}
}}

\section{Methods}{\describe{
\item{[, [[}{subsetting. If \code{x} is \code{flowSet}, then \code{x[i]}
  returns a \code{flowSet} object, and \code{x[[i]]} a \code{flowFrame}
  object. The semantics is similar to the behavior of the subsetting operators for lists.}
\item{fsApply}{apply fucntion on a \code{flowSet} object.}
\item{compensate}{apply a compensation matrix on a \code{flowSet} object.}
\item{transform}{apply a transformation function on a \code{flowSet}
  object.}
\item{colnames, colnames<-}{extract or replace the \code{colnames} slot.}
\item{phenoData, phenoData<-}{extract or replace the \code{phenoData}
  slot.}
\item{pData}{extract actual information from the \code{phenoData} slot.}
\item{sampleNames}{extract sample names from the \code{phenoData} object.}
\item{length}{number of \code{flowFrame} objects in the set.}
\item{keyword}{extract entries from the description slots of the
  \code{flowFrame} objects by keywords}
\item{show}{display summary.}
\item{filter}{apply a filter object on a \code{flowSet}
  object. See \code{\link[flowCore]{filter} for details}}
\item{spillover}{compute spillover matrix from a compensation set.}
\item{split}{split all \code{flowFrame} objects according to a filter or
  gate. Returns a list of \code{flowFrames}}
\item{rbind2}{combine two \code{flowSet} objects}
\item{Subset}{add description here}
}}

\section{Important note on storage and performance}{
The bulk of the data in a \code{flowSet} object is stored in an
\code{\link[base:environment]{environment}}, and is therefore not
automatically copied when the \code{flowSet} object is copied. If
\code{x} is an object of class \code{flowSet}, then the code
\preformatted{y <- x} will create an object \code{y} that contains
copies of the \code{\link[Biobase:phenoData-class]{phenoData}} and
administrative data in \code{x}, but refers to the \emph{same}
environment with the actual fluorescence data. See below for how to
create proper copies.

The reason for this is performance. The pass-by-value semantics of
function calls in \code{R} can result in numerous copies of the same
data object being made in the course of a series of nested function
calls. If the data object is large, this can result in a considerable
cost of memory and performance. \code{flowSet} objects are intended to
contain experimental data in the order of hundreds of Megabytes,
which can effectively be treated as read-only: typical tasks
are the extraction of subsets and the calculation of summary statistics.
This is afforded by the design of the \code{flowSet}
class: an object of that class contains a \code{phenoData} slot,
some administrative information, and a \emph{reference} to an
environment with the fluorescence data; when it is copied, only the
reference is copied, but not the potentially large set of fluorescence
data themselves.

However, note that subsetting operations, such as
\preformatted{y <- x[i]}
do create proper copies, including a copy of the appropriate part of the
fluorescence data, as it should be expected. Thus, to make a proper copy
of a \code{flowSet} \code{x}, use
\preformatted{y <- x[seq(along=x)]}
}

\author{F. Hahne, W. Huber}

\seealso{
\code{\link[flowCore]{flowFrame}}
}

\examples{

}

}
\keyword{classes}

