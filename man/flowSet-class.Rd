\name{flowSet-class}
\docType{class}
\alias{flowSet}
\alias{flowSet-class}
\alias{[,flowSet-method}
\alias{[[,flowSet-method}
\alias{[[<-,flowSet-method}
\alias{show,flowSet-method}
\alias{length,flowSet-method}
\alias{pData,flowSet-method}
\alias{colnames,flowSet-method}
\alias{colnames<-,flowSet-method}
\alias{phenoData,flowSet-method}
\alias{phenoData<-,flowSet,phenoData-method}
\alias{plot,flowSet,missing-method}
\alias{hist,flowSet-method}


\title{'flowSet': a class for storing raw data from
a quantitative cell-based assay}

\description{This class is a container for a set of
  \code{\link[flowCore]{flowFrame}} objects}


\section{Creating Objects}{
Objects can be created using the function
\code{\link[flowCore]{readflowSet}} or via\cr
\code{  new('flowSet,}\cr
\code{    frames = ....,   # environment with flowFrames}\cr
\code{    phenoData = .... # object of class phenoData}\cr
\code{    colnames = ....  # object of class character}\cr
\code{  )}\cr
}

\section{Slots}{\describe{
\item{\code{frames}:}{An \code{\link[base:environment]{environment}}
containing one or more \code{\link[flowFrame-class]{flowFrame}} objects.}
\item{\code{phenoData}:}{A
\code{\link[Biobase:phenoData-class]{phenoData}}. Each row
corresponds to one of the flowFrames in the \code{frames} slot.
It is mandatory that the pData has column named \code{name}}
\item{\code{colnames}:}{A \code{character} object with the
(common) column names of all the data matrices in the flowFrames.}
}}

\section{Methods}{\describe{
\item{[, [[}{subsetting. If \code{x} is \code{flowSet}, then \code{x[i]}
returns a \code{flowSet} object, and \code{x[[i]]} a \code{flowFrame}
object. The semantics is similar to the behavior of the subsetting
operators for lists.}
\item{colnames, colnames<-}{extract or replace the \code{colnames} slot.}
\item{phenoData, phenoData<-}{extract or replace the \code{phenoData} slot.}
\item{show}{display summary.}
\item{plot}{Scatterplot of one or all (consecutively) \code{\link[flowFrame]{flowFrame}}
objects.  The additional argument \code{gate} can be used to plot
subsets of the data defined by an object of class
\code{\link[flowCore]{filter}} or \code{\link[flowCore]{filterSet}}.}
\item{hist}{Draw histogram of the data. The additional argument
\code{variable} can be used to subset to one variable prior to plotting.}
}}

\section{Important note on storage and performance}{
The bulk of the data in a \code{flowSet} object is stored in an
\code{\link[base:environment]{environment}}, and is therefore not
automatically copied when the \code{flowSet} object is copied. If
\code{x} is an object of class \code{flowSet}, then the code
\preformatted{y <- x} will create a an object \code{y} that contains
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

\author{Florian Hahne}

\seealso{
\code{\link[flowCore]{flowFrame}}
}

\examples{

}

}
\keyword{classes}

