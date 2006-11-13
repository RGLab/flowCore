\name{fcsFrame-class}
\docType{class}
\alias{fcsFrame-class}
\alias{fcsFrame}
\alias{[,fcsFrame-method}
\alias{description,fcsFrame-method}
\alias{description<-,fcsFrame,character-method}
\alias{exprs,fcsFrame-method}
\alias{exprs<-,fcsFrame,matrix-method}
\alias{colnames,fcsFrame-method}
\alias{colnames<-,fcsFrame-method}
\alias{show,fcsFrame-method}
\alias{plot,fcsFrame,missing-method}
\alias{$.fcsFrame}
\alias{gate,fcsFrame-method}
\alias{gate<-,fcsFrame,gate-method}
\alias{gate<-,fcsFrame,gateSet-method}
\alias{ncol,fcsFrame-method}
\alias{nrow,fcsFrame-method}



\title{'fcsFrame': a class for storing observed quantitative
  properties from a population of cells, most likely from a FACS run or, 
  alternatively, from automated microscopy}

\description{This class represents the data contained in a
  FCS 3.0 file or similar data structures.}


\section{Creating Objects}{
Objects can be created using\cr
\code{  new('fcsFrame,}\cr
\code{    exprs       = ...., # Object of class matrix}\cr
\code{    description = ....  # Object of class character}\cr
\code{  )}\cr
or the function \code{\link[prada]{readFCS}}.}

\section{Slots}{\describe{
\item{\code{exprs}:}{Object of class \code{matrix} containing
      the measured intensities. Rows correspond to cells, columns to the
      different channels. The \code{colnames} attribute of the matrix is
      supposed to hold the names or identifiers for the channels. The
      \code{rownames} attribute would usually not be set.}
\item{\code{description}:}{A named character vector
  containing the experiment description as key-value pairs.}
}}

\section{Methods}{\describe{
\item{[}{subsetting. Returns an object of class \code{fcsFrame}.
  The subsetting is applied to the \code{exprs} slot, while the
  \code{description} slot is unchanged.}
\item{exprs, exprs<-}{extract or replace the intensities.}
\item{description, description<-}{extract or replace the description.}
\item{show}{display summary.}
\item{plot}{scatterplot for \code{fcsFrame} objects. The additional
argument \code{gate} can be used to plot subsets of the data defined by
either an object of class \code{\link[prada]{gate}} or by a character
vector giving the name of one of the gates in the list.}
\item{ncol,nrow}{extract the dimensions of the data matrix.}
}}

\details{
Although objects of class \code{fcsFrame} can be used to hold arbitrary
data of cell populations, the main focus lies on flow-cytometry data.
 
FCS 3.0 is the Data File Standard for Flow Cytometry, Version
FCS3.0. See the vignette of this package for additional information on
using the object system for handling of flow-cytometry data. 
}
			
\author{Florian Hahne, Wolfgang Huber}

\seealso{}

\examples{
\dontrun{
intens <- matrix(runif(100), ncol=4)
colnames(intens) <- c("FL1-H", "FL2-H", "FL3-H", "FL4-H")

a <- new("fcsFrame",
          exprs=intens,
          description=c(name="example data", date=date()))

description(a)
dim(exprs(a))           

a[1:3, -4]

plot(a)

}
}

\keyword{classes}

