\name{compensation-class}

\docType{class}

\alias{compensation-class}
\alias{compensation}
\alias{compensate,flowFrame,compensation-method}
\alias{identifier,compensation-method}
\alias{identifier<-,compensation,character-method}
\alias{show,compensation-method}
\alias{compensate}
\alias{compensate,flowFrame-method}
\alias{compensate,flowSet-method}
\alias{compensate,matrix-method}

\title{Class "compensation"}

\description{

  Class and methods to compensate for spillover between channels by
  applying a spillover matrix to a \code{flowSet} or a \code{flowFrame}
  assuming a simple linear combination of values.

}


\section{Objects from the Class}{
  
  Objects should be created using the contructor
  \code{compensation()}. See the \code{Usage} and \code{Arguments}
  sections for details.
  
}

\section{Slots}{
  
  \describe{
    
    \item{\code{spillover}:}{Object of class \code{matrix}; the
      spillover matrix. Note that you can switch from a spillover to a
      compensation matrix representation using the \code{inv} argument
      of either the constructor or the \code{compensate} methods. }
    
    \item{\code{compensationId}:}{Object of class \code{character}. An
      identifier for the object. }

     \item{\code{parameters}:}{Object of class \code{parameters}. The
       flow parameters for which the compensation is defined. This can
       also be objects of class \code{\linkS4class{transform}}, in which
       case the compensation is peformed on the compensated parameters.}
    
  }
}


\usage{

compensation(\dots, spillover, inv=FALSE,
compensationId="defaultCompensation")

compensate(x, spillover, \dots)

}

\arguments{

  \item{spillover}{ The spillover or compensation matrix. }

  \item{inv}{ A logical to swich between the two matrix
    representation (spillover or compensation matrix). }

  \item{compensationId}{ The identifier for the compensation object. }

  \item{\dots}{Further arguments. For the \code{compensate} methods,
    this is only \code{inv}, a logical that forces the function to use
    the inverse of the spillover matrix.

    The constructor is designed to be useful in both programmatic and
    interactive settings, and \dots serves as a container for possible
    arguments. The following combinations of values are allowed:

    Elements in \dots are \code{character} scalars of parameter names or
    \code{\linkS4class{transform}} objects and the colnames in
    \code{spillover} match to these parameter names.
  
    The first element in \dots is a \code{character} vector of parameter
    names or a list of \code{character} scalars or
    \code{\linkS4class{transform}} objects and the colnames in
    \code{spillover} match to these parameter names.

    Argument \code{spillover} is missing and the first element in \dots
    is a \code{matrix}, in which case it is assumed to be the spillover
    matrix.

    \dots is missing, in which case all parameter names are taken from the
    colnames of \code{spillover}.
  
  }
}

\value{

  A \code{compensation} object for the constructor.

  A \code{\linkS4class{flowFrame}} or \code{\linkS4class{flowSet}} for
  the \code{compensate} methods.
  
}


\section{Methods}{
  
  \describe{
    
    \item{compensate}{\code{signature(x = "flowFrame", spillover =
	"compensation")}: Apply the compensation defined in a
      \code{compensation} object on a \code{\linkS4class{flowFrame}}.
      This returns a compensated \code{flowFrame}.
      
      \emph{Usage:}
      
      \code{   compensate(flowFrame, compensation)}

    }

    \item{compensate}{\code{signature(x = "flowFrame", spillover =
	"matrix")}: Apply a compensation matrix to a
	\code{\linkS4class{flowFrame}}.  This returns a compensated
	\code{flowFrame}.
      
      \emph{Usage:}

      \code{   compensate(flowFrame, matrix)}

    }

    \item{compensate}{\code{signature(x = "flowFrame", spillover =
	"data.frame")}:Try to coerce the \code{data.frame} to a
	\code{matrix} and apply that to a
	\code{\linkS4class{flowFrame}}.  This returns a compensated
	\code{flowFrame}.
      
      \emph{Usage:}

      \code{   compensate(flowFrame, data.frame)}

    }
    
    \item{identifier, identifier<-}{\code{signature(object =
      "compensation")}: Accessor and replacement methods for the
    \code{compensationId} slot.

    \emph{Usage:}
      
      
      \code{   identifier(compensation)}

      \code{   identifier(compensation) <- value}
      
    }
    
    \item{show}{\code{signature(object = "compensation")}: Print details
      about the object.
      
      \emph{Usage:}

      This method is automatically called when the object is printed on
      the screen.
      
    }
    
  }
}

\details{

  The essential premise of compensation is that some fluorochromes may
  register signals in detectors that do not correspond to their primary
  detector (usually a photomultiplier tube). To compensate for this
  fact, some sort of standard is used to obtain the background signal
  (no dye) and the amount of signal on secondary channels for each
  fluorochrome relative to the signal on their primary channel.

  To calculate the spillover percentage we use either the mean or the
  median (more often the latter) of the secondary signal minus the
  background signal for each dye to obtain \code{n} by \code{n} matrix,
  \code{S}, of so-called spillover values, expressed as a percentage of
  the primary channel. The observed values are then considered to be a
  linear combination of the true fluorescence and the spillover from
  each other channel so we can obtain the true values by simply
  multiplying by the inverse of the spillover matrix.

  The spillover matrix can be obtained through several means. Some flow
  cytometers provide a spillover matrix calculated during acquisition,
  possibly by the operator, that is made available in the metadata of
  the flowFrame.  While there is a theoretical standard keyword
  \code{\$SPILL} it can also be found in the \code{SPILLOVER} or
  \code{SPILL} keyword depending on the cytometry. More commonly the
  spillover matrix is calculated using a series of compensation cells or
  beads collected before the experiment. If you have set of FCS files
  with one file per fluorochrome as well as an unstained FCS file you
  can use the \code{\link[flowCore]{spillover}} method for
  \code{\link[flowSet-class]{flowSets}} to automatically calculate a
  spillover matrix.

  The \code{compensation} class is essentially a wrapper around a
  \code{matrix} that allows for transformed parameters and method
  dispatch.
 
}


\seealso{
  
  \code{\link[flowCore]{spillover}}
  
}

\author{F.Hahne, B. Ellis, N. Le Meur}

\examples{
## Read sample data and a sample spillover matrix
samp   <- read.flowSet(path=system.file("extdata", "compdata", "data", package="flowCore"))
cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
comp.mat

## compensate using the spillover matrix directly
summary(samp)
samp <- compensate(samp, comp.mat)
summary(samp)

## create a compensation object and compensate using that
comp <- compensation(comp.mat)
compensate(samp, comp)

}

\keyword{methods}  
\keyword{classes}

