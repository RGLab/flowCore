\name{Subset}
\alias{Subset}
\alias{Subset,flowFrame,filter-method}
\alias{Subset,flowSet,ANY-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subset a flowFrame or a flowSet}
\description{
An equivalent of a \code{\link{subset}} function for \code{\link{flowFrame}} or
a \code{\link{flowSet}} object.
}
\usage{
Subset(x, subset, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The flow object, frame or set, to subset.}
  \item{subset}{A filter object or, in the case of \code{flowSet} subsetting, a named list of filters.}
  \item{\dots}{Like the original \code{\link{subset}} function, you can also select columns.}
}
\details{
The \code{Subset} method is the recommended method for obtaining a \code{flowFrame} that only contains events
consistent with a particular filter. It is functionally equivalent to \code{frame[as(filter(frame,subset),"logical"),]} 
when used in the \code{flowFrame} context. Used in the \code{flowSet} context, it is equivalent to using \code{fsApply}
to apply the filtering operation to each \code{flowFrame}. 

Additionally, using \code{Subset} on a \code{flowSet} can also take a named \code{list} as the subset. In this case, 
the names of the list object should correspond to the \code{sampleNames} of the flowSet, allowing a different filter 
to be applied to each frame. If not all of the names are used or excess names are present a warning will be generated
but the valid filters will be applied for the rare instances where this is the intended operation. Note that a 
\code{\link{filter}} operation will generate a list of \code{\linkS4class{filterResult}} objects that can be used directly with
\code{Subset} in this manner.
}
\value{
  Depending on the original context, either a \code{\link{flowFrame}} or
a \code{\link{flowSet}}.
}
\author{B. Ellis}
\seealso{\code{\link{split}},\code{\link{subset}}}
\examples{
sample = read.flowSet(path=system.file("extData",package="flowCore"),patt="0877408774")
result = filter(sample,rectangleGate("FSC-H"=c(-Inf,1024)))
result
Subset(sample,result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

