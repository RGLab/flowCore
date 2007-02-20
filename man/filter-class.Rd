\name{filter-class}
\docType{class}
\alias{filter-class}
\alias{filtergate,filter-class}
\alias{rectangleGate,filter-class}
\alias{polygonGate,filter-class}
\alias{polytopeGate,}
\alias{ellipsoidGate,filter-class}
\alias{norm2Filter,filter-class}
\alias{decisionTreeGate,filter-class}
\alias{booleanGate,filter-class}
\alias{filter, filter-method}
\title{A class for representing filtering operations to be applied to flow data.}
\description{
The \code{filter} class is the base class for all filter/gating objects in \code{flowCore}. In general you will want to
subclass or create a more specific filter.
}

\section{Slots}{\describe{
\item{\code{filterId}:}{A character vector that identifies this \code{filter}. This is typically user specified but can be
automatically deduced by certain filter operations, particularly boolean and set operations.}
\item{\code{parameters}:}{A character vector the identifies the name of the parameters that this filter operators upon.}
}}

\section{Methods}{\describe{
\item{\code{\%in\%}}{Used in the usual way this returns a vector of values that identify which events were accepted by the filter. A single filter may encode several populations so this can return either a \code{logical} vector, a \code{factor} vector or a \code{numeric} vector of probabilities that the event is accepted by the filter. Minimally, you must implement this method when
creating a new filter}
\item{\code{&}, \code{|}, \code{!}}{Two filters can be composed using the usual boolean operations returning a \code{filter} class
of a type appropriate for handling the operation. These methods attempt to guess an appropriate \code{filterId} for the new \code{filter}}
\item{\code{\%subset\%}, \code{\%&\%}}{Defines a filter as being a subset of another filter. For deterministic filters the results
will typically be equivalent to using an \code{&} operation to compose the two filters, though summary methods will use subset 
semantics when calculating proportions. Additionally, when the filter is data driven, such as \filter{norm2Filter}, the subset
semantics are applied to the data used to fit the filter possibly resulting in quite different, and usually more desirable, results.}
}
\item{\code{\%on\%}}{Used in conjunction with a \code{\link{transformList}} to create a \code{transformFilter}. This filter is
similar to the subset filter in that the filtering operation takes place on transformed values rather than the original values.}
\item{\code{filter}}{A more formal version of \code{\%in\%}, this method returns a \code{\link{filterResult}} object that can be
used in subsequent filter operations as well as providing more metadata about the results of the filtering operation}}
\item{\code{summarizeFilter}}{When implementing a new filter this method is used to update the \code{filterDetails} slot of a
 \code{filterResult}. It is optional and typically only needs to be implemented for data-driven filters.}


\details{}

\seealso{\code{\link{transform}},\code{\link{filter}}}

\author{B. Ellis, P.D. Haaland and N. Le Meur}

\keyword{classes}