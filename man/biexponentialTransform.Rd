\name{biexponentialTransform}
\alias{biexponentialTransform}

\title{Compute a transform using the 'biexponential' function}
\description{
  The 'biexponential' is an over-parameterized inverse of the hyperbolic sine.
  The function to be inverted takes the form 
	biexp(x) = a*exp(b*(x-w))-c*exp(-d*(x-w))+f with default parameters 
	selected to correspond to the hyperbolic sine.
}
\usage{
biexponentialTransform(transformationId="defaultBiexponentialTransform", a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0, tol = .Machine$double.eps^0.25, maxit = as.integer(5000))
}

\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used by the transform/filter integration routines. }
  \item{a}{ See the function description above. Defaults to 0.5 }
  \item{b}{ See the function description above. Defaults to 1.0 }
  \item{c}{ See the function description above. Defaults to 0.5 (the same as \code{a}) }
  \item{d}{ See the function description above. Defaults to 1 (the same as \code{b}) }
  \item{f}{ A constant bias for the intercept. Defaults to 0. }
  \item{w}{ A constant bias for the 0 point of the data. Defaults to 0. }
  \item{tol}{ A tolerance to pass to the inversion routine (\code{\link{uniroot}} usually) }
  \item{maxit}{ A maximum number of iterations to use, also passed to \code{\link{uniroot}} }
}

\value{
	Returns values giving the inverse of the biexponential within a certain tolerance. This function should be used with care as numerical inversion routines often have problems with the inversion process due to the large range of values that are essentially 0. Do not be surprised if you end up with population splitting about \code{w} and other odd artifacts.
}

\author{ B. Ellis, N Gopalakrishnan }
\seealso{ \code{\link{transform}} }
\examples{
# Construct some "flow-like" data which tends to be hetereoscedastic.
data(GvHD)
biexp  <- biexponentialTransform("myTransform")
after.1 <- transform(GvHD, `FSC-H`= biexp(`FSC-H`))

biexp  <- biexponentialTransform("myTransform",w=10)
after.2 <- transform(GvHD, `FSC-H`= biexp(`FSC-H`))

opar = par(mfcol=c(3, 1))
plot(density(exprs(GvHD[[1]])[, 1]), main="Original")
plot(density(exprs(after.1[[1]])[, 1]), main="Standard Transform")
plot(density(exprs(after.2[[1]])[, 1]), main="Shifted Zero Point")
}

\keyword{methods}

