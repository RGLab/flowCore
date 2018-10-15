## FCSTrans Authors
## Authors: Yue Liu and Yu "Max" Qian
# Contact: yu.qian@utsouthwestern.edu or yliu0@yahoo.com
 
# set output to 0 when input is less than cutoff value
ipfloor <- function (x, cutoff = 0, target = 0) {
  ifelse(x <= cutoff, target, x)
}

# set output to 0 when input is less than cutoff value
ipceil <- function (x, cutoff = 0, target = 0) {
  ifelse(x >= cutoff, target, x)
}

# calculation core of iplogicle
iplogicore <- function (x, w, r, d, scale, rescale = TRUE, maxit = 5000,
                        tol = sqrt(.Machine$double.eps)) {
  
  x <- as.numeric(x)
  d <- d * log(10)
  scale <- ifelse(rescale, scale / d, 1)

  if (w == 0) {
    p <- 1
  } else {
    p <- uniroot(f = function(p) { -w + 2 * p * log(p) / (p + 1)},
            interval = c(.Machine$double.eps, 2 * (w + d)))$root
  }
  a <- r * exp(-(d - w))
  b <- 1
  c <- r * exp(-(d - w)) * p^2
  d <- 1 / p
  f <- a * (p^2 - 1)
  y <- biexponential_transform(x, a, b, c, d, f, w, tol, maxit)

  sapply(y * scale, ipfloor)
}

# function for calculating w 
iplogiclew <- function (w, cutoff = -111, r = 262144, d = 4.5, scale = 1) {
  if (w > d) {
    w <- d
  }
  iplogicore(cutoff, w, r, d, scale) - sqrt(.Machine$double.eps)
}

# import logicle function - convert fluorescent marker values to channel output
# rescale parameter will transform the data to cover the range, otherwise it will
# be on the scale of the actual transformation.
iplogicle <- function (x, r = 2^18, d = 4.5, range = 4096, cutoff = -111,
                       w = NULL, rescale = TRUE) {
  if (is.null(w)) {
    w <- uniroot(iplogiclew, interval = c(0, d), cutoff = cutoff)$root
  } else if (w > d) {
    stop("Negative range decades must be smaller than total number of decades")
  }

  print(paste("params: r=", r, "d=", d, "range=", range, "cutoff=", cutoff, "w=", w))
  iplogicore(x, w, r, d, range, rescale = rescale)
}

# FlowCore Transformation using FCSTrans
#' Computes a transform using the 'iplogicle' function
#' 
#' Transforms FCS data using the iplogicle function from FCSTrans by Quian et
#' al.  The core functionality of FCSTrans has been imported to produce
#' transformed FCS data rescaled and truncated as produced by FCSTrans. The
#' \code{w} parameter is estimated by \code{iplogicle} automatically, then
#' makes a call to \code{iplogicore} which in turn uses the logicle transform
#' code of Wayne Moore.
#' 
#' For the details of the FCSTrans transformation, we recommend the excellent
#' Supplementary File that accompanies Quian et al. (2012):
#' \url{http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22037/suppinfo}
#' 
#' @usage 
#' FCSTransTransform(transformationId = "defaultFCSTransTransform", 
#'                   channelrange = 2^18, channeldecade = 4.5, 
#'                   range = 4096, cutoff = -111, w = NULL, rescale = TRUE)
#' 
#' @param transformationId A name to assign to the transformation. Used by the
#' transform/filter routines.
#' @param channelrange is the range of the data. By default, 2^18 = 262144.
#' @param channeldecade is the number of logarithmic decades. By default, it is
#' set to 4.5.
#' @param range the target resolution. The default value is 2^12 = 4096.
#' @param cutoff a threshold below which the logicle transformation maps values
#' to 0.
#' @param w the logicle width. This is estimated by \code{iplogicle} by
#' default. Details can be found in the Supplementary File from Quian et al.
#' @param rescale logical parameter whether or not the data should be rescaled
#' to the number of channels specified in \code{range}. By default, the value
#' is \code{TRUE} but can be set to FALSE if you want to work on the
#' transformed scale.
#' 
#' 
#' @author Wayne Moore, N Gopalakrishnan
#' @seealso 
#' \code{\link[flowCore]{inverseLogicleTransform}},
#' \code{\link[flowCore]{estimateLogicle} },
#' \code{\link[flowCore]{logicleTransform}}
#' @references Y Quian, Y Liu, J Campbell, E Thompson, YM Kong, RH Scheuermann;
#' FCSTrans: An open source software system for FCS file conversion and data
#' transformation. Cytometry A, 2012
#' @keywords methods
#' @examples
#' 
#' data(GvHD)
#' samp <- GvHD[[1]] 
#' ## User defined logicle function
#' lgcl <- transformList(c('FL1-H', 'FL2-H'), FCSTransTransform())
#' after <- transform(samp, lgcl)
#' 
#' @export
FCSTransTransform <- function(transformationId = "defaultFCSTransTransform",
                              channelrange = 2^18,
                              channeldecade = 4.5,
                              range = 4096, cutoff = -111, w = NULL,
                              rescale = TRUE) {
  k <- new("transform", .Data = function(x) {
    x <- iplogicle(x, r = channelrange, d = channeldecade, range = range,
                   cutoff = cutoff, w = w, rescale = rescale)
  })
  k@transformationId <- transformationId;
  k
}

