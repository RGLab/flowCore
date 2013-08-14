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
  y <- .Call("biexponential_transform", x, a, b, c, d, f, w, tol, maxit)

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

