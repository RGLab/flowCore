source("blendFuns.R")

## A function that blends two color palettes of varying intensities.
## The first argument is one of the blending functions ('bl_xyz').
## Additional arguments can be passed to them (e.g. alpha)
blendCols <- function(blendFun=bl_addition, range1=c("blue", "red"),
                      range2=c("green", "yellow"), n=20, ...){

  ## create ramps for range1 and range2
  rr1 <- colorRampPalette(range1)(n)
  rr2 <- colorRampPalette(range2)(n)
  
  ## create color ramps from each color to white
  colMat1 <- t(sapply(rr1, function(x) colorRampPalette(c(x, "white"))(n)))
  colMat2 <- t(sapply(rr2, function(x) colorRampPalette(c("white", x))(n)))

  ## apply a belding function to the color matrices
  addColMat <- blendWrapper(colMat1, colMat2, blendFun, ...)

  ## plot background, forground and merged colors
  par(mfrow=c(1,3))
  plot(0, xlim=c(0, n), ylim=c(0, n), type="n", ann=FALSE, axes=FALSE)
  mtext("fg")
  
  rect(rep(0:((n)-1), each=n), 0:((n)-1), rep(1:(n), each=n),
       1:(n), col=colMat1[(n):1,], lty=0)
  box()
  plot(0, xlim=c(0, n), ylim=c(0, n), type="n", ann=FALSE, axes=FALSE)
  mtext("bg")
  rect(rep(0:((n)-1), each=n), 0:((n)-1), rep(1:(n), each=n),
       1:(n), col=colMat2[(n):1,], lty=0)
  box()
  plot(0, xlim=c(0, n), ylim=c(0, n), type="n", ann=FALSE, axes=FALSE)
  mtext("blend")
  rect(rep(0:((n)-1), each=n), 0:((n)-1), rep(1:(n), each=n),
       1:(n), col=addColMat[(n):1,], lty=0)
  box()
}

## A helper function for smoothScatterMult to compute densities
smoothScatterMultCalcDensity <- function (x, nbin, bandwidth, range)
{
    require(KernSmooth)
    if (length(nbin) == 1)
        nbin <- c(nbin, nbin)
    if (!is.numeric(nbin) || (length(nbin) != 2))
        stop("'nbin' must be numeric of length 1 or 2")
    if (missing(bandwidth) | is.null(bandwidth)) {
        bandwidth <- diff(apply(x, 2, quantile, probs = c(0.05,
            0.95), na.rm = TRUE))/25
    }
    else {
        if (!is.numeric(bandwidth))
            stop("'bandwidth' must be numeric")
    }
    rv <- bkde2D(x, gridsize = nbin, bandwidth = bandwidth, range.x=range)
    rv$bandwidth <- bandwidth
    return(rv)
}

## A helper function to sequentially blend all columns of a color matrix
mergeCols <- function(colors, mFun, dims, bfArgs){
 if(!is(colors, "matrix"))
   stop("'colors' must be matrix containing color values")
 d <- dim(colors)
 mCols <- do.call(blendWrapper, args=c(list(background=colors[,1],
                                    foreground=colors[,2], blendFun=mFun),
                                    bfArgs))
 i <- 3
 while(i<=d[2]){
   mCols <- do.call(blendWrapper, args=c(list(background=mCols,
                                    foreground=colors[,i], blendFun=mFun),
                                    bfArgs))
   i <- i+1
 }
 dim(mCols) <- dims
 return(mCols)
}

## A helper function to compute densities and corresponding colors
computeDensCols <- function(z, ramp, nbin, transformation,
                            bandwidth, range){
  z <- matrix(z, ncol=2)
  map <- smoothScatterMultCalcDensity(z, nbin, bandwidth=bandwidth,
                                      range=range)
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim = dim(dens))
  colpal <- as.integer(1 + (length(dens) - 1) * dens/max(dens))
  cols <- rep(as.character(NA), prod(dim(dens)))
  cols <- matrix(ramp(length(dens))[colpal], ncol=nbin)
  return(list(density=dens, color=cols, xm=xm, ym=ym))
}

## An addition to the smoothScatter function that takes as additional
## argument a factor according to which the rows of the data matrix
## are classified. FIXME: Make smoothScatter a generic and add mthods.
smoothScatterMult <- function(x, y = NULL, fact, nbin = 128, bandwidth=NULL, 
                              nrpoints = 100, blendFun=bl_alpha, bfArgs=list(),
                              transformation = function(x) x^0.25,
                              xlab = NULL, ylab = NULL, postPlotHook = box,
                              pch = ".", cex = 1, ...){

  ## argument validation
  require(RColorBrewer)
  if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) != 1))
    stop("'nrpoints' should be numeric scalar with value >= 0.")

  ## create xy.coords object
  xlabel <- if (!missing(x))
    deparse(substitute(x))
  ylabel <- if (!missing(y))
    deparse(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab))
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab))
    xy$ylab
  else ylab
  odims <- if(length(nbin)==2)
    nbin
  else rep(nbin,2)

  ## do factorization
  sel <- !(is.na(xy$x) | is.na(xy$y))
  if(missing(fact))
    fact <- factor(rep(1, length(xy$x)))
  else
    fact <- as.factor(fact)
  fact <- as.factor(fact)
  if(length(fact) != length(xy$x))
    stop("'fact' must be factor vector mapping to the data matrix")
  range <- lapply(xy[1:2], range, na.rm=TRUE)
  nf <- length(levels(fact))

  ## the default colramps (FIXME: Needs to be generalized for n>3 ramps)
  colramps <- c(colorRampPalette(c("white", brewer.pal(9, "Reds"))),
                colorRampPalette(c("white", brewer.pal(9, "Blues"))),
                colorRampPalette(c("white", brewer.pal(9, "Greens"))))
  
  ## split data matrix according to factors
  tmp <- cbind(xy$x, xy$y)[sel,]
  xySplit <- split(tmp, fact[sel])
  

  ## compute densities and colors for each item 
  densCols <- mapply(computeDensCols, xySplit, ramp=colramps[1:nf],
                     MoreArgs=list(nbin=nbin, transformation=transformation,
                       bandwidth=bandwidth, range=range), SIMPLIFY=FALSE)
  densities <- sapply(densCols, function(x) x$density)
  colors <- sapply(densCols, function(x) x$color)
  xm <- densCols[[1]]$xm
  ym <- densCols[[1]]$ym

  ## merge the overlapping colors
  if(ncol(colors)>1)
   mCols <- mergeCols(colors, mFun=blendFun, dims=odims, bfArgs=bfArgs)
  else
    mCols <- matrix(colors, nrow=odims[1], ncol=odims[2])

  ## create plot
  xdelta <- diff(xm)
  ydelta <- diff(ym)
  xm1 <- xm[-length(xm)]
  ym1 <- ym[-length(ym)]
  mCols <- mCols[-odims[1], -odims[2]]
  plot(xm,ym, type="n", xlab=xlab, ylab=ylab, ...)
  rect(rep(xm1, each=odims[1]-1), rep(ym1, odims[2]-1),
       rep(xm1+xdelta, each=odims[1]-1), rep(ym1+ydelta, odims[2]-1),
       col=t(mCols), lty=0)

  ## post plot hook
  if (!is.null(postPlotHook))
    postPlotHook()

  ## add points to plot
  if (nrpoints != 0) {
    dens <- rowSums(densities)
    dim(dens) <- odims
    stopifnot(length(xm) == nrow(dens), length(ym) == ncol(dens))
    ixm <- round((xy$x - xm[1])/(xm[length(xm)] - xm[1]) *
                 (length(xm) - 1))
    iym <- round((xy$y - ym[1])/(ym[length(ym)] - ym[1]) *
                 (length(ym) - 1))
    idens <- dens[1 + iym * length(xm) + ixm]
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    sel <- order(idens, decreasing = FALSE)[1:nrpoints]
    points(xy$x[sel], xy$y[sel], pch = pch, cex = cex, col = "black")
  }
}

