## These functions copied from CRAN package feature_1.2.10
## See: http://cran.fhcrc.org/web/packages/feature/index.html

SignifFeatureData <- function (x, d, dest, SignifFeature)
{
  n <- nrow(x)
  x.ind <- matrix(0, ncol = d, nrow = n)
  for (j in 1:d) x.ind[, j] <- findInterval(x[, j], dest$x.grid[[j]])
  return(SignifFeature[x.ind])
}

SignifFeatureRegion <- function (n, d, gcounts, gridsize, dest, bandwidth, signifLevel,
  range.x, grad = TRUE, curv = TRUE, neg.curv.only = TRUE)
{
  h <- bandwidth
  ESS <- n * dest$est * prod(h) * (sqrt(2 * pi)^d)
  SigESS <- ESS >= 5
  Sig.scalar <- array(NA, dim = gridsize)
  Sig2.scalar <- array(NA, dim = gridsize)
  dest$est[dest$est < 0] <- 0
  Sig.scalar <- 1/2 * (2 * sqrt(pi))^(-d) * n^(-1) * prod(h)^(-1) *
    dest$est
  if (d == 1)
    Sig2.scalar <- (8 * sqrt(pi) * n * prod(h))^(-1) * dest$est
  else if (d == 2)
    Sig2.scalar <- (16 * pi * n * prod(h))^(-1) * dest$est
  else if (d == 3)
    Sig2.scalar <- (32 * pi^(3/2) * n * prod(h))^(-1) * dest$est
  else if (d == 4)
    Sig2.scalar <- (64 * pi^2 * n * prod(h))^(-1) * dest$est
  matrix.sqrt <- function(A) {
    sva <- svd(A)
    if (min(sva$d) >= 0)
      Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
    else stop("Matrix square root is not defined")
    return(Asqrt)
  }
  if (d > 1) {
    WaldGrad <- array(NA, dim = gridsize)
    WaldCurv <- array(NA, dim = gridsize)
    local.mode <- array(FALSE, dim = gridsize)
  }
  if (d == 1) {
    if (grad) {
      obj1 <- drvkde(gcounts, drv = 1, bandwidth = h, binned = TRUE,
        range.x = range.x, se = FALSE)
      fhat1 <- obj1$est
      Sig.inv12 <- 1/sqrt(Sig.scalar * h^(-2))
      WaldGrad <- (Sig.inv12 * fhat1)^2
    }
    if (curv) {
      obj2 <- drvkde(gcounts, drv = 2, bandwidth = h, binned = TRUE,
        range.x = range.x, se = FALSE)
      fhat2 <- obj2$est
      Sig2.inv12 <- 1/sqrt(Sig2.scalar * 3 * h^(-4))
      lambda1 <- Sig2.inv12 * fhat2
      WaldCurv <- lambda1^2
      local.mode <- (lambda1 < 0)
    }
  }
  if (d == 2) {
    if (grad) {
      obj10 <- drvkde(gcounts, drv = c(1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj01 <- drvkde(gcounts, drv = c(0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat10 <- obj10$est
      fhat01 <- obj01$est
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) if (SigESS[i1,
        i2]) {
        Sig.inv12 <- 1/sqrt(Sig.scalar[i1, i2] * h^(-2))
        WaldGrad[i1, i2] <- sum((Sig.inv12 * c(fhat10[i1,
          i2], fhat01[i1, i2]))^2)
      }
    }
    if (curv) {
      Sig2.mat <- matrix(c(3/h[1]^4, 0, 1/(h[1]^2 * h[2]^2),
        0, 1/(h[1]^2 * h[2]^2), 0, 1/(h[1]^2 * h[2]^2),
        0, 3/h[2]^4), nrow = 3, ncol = 3)
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))
      obj20 <- drvkde(gcounts, drv = c(2, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj11 <- drvkde(gcounts, drv = c(1, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj02 <- drvkde(gcounts, drv = c(0, 2), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat20 <- obj20$est
      fhat11 <- obj11$est
      fhat02 <- obj02$est
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) if (SigESS[i1,
        i2]) {
        Sig2.inv12 <- sqrt(1/Sig2.scalar[i1, i2]) * matrix.sqrt(Sig2.mat.inv)
        fhat.temp <- Sig2.inv12 %*% c(fhat20[i1, i2],
          fhat11[i1, i2], fhat02[i1, i2])
        WaldCurv[i1, i2] <- sum(fhat.temp^2)
      }
      lambda1 <- ((fhat20 + fhat02) - sqrt((fhat20 - fhat02)^2 +
          4 * fhat11^2))/2
      lambda2 <- ((fhat20 + fhat02) + sqrt((fhat20 - fhat02)^2 +
          4 * fhat11^2))/2
      local.mode <- (lambda1 < 0) & (lambda2 < 0)
    }
  }
  if (d == 3) {
    if (grad) {
      obj100 <- drvkde(gcounts, drv = c(1, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj010 <- drvkde(gcounts, drv = c(0, 1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj001 <- drvkde(gcounts, drv = c(0, 0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat100 <- obj100$est
      fhat010 <- obj010$est
      fhat001 <- obj001$est
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) for (i3 in 1:gridsize[3]) if (SigESS[i1,
        i2, i3]) {
        Sig.inv12 <- 1/sqrt(Sig.scalar[i1, i2, i3] *
            h^(-2))
        WaldGrad[i1, i2, i3] <- sum((Sig.inv12 * c(fhat100[i1,
          i2, i3], fhat010[i1, i2, i3], fhat001[i1, i2,
            i3]))^2)
      }
    }
    if (curv) {
      obj200 <- drvkde(gcounts, drv = c(2, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj110 <- drvkde(gcounts, drv = c(1, 1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj101 <- drvkde(gcounts, drv = c(1, 0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj020 <- drvkde(gcounts, drv = c(0, 2, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj011 <- drvkde(gcounts, drv = c(0, 1, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj002 <- drvkde(gcounts, drv = c(0, 0, 2), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat200 <- obj200$est
      fhat110 <- obj110$est
      fhat101 <- obj101$est
      fhat020 <- obj020$est
      fhat011 <- obj011$est
      fhat002 <- obj002$est
      Sig2.mat <- matrix(c(3/h[1]^4, 0, 0, 1/(h[1] * h[2])^2,
        0, 1/(h[1] * h[3])^2, 0, 1/(h[1] * h[2])^2, 0,
        0, 0, 0, 0, 0, 1/(h[1] * h[3])^2, 0, 0, 0, 1/(h[1] *
            h[2])^2, 0, 0, 3/h[2]^4, 0, 1/(h[2] * h[3])^2,
        0, 0, 0, 0, 1/(h[2] * h[3])^2, 0, 1/(h[1] * h[3])^2,
        0, 0, 1/(h[2] * h[3])^2, 0, 3/h[3]^4), nrow = 6,
        ncol = 6)
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) for (i3 in 1:gridsize[3]) if (SigESS[i1,
        i2, i3]) {
        Sig2.inv12 <- sqrt(1/Sig2.scalar[i1, i2, i3]) *
          matrix.sqrt(Sig2.mat.inv)
        fhat.temp <- Sig2.inv12 %*% c(fhat200[i1, i2,
          i3], fhat110[i1, i2, i3], fhat101[i1, i2, i3],
          fhat020[i1, i2, i3], fhat011[i1, i2, i3], fhat002[i1,
            i2, i3])
        D2.mat <- matrix(c(fhat200[i1, i2, i3], fhat110[i1,
          i2, i3], fhat101[i1, i2, i3], fhat110[i1, i2,
            i3], fhat020[i1, i2, i3], fhat011[i1, i2, i3],
          fhat101[i1, i2, i3], fhat011[i1, i2, i3], fhat002[i1,
            i2, i3]), nrow = 3)
        lambda <- eigen(D2.mat, symmetric = TRUE, only.values = TRUE)$values
        WaldCurv[i1, i2, i3] <- sum(fhat.temp^2)
        local.mode[i1, i2, i3] <- all(lambda < 0)
      }
    }
  }
  if (d == 4) {
    if (grad) {
      obj1000 <- drvkde(gcounts, drv = c(1, 0, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0100 <- drvkde(gcounts, drv = c(0, 1, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0010 <- drvkde(gcounts, drv = c(0, 0, 1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0001 <- drvkde(gcounts, drv = c(0, 0, 0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat1000 <- obj1000$est
      fhat0100 <- obj0100$est
      fhat0010 <- obj0010$est
      fhat0001 <- obj0001$est
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) for (i3 in 1:gridsize[3]) for (i4 in 1:gridsize[4]) if (SigESS[i1,
        i2, i3, i4]) {
        Sig.inv12 <- 1/sqrt(Sig.scalar[i1, i2, i3, i4] *
            h^(-2))
        WaldGrad[i1, i2, i3, i4] <- sum((Sig.inv12 *
            c(fhat1000[i1, i2, i3, i4], fhat0100[i1, i2,
              i3, i4], fhat0010[i1, i2, i3, i4], fhat0001[i1,
                i2, i3, i4]))^2)
      }
    }
    if (curv) {
      obj2000 <- drvkde(gcounts, drv = c(2, 0, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj1100 <- drvkde(gcounts, drv = c(1, 1, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj1010 <- drvkde(gcounts, drv = c(1, 0, 1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj1001 <- drvkde(gcounts, drv = c(1, 0, 0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0200 <- drvkde(gcounts, drv = c(0, 2, 0, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0110 <- drvkde(gcounts, drv = c(0, 1, 1, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0101 <- drvkde(gcounts, drv = c(0, 1, 0, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0020 <- drvkde(gcounts, drv = c(0, 0, 2, 0), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0011 <- drvkde(gcounts, drv = c(0, 0, 1, 1), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      obj0002 <- drvkde(gcounts, drv = c(0, 0, 0, 2), bandwidth = h,
        binned = TRUE, range.x = range.x, se = FALSE)
      fhat2000 <- obj2000$est
      fhat1100 <- obj1100$est
      fhat1010 <- obj1010$est
      fhat1001 <- obj1001$est
      fhat0200 <- obj0200$est
      fhat0110 <- obj0110$est
      fhat0101 <- obj0101$est
      fhat0020 <- obj0020$est
      fhat0011 <- obj0011$est
      fhat0002 <- obj0002$est
      Sig2.mat <- matrix(c(3/h[1]^4, 0, 0, 0, 1/(h[1] *
          h[2])^2, 0, 0, 1/(h[1] * h[3])^2, 0, 1/(h[1] *
              h[4])^2, 0, 1/(h[1] * h[2])^2, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1/(h[1] * h[3])^2, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1/(h[1] * h[4])^2, 0, 0, 0,
        0, 0, 0, 1/(h[1] * h[2])^2, 0, 0, 0, 3/h[2]^4,
        0, 0, 1/(h[2] * h[3])^2, 0, 1/(h[2] * h[4])^2,
        0, 0, 0, 0, 0, 1/(h[2] * h[3])^2, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1/(h[2] * h[4])^2, 0, 0, 0,
        1/(h[1] * h[3])^2, 0, 0, 0, 1/(h[2] * h[3])^2,
        0, 0, 3/h[3]^4, 0, 1/(h[3] * h[4])^2, 0, 0, 0,
        0, 0, 0, 0, 0, 1/(h[3] * h[4])^2, 0, 1/(h[1] *
            h[4])^2, 0, 0, 0, 1/(h[2] * h[4])^2, 0, 0,
        1/(h[3] * h[4])^2, 0, 3/h[4]^4), nrow = 10, ncol = 10)
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))
      for (i1 in 1:gridsize[1]) for (i2 in 1:gridsize[2]) for (i3 in 1:gridsize[3]) for (i4 in 1:gridsize[4]) if (SigESS[i1,
        i2, i3, i4]) {
        Sig2.inv12 <- sqrt(1/Sig2.scalar[i1, i2, i3,
          i4]) * matrix.sqrt(Sig2.mat.inv)
        fhat.temp <- Sig2.inv12 %*% c(fhat2000[i1, i2,
          i3, i4], fhat1100[i1, i2, i3, i4], fhat1010[i1,
            i2, i3, i4], fhat1001[i1, i2, i3, i4], fhat0200[i1,
              i2, i3, i4], fhat0110[i1, i2, i3, i4], fhat0101[i1,
                i2, i3, i4], fhat0020[i1, i2, i3, i4], fhat0011[i1,
                  i2, i3, i4], fhat0002[i1, i2, i3, i4])
        D2.mat <- matrix(c(fhat2000[i1, i2, i3, i4],
          fhat1100[i1, i2, i3, i4], fhat1010[i1, i2,
            i3, i4], fhat1001[i1, i2, i3, i4], fhat1100[i1,
              i2, i3, i4], fhat0200[i1, i2, i3, i4], fhat0110[i1,
                i2, i3, i4], fhat0101[i1, i2, i3, i4], fhat1010[i1,
                  i2, i3, i4], fhat0110[i1, i2, i3, i4], fhat0020[i1,
                    i2, i3, i4], fhat0011[i1, i2, i3, i4], fhat1001[i1,
                      i2, i3, i4], fhat0101[i1, i2, i3, i4], fhat0011[i1,
                        i2, i3, i4], fhat0002[i1, i2, i3, i4]), nrow = 4)
        WaldCurv[i1, i2, i3, i4] <- sum(fhat.temp^2)
        lambda <- eigen(D2.mat, symmetric = TRUE, only.values = TRUE)$values
        local.mode[i1, i2, i3, i4] <- all(lambda < 0)
      }
    }
  }
  if (grad) {
    pval.Grad <- 1 - pchisq(WaldGrad, d)
    pval.Grad.ord <- pval.Grad[order(pval.Grad)]
    num.test <- sum(!is.na(pval.Grad.ord))
    if (num.test >= 1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) -
          num.test))
    else num.test.seq <- rep(NA, prod(gridsize))
    reject.nonzero <- ((pval.Grad.ord <= signifLevel/(num.test +
        1 - num.test.seq)) & (pval.Grad.ord > 0))
    reject.nonzero.ind <- which(reject.nonzero)
    SignifGrad <- array(FALSE, dim = gridsize)
    SignifGrad[which(pval.Grad == 0, arr.ind = TRUE)] <- TRUE
    for (i in reject.nonzero.ind) SignifGrad[which(pval.Grad ==
        pval.Grad.ord[i], arr.ind = TRUE)] <- TRUE
  }
  if (curv) {
    pval.Curv <- 1 - pchisq(WaldCurv, d * (d + 1)/2)
    pval.Curv.ord <- pval.Curv[order(pval.Curv)]
    num.test <- sum(!is.na(pval.Curv.ord))
    if (num.test >= 1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) -
          num.test))
    else num.test.seq <- rep(NA, prod(gridsize))
    reject.nonzero <- ((pval.Curv.ord <= signifLevel/(num.test +
        1 - num.test.seq)) & (pval.Curv.ord > 0))
    reject.nonzero.ind <- which(reject.nonzero)
    SignifCurv <- array(FALSE, dim = gridsize)
    SignifCurv[which(pval.Curv == 0, arr.ind = TRUE)] <- TRUE
    for (i in reject.nonzero.ind) SignifCurv[which(pval.Curv ==
        pval.Curv.ord[i], arr.ind = TRUE)] <- TRUE
    if (neg.curv.only)
      SignifCurv <- SignifCurv & local.mode
  }
  if (grad & !curv)
    return(list(grad = SignifGrad))
  else if (!grad & curv)
    return(list(curv = SignifCurv))
  else if (grad & curv)
    return(list(grad = SignifGrad, curv = SignifCurv))
}

dfltBWrange <- function (x, tau) {
  d <- ncol(x)
  if (d == 1)
    x <- as.matrix(x)
  r <- 2
  cmb.fac.upp <- (4/((d + 2 * r + 2) * nrow(x)))^(1/(d + 2 *
      r + 4))
  r <- 0
  cmb.fac.low <- (4/((d + 2 * r + 2) * nrow(x)))^(1/(d + 2 *
      r + 4))
  st.devs <- apply(x, 2, sd)
  IQR.vals <- apply(x, 2, IQR)/(qnorm(3/4) - qnorm(1/4))
  sig.hats <- apply(cbind(st.devs, IQR.vals), 1, min)
  range.h <- list()
  for (id in 1:d) {
    h.upp <- cmb.fac.upp * sig.hats[id]
    h.low <- 0.1 * cmb.fac.low * sig.hats[id]
    range.h[[id]] <- c(h.low, h.upp)
  }
  return(range.h)
}

featureSignif <- function (x, bw, gridsize, scaleData = FALSE, addSignifGrad = TRUE,
    addSignifCurv = TRUE, signifLevel = 0.05)
{
    tau <- 5
    if (is.vector(x)) {
        d <- 1
        n <- length(x)
        names.x <- deparse(substitute(x))
        if (scaleData)
            x <- (x - min(x))/(max(x) - min(x))
    }
    else {
        d <- ncol(x)
        n <- nrow(x)
        names.x <- colnames(x)
        if (is.null(names.x)) {
            names.xx <- deparse(substitute(x))
            names.xx <- strsplit(names.xx, "\\[")[[1]][1]
            names.x <- paste(names.xx, "[,", 1:d, "]", sep = "")
        }
        if (scaleData)
            for (i in 1:d) x[, i] <- (x[, i] - min(x[, i]))/(max(x[,
                i]) - min(x[, i]))
    }
    x <- as.matrix(x)
    if (d > 4)
        stop("Feature significance currently only available for 1- to 4-dimensional data")
    if (missing(gridsize)) {
        if (d == 1)
            gridsize <- 401
        if (d == 2)
            gridsize <- rep(151, 2)
        if (d == 3)
            gridsize <- rep(51, 3)
        if (d == 4)
            gridsize <- rep(21, 4)
    }
    if (missing(bw)) {
        bw.range <- dfltBWrange(as.matrix(x), tau)
        bw <- matrix(unlist(bw.range), nrow = 2, byrow = FALSE)
        dfltCounts.out <- dfltCounts(x, gridsize, apply(bw, 2,
            max))
        h.low <- bw[1, ]
        h.upp <- bw[2, ]
        hmix.prop <- 1/4
        h.init <- h.low^(hmix.prop) * h.upp^(1 - hmix.prop)
        h <- h.init
    }
    else {
        dfltCounts.out <- dfltCounts(x, gridsize, bw)
        h <- bw
    }
    gcounts <- dfltCounts.out$counts
    range.x <- dfltCounts.out$range.x
    dest <- drvkde(gcounts, drv = rep(0, d), bandwidth = h, binned = TRUE,
        range.x = range.x, se = FALSE, gridsize = gridsize)
    dest$est[dest$est < 0] <- 0
    SignifFeatureRegion.mat <- SignifFeatureRegion(n, d, gcounts,
        gridsize, dest, h, signifLevel, range.x, grad = addSignifGrad,
        curv = addSignifCurv)
    ESS <- n * dest$est * prod(h) * (sqrt(2 * pi)^d)
    SigESS <- ESS >= 5
    SignifGradRegion.mat <- SignifFeatureRegion.mat$grad & SigESS
    SignifGradData.mat <- SignifFeatureData(x, d, dest, SignifGradRegion.mat)
    SignifGradDataPoints <- x[SignifGradData.mat, ]
    SignifCurvRegion.mat <- SignifFeatureRegion.mat$curv & SigESS
    SignifCurvData.mat <- SignifFeatureData(x, d, dest, SignifCurvRegion.mat)
    SignifCurvDataPoints <- x[SignifCurvData.mat, ]
    feat <- c(list(x = x, names = names.x, bw = h, fhat = dest),
        SignifFeatureRegion.mat, list(gradData = SignifGradData.mat,
            gradDataPoints = SignifGradDataPoints, curvData = SignifCurvData.mat,
            curvDataPoints = SignifCurvDataPoints))
    class(feat) <- "fs"
    return(feat)
}
