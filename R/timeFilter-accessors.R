## =========================================================================##
## =========================================================================##
##                       Methods for curv1Filter object                     ##
## =========================================================================##
## =========================================================================##




## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("timeFilter"),function(object)
          {
            msg <- paste("A time filter named '",object@filterId,
                         "' with settings:\n  bandwidth=",
                         object@bandwidth, sep="")
            cat(msg)
            if(length(object@binSize))
              cat("\n  binSize=", object@binSize, sep="")
            if(length(object@timeParameter))
              cat("\n  timeParameter=", object@timeParameter, sep="")
            cat("\n")
            invisible(msg)
          })



## ==========================================================================
## Filtering Method -- Strictly, this is not a logical filter because
## there might be multiple stretches over time that behave strangly, however
## we treat it as one since there is not real use case for keeping them
## separate. If one wants to restrict the time filtering it can always
## be combined with rectangleGates...
## ---------------------------------------------------------------------------
setMethod("%in%",signature("flowFrame","timeFilter"),function(x,table)
      {
          ## We first bin the data and compute summary statistics
          ## for each bin. These are then used to identify stretches
          ## of unusual data distribution over time. The time parameter
          ## has to be guessed if not explicitely given before.
          ex <- exprs(x)
          if(!length(table@timeParameter))
              time <- findTimeChannel(ex) 
          param <- parameters(table)
          bw <- table@bandwidth
          bs <- table@binSize
          if(!length(bs))
              bs <- min(max(1, floor(nrow(x)/100)), 500)
          binned <- prepareSet(ex, param, time, bs,
                               locM=median, varM=mad)
          
          ## Standardize to compute meaningful scale-free scores.
          ## This is done by substracing the mean values of each
          ## bin and divide by the mean bin variance.
          med <- median(binned$smooth[,2], na.rm=TRUE)
          gvars <- mean(binned$variance)
          stand <- abs(binned$smooth[,2]-med)/(gvars*bw)
          outBins <- which(stand > 1)
          bins <- c(-Inf, binned$bins, Inf)

          ## we can treat adjacend regions as one
          ## FIXME: There must be a more elegant way to do that
          if(length(outBins)){
              ld <- length(outBins)
              db <- c(diff(c(0, outBins)),2)
              br <- tr <- NULL
              adj <- FALSE
              for(i in seq_len(ld)){
                  adj <- db[i]==1
                  if(!adj){
                      br <- c(br, outBins[i])
                      if(db[i+1]!=1)
                          tr <- c(tr, outBins[i]+1)
                  }else if(db[i+1]!=1)
                      tr <- c(tr, outBins[i]+1) 
              }
              if(db[1]==1)
                  br <- c(outBins[1], br)
            
              ## Now generate rectangle gates over the identified
              ## regions and use them for the filtering
              ## FIXME: This step is notoriously slow. Is there a
              ## faster way to do that?
              gates <- mapply(function(l,r)
                              rectangleGate(.gate=matrix(c(l,r),ncol=1,
                                              dimnames=list(NULL, time))),
                              bins[pmax(1, br-1)], bins[pmin(length(bins),
                                                             tr+2)])
              tmp <- filter(x, !gates[[1]])
              if(length(gates)>1)
                  for(i in 2:length(gates))
                      tmp <- filter(x, tmp & !gates[[i]])
              return(as.logical(tmp@subSet) )
          }else
              return(rep(TRUE, nrow(x)))   
        })




## FIXME: The following was copied from flowViz, eventually it
## should live only here:
## Run over a cytoFrame, bin the values according to the time domain
## and compute a location measure locM as well as variances varM for each bin.
## The result of this function will be the input to the plotting functions
## and the basis for the quality score.
prepareSet <- function(x, parm, time, binSize, locM=median, varM=mad){
    xx <- x[, time]
    ord <- order(xx)
    xx <- xx[ord]
    yy <- x[ord, parm]
    lenx <- length(xx)
    nrBins <- floor(lenx/binSize)
    ## how many events per time tick
    nr <- min(length(unique(xx)), max(51, nrBins))
    timeRange <- seq(min(xx), max(xx), len=nr)
    hh <- hist(xx, timeRange, plot = FALSE)
    freq <- hh$counts
    expEv <- length(xx)/(nr-1)
    ## time parameter is already binned or very sparse events
    ux <- unique(xx)
    if(length(ux) < nrBins){
        tmpy <- split(yy, xx)
        yy <- sapply(tmpy, locM, na.rm=TRUE)
        xx <- unique(xx)
        binSize <- 1
    }else{
        ## bin values in nrBins bins
        if(lenx > binSize){
            cf <- c(rep(1:nrBins, each=binSize),
                    rep(nrBins+1, lenx-nrBins*binSize))
            stopifnot(length(cf) == lenx)
            tmpx <- split(xx,cf)
            tmpy <- split(yy,cf)
            yy <- sapply(tmpy, locM, na.rm=TRUE)
            xx <- sapply(tmpx, mean, na.rm=TRUE)
        }else{
            ## very little events
            warning("Low number of events", call.=FALSE)
            tmpy <- split(yy,xx)
            yy <- sapply(tmpy, locM, na.rm=TRUE)
            xx <- unique(xx)
            binSize <- 1
        }
    }
    var <- sapply(tmpy, varM, na.rm=TRUE)
    ## avoid 0 variance estimates created by mad
    zv <- which(var==0)
    if(length(zv))
       var[zv] <- mean(sapply(tmpy, sd, na.rm=TRUE), na.rm=TRUE)
    return(list(smooth=cbind(xx,yy), variance=var, binSize=binSize,
                frequencies=cbind(timeRange[-1], freq),
                expFrequency=expEv, bins=unique(xx)))
}


findTimeChannel <- function(xx)
{
    time <- grep("^Time$", colnames(xx), value=TRUE,
                 ignore.case=TRUE)[1]
    if(!length(time)){
        cont <- apply(xx, 2, function(y) all(sign(diff(y)) >= 0))
        time <- names(which(cont))
    }
    if(!length(time))
        stop("No time domain recording for this data.\n",
             "Please define time parameter in the filter's",
             "'timeParameter' slot.", call.=FALSE)
    return(time)
}

