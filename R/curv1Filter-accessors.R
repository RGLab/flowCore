## =========================================================================##
## =========================================================================##
##                       Methods for curv1Filter object                     ##
## =========================================================================##
## =========================================================================##




## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("curv1Filter"),function(object) {
	msg = paste("A curv1 filter named '",object@filterId,"' with settings:",
        "\n  bwFac=", object@bwFac, "\n  gridsize=",
        paste(object@gridsize, collapse=",", sep=""), sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})



## ==========================================================================
## Filtering Method -- we are not a logical filter so we return a vector
## of indices indicating a population. Additional information about the
## filter result (boundaries of regions, fSObj) are stored as
## attributes of the subSet vector.
## ---------------------------------------------------------------------------
setMethod("%in%",signature("flowFrame","curv1Filter"),function(x,table)
      {
          ## We accomplish the actual filtering via Matt Wands feature software
          param <- table@parameters
          values <- exprs(x)[, param]
          bwFac <- table@bwFac
          gridsize <- table@gridsize

          ## Compute normal scale bandwidth (second derivative).
          st.dev <- sqrt(var(values))
          Q1.val <- quantile(values,1/4) ; Q3.val <- quantile(values,3/4)
          IQR.val <- (Q3.val - Q1.val)/(qnorm(3/4) - qnorm(1/4))
          bwNS <- min(st.dev,IQR.val)*(4/(7*length(values)))^(1/9)
          
          ## Obtain significant high curvature intervals.
          fSObj <- featureSignif(values, bw=bwFac*bwNS,
                                 addSignifCurvRegion=TRUE,
                                 gridsize=gridsize, plotFS=FALSE)
          xGrid <- unlist(fSObj$fhat$x.grid)
          hiCurvIndic <- as.numeric(fSObj$curv)
          diffGrid <- diff(c(0,hiCurvIndic,0))
          lowInds <- (1:length(diffGrid))[diffGrid==1]
          uppInds <- (1:length(diffGrid))[diffGrid==-1]
          lowLims <- (xGrid[lowInds] + xGrid[lowInds-1])/2
          uppLims <- (xGrid[uppInds] + xGrid[uppInds-1])/2
          lims <- lapply(1:length(lowLims), function(i)
              c(lowLims[i],uppLims[i]))
          
          ## Determine filter member indicator
          indices <- rep(0, length(values))
          for(i in seq(along=lims))
              indices[values>=lims[[i]][1] & values <= lims[[i]][2]] <- i
          result <- factor(indices)
          attr(result,'boundaries') <- lims
          attr(result,'fSObj') <- fSObj
          result
})



## ==========================================================================
## summarize results of a curv2 filtering operation
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",signature("filterResult","curv1Filter"),
          function(result,filter) {
	ret = callNextMethod()
	ret$boundaries = attr(result@subSet, "boundaries")
	ret
})
