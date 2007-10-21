## =========================================================================##
## =========================================================================##
##                       Methods for curv2Filter object                     ##
## =========================================================================##
## =========================================================================##




## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("curv2Filter"),function(object) {
	msg = paste("A curv2 filter named '",object@filterId,"' with settings:",
        "\n  bwFac=", object@bwFac, "\n  gridsize=",
        paste(object@gridsize, collapse=",", sep=""), sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})



## ==========================================================================
## Filtering Method -- we are not a logical filter so we return a vector
## of indices indicating a population. Additional information about the
## filter result (polygon vertices of populations, fSObj) are stored as
## attributes of the subSet vector.
## ---------------------------------------------------------------------------
setMethod("%in%",signature("flowFrame","curv2Filter"),function(x,table)
      {
          ## We accomplish the actual filtering via Matt Wands feature software
          param <- table@parameters
          values <- exprs(x)[, param]
          bwFac <- table@bwFac
          gridsize <- table@gridsize

          ## Compute normal scale bandwidths.
          st.devs <- sqrt(apply(values, 2, var))
          Q1.vals <- apply(values, 2, quantile, 1/4)
          Q3.vals <- apply(values, 2, quantile, 3/4)
          corr.fac <- qnorm(3/4) - qnorm(1/4)
          IQR.vals <- (Q3.vals - Q1.vals)/corr.fac
          sig.hats <- apply(cbind(st.devs, IQR.vals), 1, min)
          samp.size.fac <- nrow(values)^(-1/6)
          bwNS <- samp.size.fac*sig.hats
          
          ## Obtain significant high curvature regions.
          fSObj <- featureSignif(values, bw=bwFac*bwNS, addSignifCurvRegion=TRUE,
                           gridsize=gridsize, plotFS=FALSE)
          contourLinesObj <- contourLines(fSObj$fhat$x[[1]], fSObj$fhat$x[[2]],
                                    fSObj$curv, levels=0.5)

          ## Determine filter member indicator
          filterInds <- rep(0,nrow(values))
          for (i in seq(along=contourLinesObj)){
              vertices <- cbind(contourLinesObj[[i]]$x, contourLinesObj[[i]]$y)
              filterInds[as.logical(flowCore:::inpolygon(values,vertices))] <- i
          }
          
          result <- factor(filterInds)
          attr(result,'polygons') = contourLinesObj
          attr(result,'fSObj')    = fSObj
          result
})



## ==========================================================================
## summarize results of a curv2 filtering operation
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",signature("filterResult","curv2Filter"),
          function(result,filter) {
	ret = callNextMethod()
	ret$polygons = attr(result@subSet, "polygon")
	ret
})
