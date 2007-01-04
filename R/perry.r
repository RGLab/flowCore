## last changed January 2, 2007
##   
## I moved my plot method here as it was getting overwritten by the core plot method
## if the consensus is that this kind of convenience function doesn't belong in the core 
## package, I'll take it local
## Last changed November 3, 2006
##  pdh -- fixed an error because I wasn't checking for xlim and ylim when y was missing.

## ==========================================================================
## Basic plot for fcsFrame object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setGeneric("flowPlot",function(x,...) standardGeneric("flowPlot"))
setMethod("flowPlot",
          ## basic plot without a gate specified
          signature(x="flowFrame"),
          function(x, y, filter=NULL, plotParameters=c("FSC-H","SSC-H"),
          	parent,colParent="Grey",colSubSet="Blue",
          	showFilter=TRUE,gate.fill="transparent",gate.border="black",
          	xlab,ylab,xlim,ylim,...){
            data <- x@exprs
            ncells <- nrow(data)
            if(missing(xlab)){
              if(is.character(plotParameters)){
                xlab <-  plotParameters[1]
              }
              else {
                xlab <-  colnames(data)[plotParameters[1]]
              }
            }
            if(missing(ylab)){
              if(is.character(plotParameters)){
                ylab <-  plotParameters[2]
              }
              else {
                ylab <-  colnames(data)[plotParameters[2]]
              }
            }
            ## check for the parent. If it is a filterResult get the parentSet
            if(missing(parent)) {
   				## if the parent is missing assume that it is all cells
				selectParent = rep(TRUE,ncells)
 			}
			else {  	
				if(class(parent) != "logicalFilterResult") {
					stop("parent must be of class filterResult.")
				}
				else {
					if(is.null(parent@subSet)) {
						stop("The filterResult parent must have a subSet slot.")
					}
					else {
						selectParent = parent@subSet
					}  	
				}
			}
            ## the default limits are going to be 0,1 or the range of the data
            if(missing(xlim)) {
              	if(max(data[selectParent,plotParameters[1]]) <= 1.0) {
              		xlim = c(0,1)
              	}
              	else {
              		xlim <-  range(data[selectParent,plotParameters[1]])
              	}
              }
              if(missing(ylim)) {
              	if(max(data[selectParent,plotParameters[2]]) <= 1.0) {
              		ylim = c(0,1)
              	}
              	else {
                	ylim <-  range(data[selectParent,plotParameters[2]])
                }
              }
            if(missing(y)) {             
              
              smoothScatter(data[selectParent,plotParameters[1]],data[selectParent,plotParameters[2]],nrpoints=50,
                            xlim=xlim,
                            ylim=ylim,
                            xlab=xlab,
                            ylab=ylab,
                            ...)       
            }
            else {
              ## check for y and if it is a filterResult, get the subSet
              if(class(y) != "logicalFilterResult") {
                stop("y must be of class filterResult.")
              }
              else {
                if(is.null(y@subSet)) {
                  stop("The filterResult y must have a subSet slot.")
                }
                else {
                  selectChild <-  y@subSet
                }
              }
              ## check for the proper relationship between the parent and child
              if(any(selectChild>selectParent)) {
                stop("The child population must be contained in the parent population.")
              }
              selectUnfiltered <- ((selectParent & !selectChild))
              
              
              plot(data[selectChild,plotParameters[1]],data[selectChild,plotParameters[2]],
                   xlab=xlab,
                   ylab=ylab,
                   xlim=xlim,
                   ylim=ylim,
                   type="n",
                   ...)
              ## filtered population
              points(data[selectChild,plotParameters[1]],data[selectChild,plotParameters[2]],col=colSubSet,...)
              ## unfiltered population
              points(data[selectUnfiltered,plotParameters[1]],data[selectUnfiltered,plotParameters[2]],col=colParent,...)
              if(showFilter){
              	yf = y@filterDetails$filter
              	         #    		browser()
              	if(class(yf) == "subsetFilter") {
              		yf = yf@left
              	}
             	if(class(yf) == "rectangleGate") {
             		## a simple range gate
             		if(length(yf@parameters)==1){
             			if(yf@parameters==plotParameters[1]) {
             				abline(v=yf@min,col=gate.border)
             				abline(v=yf@max,col=gate.border)
             			}
             			else if(yf@parameters==plotParameters[2]) {
             				abline(h=yf@min,col=gate.border)
             				abline(h=yf@max,col=gate.border)
             			}
             		}
             		## the classic rectangle
             		else if(length(yf@parameters)==2) {
             			rectends = c(yf@min[plotParameters[1]],yf@min[plotParameters[2]],yf@max[plotParameters[1]],yf@max[plotParameters[2]])
             			## draw the rectangles out to the boundaries if -Inf or Inf appears

             			if(rectends[1] == -Inf) {
             				rectends[1] = xlim[1]
             			}
             			if(rectends[2] == -Inf) {
             				rectends[2] = ylim[1]
             			}
             			if(rectends[3] == Inf) {
             				rectends[3] = xlim[2]
             			}
             			if(rectends[4] == -Inf) {
             				rectends[4] = ylim[2]
             			}
             			rect(rectends[1],rectends[2],rectends[3],rectends[4],col=gate.fill,border=gate.border,...)	
             		}
				}
 
              }
            }
          })
