## Last changed November 3, 2006
##  pdh -- fixed an error because I wasn't checking for xlim and ylim when y was missing.

## ==========================================================================
## Basic plot for fcsFrame object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
          ## basic plot without a gate specified
          signature(x="flowFrame",y="missing"),
          function(x, y, filter=NULL, plotParameters=c("FSC-H","SSC-H"),parent,colParent="Grey",colSubSet="Blue",xlab,ylab,xlim,ylim,...){
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
            if(missing(y)) {
              ## check for the parent. If it is a filterResult get the parentSet
              parentSet <-  getParentSet(parent,ncells)
              selectParent <- (parentSet==1)
              if(missing(xlim)) {
                xlim <- range(data[selectParent,plotParameters[1]])
              }
              if(missing(ylim)) {
                ylim <-  range(data[selectParent,plotParameters[2]])
              }
              
              smoothScatter(data[selectParent,plotParameters[1]],data[selectParent,plotParameters[2]],nrpoints=50,
                            xlim=xlim,
                            ylim=ylim,
                            xlab=xlab,
                            ylab=ylab,
                            ...)       
            }
            else {
              ## check for y and if it is a filterResult, get the subSet
              if(class(y) != "filterResult") {
                stop("y must be of class filterResult.")
              }
              else {
                if(is.null(y@subSet)) {
                  stop("The filterResult y must have a subSet slot.")
                }
                else {
                  subSet <-  y@subSet
                }
              }
              ## check for the parent. If it is a filterResult get the parentSet
              parentSet <- getParentSet(parent,ncells)
              
              ## check for the proper relationship between the parent and child
              if(any(subSet>parentSet)) {
                stop("The child population must be contained in the parent population.")
              }
              selectParent <-  (parentSet==1)
              selectSubSet <-  (subSet==1)
              selectUnfiltered <- ((parentSet-subSet) ==1)
              if(missing(xlim)) {
                xlim = range(data[selectParent,plotParameters[1]])
              }
              if(missing(ylim)) {
                ylim = range(data[selectParent,plotParameters[2]])
              }
              
              plot(data[selectSubSet,plotParameters[1]],data[selectSubSet,plotParameters[2]],
                   xlab=xlab,
                   ylab=ylab,
                   xlim=xlim,
                   ylim=ylim,
                   type="n",
                   ...)
              ## filtered population
              points(data[selectSubSet,plotParameters[1]],data[selectSubSet,plotParameters[2]],col=colSubSet,...)
              ## unfiltered population
              points(data[selectUnfiltered,plotParameters[1]],data[selectUnfiltered,plotParameters[2]],col=colParent,...)
              if(length(filter)!=0){
                ln <- filter@boundaries
                lines(rbind(ln,ln[1,]),col="green",lwd=2)  
              }
            }
          })
