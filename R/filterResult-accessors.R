#We can convert a factor, logical or a numeric into a filterResult by selecting a 
#secific filterResult type. This is done through the standard R coercion techniques.
setAs("factor", "filterResult",function(from) new("multipleFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("logical","filterResult",function(from) new("logicalFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("numeric","filterResult",function(from) new("randomFilterResult",parameters=character(0),filterId="",subSet=from))

setAs("filterResult","logical",function(from) stop("Unable to convert to a logical vector"))
setAs("logicalFilterResult","logical",function(from) from@subSet)
setAs("randomFilterResult","logical",function(from) runif(length(from@subSet))<from@subSet)

##Allow us to compare filterResults and flowFrames. This lets us check (and warn or stop)
##that a particular flowFrame generated a filterResult allowing us to use it for further processing.
setMethod("identifier", signature="filterResult",
          definition=function(object) object@frameId)
setMethod("==",signature("flowFrame","filterResult"),
          definition=function(e1,e2) {
            i1 = identifer(e1)
            i2 = identifier(e2)
            (length(i1) == 0 || length(i2) == 0 || i1 == i2)
          })
#Does S4 do this for us automagically? I don't know!
setMethod("==",signature("filterResult","flowFrame"),
          definition=function(e1,e2) e2==e1)





## Last changed November 3, 2006
##  pdh -- fixed an error because I wasn't checking for xlim and ylim when y was missing.

## ==========================================================================
## Basic plot for fcsFrame object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
          ## basic plot without a gate specified
          signature(x="flowFrame"),
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
