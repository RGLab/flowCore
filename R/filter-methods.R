## last changed November 2, 2006
## pdh: rectangleGate was still using cytoFrame. I changed it to fcsFrame
##      made the filterDetails more uniform among the methods
##      got the right filename for a fcsFrame to go in the msg



## ==========================================================================
## Apply polygon filter on fcsFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="polygonGate",flowSet="fcsFrame",parent="ANY"),
          definition=function(filter,flowSet,parent) {
              selectNew <-polygonFiltering(filter,flowSet,parent)
              msg <- paste("polygonGate applied on ",
                           deparse(substitute(flowSet)),
                           " (file:",basename(flowSet@description["$FIL"]),
                           ") a ", class(flowSet)," object", sep="")
       			out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          } 
          )
## ==========================================================================



  

## ==========================================================================
##  Apply rectangular filter on fcsFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="rectangleGate",flowSet="fcsFrame",parent="ANY"),
          definition=function(filter,flowSet,parent) {
              selectNew <-rectangleFiltering(filter,flowSet,parent)
               msg <- paste("rectangleGate applied on ",
                            deparse(substitute(flowSet)),
                            " (file:",basename(flowSet@description["$FIL"]),
                            ") a ", class(flowSet)," object", sep="")
               out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          })
## ==========================================================================




## ==========================================================================
##  Apply norm2Filter on fcsFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="norm2Filter",flowSet="fcsFrame",parent="ANY"),
          definition=function(filter,flowSet,parent) {
              selectNew <-normFiltering(filter,flowSet,parent)
              msg <- paste("norm2Filter applied on ",
                            deparse(substitute(flowSet)),
                            " (file:",basename(flowSet@description["$FIL"]),
                            ") a ", class(flowSet)," object", sep="")
               out = list(msg,filter=filter,mu=selectNew$mu,S=selectNew$S)
              new("filterResult", subSet=selectNew$sel,filterDetails=out)
          })
## ==========================================================================

  
