## last changed November 2, 2006
## pdh: rectangleGate was still using cytoFrame. I changed it to flowFrame
##      made the filterDetails more uniform among the methods
##      got the right filename for a flowFrame to go in the msg

## ==========================================================================
## Apply polygon filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="polygonGate",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-polygonFiltering(filter,flowObject,parent)
              msg <- paste("polygonGate applied on ",
                           deparse(substitute(flowObject)),
                           " (file:",basename(flowObject@description["$FIL"]),
                           ") a ", class(flowObject)," object", sep="")
       			out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          } 
          )
## ==========================================================================



  

## ==========================================================================
##  Apply rectangular filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="rectangleGate",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-rectangleFiltering(filter,flowObject,parent)
               msg <- paste("rectangleGate applied on ",
                            deparse(substitute(flowObject)),
                            " (file:",basename(flowObject@description["$FIL"]),
                            ") a ", class(flowObject)," object", sep="")
               out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          })
## ==========================================================================




## ==========================================================================
##  Apply norm2Filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="norm2Filter",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-normFiltering(filter,flowObject,parent)
              msg <- paste("norm2Filter applied on ",
                            deparse(substitute(flowObject)),
                            " (file:",basename(flowObject@description["$FIL"]),
                            ") a ", class(flowObject)," object", sep="")
               out = list(msg,filter=filter,mu=selectNew$mu,S=selectNew$S)
              new("filterResult", subSet=selectNew$sel,filterDetails=out)
          })
## ==========================================================================

  
