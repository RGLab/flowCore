## ==========================================================================
## Transformation function for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...) {
              flowObject <- `_data`
              e <- substitute(list(...))
              transformed <- transform(as.data.frame(flowObject@exprs),...)
              newFlowFrame <- new("flowFrame", exprs=as.matrix(transformed),
                                  description=flowObject@description)
              newFlowFrame@description <- c(newFlowFrame@description,"transformation"=e)
              return(newFlowFrame)
          }
          
 )
