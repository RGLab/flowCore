
## ==========================================================================
## Apply linear transformation on fcsframe Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="linearTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <-linearTransformation(transformation,flowSet)
             msg <- paste("Linear transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:",deparse(substitute(flowset)),
                           ") a ", class(flowSet)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply quadratic transformation on fcsFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="quadraticTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <-quadraticTransformation(transformation,flowSet)
             msg <- paste("Quadratic transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:",deparse(substitute(flowset)),
                           ") a ", class(flowSet)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply LN transformation on fcsFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="lnTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <-lnTransformation(transformation,flowSet)
             msg <- paste("Natural Logarithmical (LN) transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:", deparse(substitute(flowSet)),
                           ") a ", class(flowSet)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )
## ==========================================================================
## Apply LOG transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="logTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <-logTransformation(transformation,flowSet)
             msg <- paste("Logarithmical transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:", deparse(substitute(flowSet)),
                           ") a ", class(flowSet)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )
## ==========================================================================
## Apply HyperLog transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="hyperlogTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <-hyperlogTransformation(transformation,flowSet)
             msg <- paste("hyperlog transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:", deparse(substitute(flowSet)),
                           ") a ", class(flowSet)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply biexponential transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="biexponentialTransformation",flowSet="fcsFrame"),
          definition=function(transformation,flowSet) {
            transformNew <- biexponentialTransformation(transformation,flowSet)
            return(transformNew)
          }
          )
