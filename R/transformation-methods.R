
## ==========================================================================
## Apply linear transformation on fcsframe Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="linearTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <-linearTransformation(transformation,flowObject)
             msg <- paste("Linear transformation applied on ",
                           deparse(substitute(flowset)),
                           " (file:",deparse(substitute(flowObject)),
                           ") a ", class(flowObject)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply quadratic transformation on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="quadraticTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <-quadraticTransformation(transformation,flowObject)
             msg <- paste("Quadratic transformation applied on ",
                           deparse(substitute(flowObject)),
                           " (file:",deparse(substitute(flowObject)),
                           ") a ", class(flowObject)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply LN transformation on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="lnTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <-lnTransformation(transformation,flowObject)
             msg <- paste("Natural Logarithmical (LN) transformation applied on ",
                           deparse(substitute(flowObject)),
                           " (file:", deparse(substitute(flowObject)),
                           ") a ", class(flowObject)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )
## ==========================================================================
## Apply LOG transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="logTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <-logTransformation(transformation,flowObject)
             msg <- paste("Logarithmical transformation applied on ",
                           deparse(substitute(flowObject)),
                           " (file:", deparse(substitute(flowObject)),
                           ") a ", class(flowObject)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )
## ==========================================================================
## Apply HyperLog transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="hyperlogTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <-hyperlogTransformation(transformation,flowObject)
             msg <- paste("hyperlog transformation applied on ",
                           deparse(substitute(flowObject)),
                           " (file:", deparse(substitute(flowObject)),
                           ") a ", class(flowObject)," object", sep="")
            print(msg)
            return(transformNew)
          }
          
 )

## ==========================================================================
## Apply biexponential transformation on FCS Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyTransformation",
          signature=signature(transformation="biexponentialTransformation",flowObject="flowFrame"),
          definition=function(transformation,flowObject) {
            transformNew <- biexponentialTransformation(transformation,flowObject)
            return(transformNew)
          }
          )
