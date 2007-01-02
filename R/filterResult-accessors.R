#We can convert a factor, logical or a numeric into a filterResult by selecting a 
#secific filterResult type. This is done through the standard R coercion techniques.
setAs("factor", "filterResult",function(from) new("multipleFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("logical","filterResult",function(from) new("logicalFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("numeric","filterResult",function(from) new("randomFilterResult",parameters=character(0),filterId="",subSet=from))

setAs("filterResult","logical",function(from) stop("Unable to convert to a logical vector"))
setAs("logicalFilterResult","logical",function(from) from@subSet)
setAs("randomFilterResult","logical",function(from) runif(length(from@subSet))<from@subSet)
setMethod("%in%",signature("ANY","filterResult"),function(x,table) {
	if(x != table) stop("filterResult doesn't match left-hand side.")
	as(table,"logical")
})


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



