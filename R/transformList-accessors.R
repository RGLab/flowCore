## ==========================================================================
## colnames method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",signature("transformList"),function(x, do.NULL=TRUE, prefix="col") {
	unique(sapply(x@transforms,slot,"input"))
})


## ==========================================================================
## %on% operators
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%on%",signature("transformList","flowFrame"),function(e1,e2) {
	x = exprs(e2)
	cN = colnames(x)
	for(y in e1@transforms){
		if( !(y@output %in% cN) )
			stop(y@output, "is not a variable in the flowFrame")
		x[,y@output] = y@f(x[,y@input])
	}
	exprs(e2) = x
	e2
})
#General %on% implementation for a flowSet.
setMethod("%on%",signature(e2="flowSet"),function(e1,e2) fsApply(e2,"%on%",e1=e1))


## ==========================================================================
## summary method
## So that we can get parameters back OUT of a transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",signature("transform"),function(object,..) {
	e = environment(object)
	x = ls(env=e)
	structure(lapply(x,"get",env=e),names=x)
})

setMethod("c",signature("transformList"),function(x,...,recursive=FALSE) {
	#Try to coerce everyone to a transformList first
	all.t = lapply(list(...),as,"transformList")
	params = c(sapply(x@transforms,slot,"output"),unlist(lapply(all.t,function(x) sapply(x@transforms,slot,"output"))))
	if(length(params) > length(unique(params)))
		stop("All output parameters must be unique when combining transforms")
	new("transformList",transforms=c(x@transforms,unlist(sapply(all.t,slot,"transforms"))))
})
