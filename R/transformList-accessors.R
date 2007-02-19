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
	for(y in e1@transforms){
		x[,y@input] = y@f(x[,y@output])
	}
	exprs(e2) = x
	e2
})
setMethod("%on%",signature("transformList","flowSet"),function(e1,e2) fsApply(e2,"%on%",e1=e1))


## ==========================================================================
## summary method
## So that we can get parameters back OUT of a transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",signature("transform"),function(object,..) {
	e = environment(object)
	x = ls(env=e)
	structure(lapply(x,"get",env=e),names=x)
})

