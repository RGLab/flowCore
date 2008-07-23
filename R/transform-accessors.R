setMethod("%on%",signature("filter","transform"),function(e1,e2) {
	new("transformFilter",
		filterId=paste(identifier(e1),"on transformed values"),
            transforms=parameterTransform(e2,parameters(e1)),
            filter=e1,parameters=parameters(e1))
})

setMethod("%on%",signature("transform","flowFrame"),function(e1,e2) {
	exprs(e2) = e1(exprs(e2))
	e2
})

setMethod("%on%",signature("filter","parameterTransform"),function(e1,e2) {
	new("transformFilter",
		filterId=paste(e1@filterId," on transformed values of ",paste(e2@parameters,sep=","),collapse=" "),
            transforms=e2,
            filter=e1,parameters=unique(c(e1@parameters,e2@parameters)))
})
setMethod("%on%",signature("parameterTransform","flowFrame"),function(e1,e2) {
	params = if(length(e1@parameters)) colnames(e2) else e1@parameters
	exprs(e2)[,params] = e1(exprs(e2))[,params]
	e2
})


## ==========================================================================
## summary method
## So that we can get parameters back OUT of a transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",
          signature("transform"),
          function(object, ...) {
              e <- environment(object)
              x <- ls(env=e)
              structure(lapply(x, "get", env=e), names=x)
          })

