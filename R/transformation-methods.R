## NLM Jan 17
## ==========================================================================
## Transformation function for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...) {
              e <- substitute(list(...))
              x = `_data`
              transformed <- as.matrix(transform(as.data.frame(exprs(x)),...))
              ##Add any new parameter values
              if(ncol(transformed) > ncol(x@exprs)) {
              	cnames = c(colnames(x),colnames(transformed)[-c(1:ncol(x@exprs))])
              }
              else {
              	cnames = colnames(x)
              }
              ##param.names <- colnames(transformed)
              ##newParams <- is.na(match(param.names,`_data`@parameters$name))
              ##params      = parameters(`_data`)$name
              ##if(any(newParams)) {
              ##    params <- cbind(params,
              ##                    data.frame(name=param.names))
                  ##params = cbind(params,
                  ##    data.frame(name=param.names[newParams]))
              ##}
              #colnames(transformed) <- parameters(`_data`)$name
              colnames(transformed) = cnames
              new("flowFrame",
                  exprs=transformed, 
                  parameters=parameters(x),#[,params$name],
                  description=description(x))
          })
## ==========================================================================
## Transform function for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",signature=signature(`_data`="flowSet"),function(`_data`,...) {
	fsApply(`_data`,transform,...)
})

setMethod("transform",signature(`_data`="missing"),function(...) {
	funs = list(...)
	io   = names(funs)
	#Consistency check
	if(!all(sapply(funs,is.function)))
		stop("All transforms must be functions")
	if(!all(sapply(io,is.character)))
		stop("All transforms must be named")
	new("transformList",transforms=lapply(seq(along=funs),function(i) new("transformMap",input=io[i],output=io[i],f=funs[[i]])))
})

setMethod("colnames",signature("transformList"),function(x, do.NULL=TRUE, prefix="col") {
	unique(sapply(x@transforms,slot,"input"))
})

setMethod("%on%",signature("transformList","flowFrame"),function(e1,e2) {
	x = exprs(e2)
	for(y in e1@transforms){
		x[,y@input] = y@f(x[,y@output])
	}
	exprs(e2) = x
	e2
})
setMethod("%on%",signature("transformList","flowSet"),function(e1,e2) fsApply(e2,"%on%",e1=e1))

#So that we can get parameters back OUT of a transform
setMethod("summary",signature("transform"),function(object,..) {
	e = environment(object)
	x = ls(env=e)
	structure(lapply(x,"get",env=e),names=x)
})

