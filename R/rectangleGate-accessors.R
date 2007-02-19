## ==========================================================================
## filter flowFrame object using rectangleGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature(x="flowFrame",table="rectangleGate"),
          function(x,table) {
            e <- if(length(table@parameters)==1) as.matrix(exprs(x)[,table@parameters]) else
            exprs(x)[,table@parameters]
            apply(sapply(seq(along=table@parameters), function(i) {
              !is.na(cut(e[,i],c(table@min[i],table@max[i]),labels=FALSE,
                         right=FALSE))
            }),1,all)
          })


## ==========================================================================
## show method for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature(object="rectangleGate"),function(object) {
	cat("Rectangular gate with dimensions:\n")
	for(i in seq(along=object@parameters)) {
		cat("\t")
		cat(object@parameters[i])
		cat(": (")
		cat(paste(object@min[i],object@max[i],sep=","))
		cat(")\n")
	}
})


## ==========================================================================
## summary method for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",signature(object="rectangleGate"),function(object,...) {
	structure(sapply(seq(along=object@parameters),
                         function(i) c(min=object@min[i],max=object@max[i])),
                  names=object@parameters)
})

