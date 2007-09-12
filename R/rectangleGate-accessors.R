## ==========================================================================
## filter flowFrame object using rectangleGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature(x="flowFrame",table="rectangleGate"),
          function(x,table) {
             
            e <- if(length(table@parameters)==1) as.matrix(exprs(x)[,table@parameters]) else
            exprs(x)[,table@parameters]
            apply(sapply(seq(along=table@parameters), function(i) {
                if(table@min[i] >= table@max[i]){
                    print(i)
                  e[,i]
                } else {
              !is.na(cut(e[,i],c(table@min[i],table@max[i]),labels=FALSE,
                         right=FALSE))}
            }), 1, all)
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


# Compose two rectangle gates together into a higher dimensional cube.
setMethod("*",signature(e1="rectangleGate",e2="rectangleGate"),function(e1,e2) {
	if(any(parameters(e1) %in% parameters(e2))) stop("Rectangle gate parameters overlap.")
	new("rectangleGate",parameters=c(parameters(e1),parameters(e2)),min=c(e1@min,e2@min),max=c(e1@max,e2@max))
})