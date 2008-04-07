## ==========================================================================
## Filtering Methods -- we are not a logical filter so we return a vector
## of factors indicating a population.
## ---------------------------------------------------------------------------
setMethod("%in%",signature(x="flowFrame",table="quadGate"),
          function(x,table) {
              e <-  exprs(x)[,table@parameters, drop=FALSE]
              lev <- c("tr", "tl", "br", "bl")
              factor(lev[as.integer(e[,1] <= table@boundary[1]) +
                         2 * (as.integer(e[,2] <= table@boundary[2]))+1])
          })


setMethod("summarizeFilter",signature("filterResult","quadGate"),
          function(result,filter) {
	ret = callNextMethod()
	ret$populations = levels(result@subSet)
	ret
})



## ==========================================================================
## show method for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature(object="quadGate"),function(object) {
	cat("Quadrant gate with dimensions:\n")
	for(i in seq(along=object@parameters)) {
		cat("\t")
		cat(object@parameters[i])
		cat(": ")
		cat(object@boundary[i])
		cat("\n")
	}
})

