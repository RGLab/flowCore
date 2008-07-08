## ==========================================================================
## Filtering Methods -- we are not a logical filter so we return a vector
## of factors indicating a population.
## ---------------------------------------------------------------------------
setMethod("%in%",signature(x="flowFrame",table="quadGate"),
          function(x,table) {
              e <-  exprs(x)[,table@parameters, drop=FALSE]
              ## create useful names for the populations:
              ## first choice are channel descriptors, then
              ## channel names.
              mt <- match(table@parameters, parameters(x)$name)
              desc <- pData(parameters(x))[mt, "desc"]
              noName <- which(is.na(desc) | desc=="")
              desc[noName] <- parameters(x)$name[mt][noName]
              lev <- c(sprintf("%s+%s+", desc[1], desc[2]),
                       sprintf("%s-%s+", desc[1], desc[2]),
                       sprintf("%s+%s-", desc[1], desc[2]),
                       sprintf("%s-%s-", desc[1], desc[2]))
              #lev <- c("tr", "tl", "br", "bl")
              factor(lev[as.integer(e[,1] <= table@boundary[1]) +
                         2 * (as.integer(e[,2] <= table@boundary[2]))+1],
                     levels=lev)
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

