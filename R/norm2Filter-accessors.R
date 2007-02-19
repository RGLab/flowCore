## ==========================================================================
##  filter flowFrame object using norm2Filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature("flowFrame",table="norm2Filter"),function(x,table) {
	if(length(table@parameters) != 2)
		stop("norm2 filters require exactly two parameters.")
	y = if(length(table@transformation)>0) 
		exprs(do.call("transform",c(x,table@transformation)))[,table@parameters] else
        exprs(x)[,table@parameters]
	
	if(is.na(match(table@method,c("covMcd","cov.rob"))))
		stop("Method must be either 'covMcd' or 'cov.rob'")
	cov = switch(table@method,
		covMcd = {
			if(nrow(y)>table@n) covMcd(y[sample(nrow(y),table@n),]) else covMcd(y)
		},
		cov.rob={cov.rob(y)},
		stop("How did you get here?")
		)
	W  = t(y)-cov$center
	result = exp(-.5*colSums((solve(cov$cov)%*%W)*W))>exp(-.5*table@scale.factor^2)
	attr(result,'center') = cov$center
	attr(result,'cov')    = cov$cov
	attr(result,'radius') = table@scale.factor
	result
})
setMethod("summarizeFilter",signature("filterResult","norm2Filter"),function(result,filter) {
	ret = callNextMethod()
	ret$cov    = attr(result@subSet,'cov')
	ret$center = attr(result@subSet,'center')
	ret$radius = attr(result@subSet,'radius')
	ret
})
