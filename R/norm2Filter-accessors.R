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
                  ## covMcd will be deprecated, need to use CovMcd which produces S4 output
			tmp <- if(nrow(y)>table@n) CovMcd(y[sample(nrow(y),table@n),]) else CovMcd(y)
                        list(center=tmp@center, cov=tmp@cov)
		},
		cov.rob={cov.rob(y)},
		stop("How did you get here?")
		)
	W  = t(y)-cov$center

        ## this used to be:

        ## result = exp(-.5*colSums((solve(cov$cov)%*%W)*W))>exp(-.5*table@scale.factor^2)

        ## which is bad because (1) it directly inverts cov$cov (which
        ## is probably not too bad for a 2x2 matrix) and (2) uses
        ## exp() unnecessarily (there are many rows, exp() is
        ## expensive, but redundant as it is a monomotone increasing
        ## transformation applied on both sides of the inequality).
        ## -DS (2007/04/30)

        result = colSums((qr.solve(cov$cov) %*% W) * W) < table@scale.factor^2

        ## FIXME: a long term change might be to save chol(cov$cov)
        ## rather than cov$cov in the result.  This helps in computing
        ## the gate boundaries and qr.solve above could be replaced by
        ## the equivalent of chol2inv(chol(cov$cov)).

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
