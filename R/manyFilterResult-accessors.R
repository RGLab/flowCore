## ==========================================================================
## %in% method
## ---------------------------------------------------------------------------
setMethod("%in%",signature("ANY","manyFilterResult"),function(x,table)
          stop("manyFilterResult: You must specify a subpopulation."))


## ==========================================================================
## summary method
## ---------------------------------------------------------------------------
setMethod("summary","manyFilterResult",function(object,...) {
	n = colnames(object@subSet)
	structure(sapply(n,function(i) summary(object[[i]])),class="filterSummary")
})


## ==========================================================================
## length method
## ---------------------------------------------------------------------------
setMethod("length","manyFilterResult",function(x) ncol(x@subSet))


## ==========================================================================
## names method
## ---------------------------------------------------------------------------
setMethod("names","manyFilterResult",function(x) colnames(x@subSet))


## ==========================================================================
## subsetting method
## ---------------------------------------------------------------------------
setMethod("[[",signature("manyFilterResult"),function(x,i,j,drop=FALSE) {
	if(is.numeric(i)) i = names(x)[i]
	if(length(i)!=1)
		stop("Only a single subpopulation can be selected.")
	filterDetails = filterDetails(x,i)
	filterDetails$population = i
	new("logicalFilterResult",subSet=x@subSet[,i],filterDetails=filterDetails,
		frameId=x@frameId,filterId=paste(x@filterId,"[[",i,"]]",sep=""))
})
