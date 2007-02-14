## ==========================================================================
## %in% method
## ---------------------------------------------------------------------------
setMethod("%in%",signature("ANY","multipleFilterResult"),function(x,table)
          stop("multipleFilterResult: You must specify a subpopulation."))


## ==========================================================================
## summary method
## ---------------------------------------------------------------------------
setMethod("summary","multipleFilterResult",function(object,...) {
	true = summary(object@subSet[!is.na(object@subSet)])
	count= rep(sum(!is.na(object@subSet)),nlevels(object@subSet))
	structure(list(name=levels(object@subSet),true=true,false=count-true,
                       n=count,p=true/count,q=1-(true/count)),class="filterSummary")
})


## ==========================================================================
## length method
## ---------------------------------------------------------------------------
setMethod("length","multipleFilterResult",function(x) nlevels(x@subSet))


## ==========================================================================
## names method
## ---------------------------------------------------------------------------
setMethod("names","multipleFilterResult",function(x) levels(x@subSet))


## ==========================================================================
## subsetting method
## ---------------------------------------------------------------------------
setMethod("[[",signature("multipleFilterResult"),function(x,i,j,drop=FALSE) {
	if(is.numeric(i)) i = names(x)[i]
	if(length(i)!=1)
		stop("Only a single subpopulation can be selected.")
	filterDetails = x@filterDetails
	filterDetails$population = i
	new("logicalFilterResult",parameters=x@parameters,subSet=(x@subSet==i &
                                    !is.na(x@subSet)),filterDetails=filterDetails,
		frameId=x@frameId,filterId=paste(x@filterId,"[[",i,"]]",sep=""))
})
