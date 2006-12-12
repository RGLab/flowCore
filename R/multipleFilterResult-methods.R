setMethod("%in%",signature("ANY","multipleFilterResult"),function(x,table) stop("multipleFilterResult: You must specify a subpopulation."))
setMethod("[[",signature("multipleFilterResult"),function(x,i,j,drop=FALSE) {
	i = if(!is.numeric(i)) match(i,x@populations) else i
	if(length(i)!=1)
		stop("Only a single subpopulation can be selected.")
	new("filterResult",subSet=x@subSet==i,filterDetails=x@filterDetails)
})
setMethod("summary","multipleFilterResult",function(object,...) {
	true  = tabulate(object@subSet[!is.na(object@subSet)])
	#If there is a "0" class then remove the first class
	if(min(object@subSet)==0)
		true = true[-1]
	#Add the names for niceness
	names(true) = object@populations
	count = length(x@subSet)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),name=object@populations),class="filterSummary")	
})