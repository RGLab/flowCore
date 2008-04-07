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
## length method, how many populations do we have?
## ---------------------------------------------------------------------------
setMethod("length","multipleFilterResult",function(x) nlevels(x@subSet))



## ==========================================================================
## names method, the levels of the subSet factor
## ---------------------------------------------------------------------------
setMethod("names","multipleFilterResult",function(x) levels(x@subSet))


## ==========================================================================
## subsetting method: we create a new logicalFilterResult if the index is
## of length one, otherwise we subset our multipleFilterResult directly
## ---------------------------------------------------------------------------
setMethod("[[",signature("multipleFilterResult"),function(x,i,j,drop=FALSE) {
	if(is.numeric(i))
            i <- names(x)[i]
        filterDetails <- x@filterDetails
        filterDetails$population <- i
	if(length(i)!=1){
           x@subSet <- factor(x@subSet[x@subSet %in% i & !is.na(x@subSet)])
           x@filterDetails <- filterDetails
           return(x)
       }else{
           new("logicalFilterResult",subSet=(x@subSet==i & !is.na(x@subSet)),
               filterDetails=filterDetails,
               frameId=x@frameId, filterId=paste(x@filterId,"[[",i,"]]",
                                  sep=""))
       }
    })
