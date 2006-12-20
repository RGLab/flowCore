setMethod("identifier","flowFrame",function(object) if(is.null(object@description["GUID"])) object@description["$FIL"] else object@description["GUID"])


## ==========================================================================
## accessor method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("exprs",
  signature="flowFrame", definition=function(object) object@exprs,
  valueClass="matrix")
## ==========================================================================


## ==========================================================================
## replace method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("exprs",
  signature=c("flowFrame", "matrix"), definition=function(object, value) {
    object@exprs <- value
    return(object)})
## ==========================================================================


## ==========================================================================
## accessor method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("description",
  signature="flowFrame", definition=function(object) object@description,
  valueClass="character")
## ==========================================================================


## ==========================================================================
## replace method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("description",
  signature=c("flowFrame", "character"), definition=function(object, value) {
    object@description <- value
    return(object)})
## ==========================================================================


## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
  signature="flowFrame", definition=function(x, do.NULL="missing",
  prefix="missing") colnames(exprs(x)), valueClass="character")
setMethod("featureNames","flowFrame",function(object) object@parameters$desc)
setMethod("names","flowFrame",function(x) {
	cn = colnames(x)
	fn = featureNames(x)
	if(length(fn) == length(cn)) {
		for(i in seq(along=fn)) {
			if(!is.na(fn[i])) cn[i] = paste("<",cn[i],"> ",fn[i],sep="")
		}
	}
	cn
})

## ==========================================================================


## ==========================================================================
## replace method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("colnames",
  signature=c("flowFrame", "ANY"), definition=function(x, value) {
    colnames(x@exprs) <- value
    return(x)})
## ==========================================================================


## ==========================================================================
## plot method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#setMethod("plot",
#  signature(x="flowFrame", y="missing"),
#  definition=function(x, col=densCols(exprs(x)[,1:2]), pch=20, ...){
#    values=exprs(x)
#    browser()
#    plot(values, col=col, pch=pch, ...)})
## ==========================================================================


## ==========================================================================
## the $-operator
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"$.flowFrame" <- function(x, val)
    (description(x))[val]
## ==========================================================================


## ==========================================================================
## subsetting method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[",
  signature="flowFrame", definition=function(x, i, j, ..., drop=FALSE) {
    exprs(x) <-  switch(1+missing(i)+2*missing(j),
         { exprs(x)[i, j, ..., drop=drop] },
         { exprs(x)[ , j, ..., drop=drop] },
         { exprs(x)[i,  , ..., drop=drop] },
         { exprs(x)[ ,  , ..., drop=drop] } )
	if(!missing(j)) 
		x@parameters = x@parameters[j,]
    x
  },
  valueClass="flowFrame")
setMethod("[",signature=signature("flowFrame","filterResult"),definition=function(x,i,j,...,drop=FALSE) {
	if(missing(j))
		x[as(i,"logical"),,...,drop=drop]
	else
		x[as(i,"logical"),j,...,drop=drop]
},valueClass="flowFrame")
## ==========================================================================
setMethod("[",signature=signature("flowFrame","filterResult"),definition=function(x,i,j,...,drop=FALSE) {
	if(missing(j))
		x[as(i,"logical"),,...,drop=drop]
	else
		x[as(i,"logical"),j,...,drop=drop]
},valueClass="flowFrame")
setMethod("[",signature=signature("flowFrame","filterResult"),definition=function(x,i,j,...,drop=FALSE) {
	if(missing(j))
		x[x%in%i,,...,drop=drop]
	else
		x[x%in%i,j,...,drop=drop]
},valueClass="flowFrame")

## ==========================================================================
## nrow method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("nrow",
  signature=signature("flowFrame"),
    definition=function(x) {
    return(nrow(x@exprs))})
## ==========================================================================


## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ncol",
  signature=signature("flowFrame"),
    definition=function(x) {
    return(ncol(x@exprs))})
## ==========================================================================


## ==========================================================================
## show method for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature=signature("flowFrame"),definition=function(object) {
   dm <- dim(exprs(object))
    msg <- paste("flowFrame object with ", dm[1], " cells and ", 
        dm[2], " observables:\n", paste(names(object), 
            collapse = " "), "\nslot 'description' has ", length(description(object)), 
        " elements\n", sep = "")
    cat(msg)
    return(msg)
})
## ==========================================================================

setMethod("Subset",signature("flowFrame","filter"),function(x,subset,select,...)
	if(!missing(select)) Subset(x,x %in% subset,select,...) else Subset(x,x %in% subset,...))
setMethod("Subset",signature("flowFrame","logical"),function(x,subset,select,...) 
	if(!missing(select)) x[subset & !is.na(subset),select] else x[subset,])

setMethod("split",signature("flowFrame","filter"),function(x,f,drop=FALSE,...) split(x,filter(x,f),drop,...))

#We actually filter on filterResults and multipleFilterResults
setMethod("split",signature("flowFrame","filterResult"),function(x,f,drop=FALSE,population=NULL,...) {
	if(is.null(population)) population=f@filterId
	structure(list(x[f@subSet,],x[!f@subSet,]),names=c(paste(population,"+",sep=""),paste(population,"-",sep="")))
})
setMethod("split",signature("flowFrame","multipleFilterResult"),function(x,f,drop=FALSE,prefix=NULL,...) {
	structure(lapply(seq(along=f),function(i) x[f[[i]],]),names=if(!is.null(prefix)) paste(prefix,names(f),sep="") else names(f))
})
