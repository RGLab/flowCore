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
setMethod("plot",
  signature(x="flowFrame", y="missing"),
  definition=function(x, col=densCols(exprs(x)[,1:2]), pch=20, ...){
    values=exprs(x)
    plot(values, col=col, pch=pch, ...)})
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
    x
  },
  valueClass="flowFrame")
setMethod("[",signature=signature("flowFrame","filterResult"),definition=function(x,i,j,...,drop=FALSE) {
	if(missing(j))
		x[i@subSet==1,,...,drop=drop]
	else
		x[i@subSet==1,j,...,drop=drop]
},valueClass="flowFrame")
## ==========================================================================



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
        dm[2], " observables:\n", paste(colnames(exprs(object)), 
            collapse = " "), "\nslot 'description' has ", length(description(object)), 
        " elements\n", sep = "")
    cat(msg)
    return(msg)
})
## ==========================================================================
