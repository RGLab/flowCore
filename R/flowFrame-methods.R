setMethod("identifier", signature="flowFrame",
          definition=function(object){
            oid <- object@description["GUID"]
            if(is.null(oid) || is.na(oid))
              as.vector(object@description["$FIL"])
            else
              as.vector(oid)
          })


## ==========================================================================
## accessor method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("exprs", signature="flowFrame",
          definition=function(object)
            object@exprs,
          )
## ==========================================================================


## ==========================================================================
## replace method for slot exprs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("exprs", signature=c("flowFrame", "matrix"),
                 definition=function(object, value) {
                   object@exprs <- value
                   return(object)
                 })
## ==========================================================================


## ==========================================================================
## accessor method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("description", signature="flowFrame",
          definition=function(object)
            object@description,
          )
## ==========================================================================


## ==========================================================================
## replace method for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("description", signature=c("flowFrame", "character"),
                 definition=function(object, value) {
                   object@description <- value
                   return(object)
                 })
## ==========================================================================


## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames", signature="flowFrame",
          definition=function(x, do.NULL="missing", prefix="missing")
          colnames(exprs(x))
          )

setMethod("featureNames", signature="flowFrame",
          definition=function(object)
          object@parameters$desc
          )

setMethod("names", signature="flowFrame",
          definition=function(x) {
            cn <- colnames(x)
            fn <- featureNames(x)
            if(length(fn) == length(cn)) {
              cn <- paste("<", cn, ">", sep="")
              for(i in seq(along=fn)) {
                if(!is.na(fn[i]) && fn[i]!="")
                  cn[i] <- paste(cn[i],fn[i])
              }
            }
            cn
          })

## ==========================================================================


## ==========================================================================
## replace method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("colnames", signature=c("flowFrame", "ANY"),
                 definition=function(x, value) {
                   colnames(x@exprs) <- value
                   return(x)
                 })
## ==========================================================================


## ==========================================================================
## a simple plot method without strange plot parameter and friends. It does
## the most intuitive thing: take a flowFrame and plot the first two columns
## of the data. If you want to plot other columns, subset the frame before
## plotting. For the more complicated stuff there may well be more elaborate
## methods but I advice to separate functionality. Let R do the method
## dispatching instead of putting everything into one "do-it-all" method.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot", signature(x="flowFrame", y="missing"),
          definition=function(x, ...){
            values=exprs(x)
            geneplotter:::smoothScatter(values[,1:2], ...)
          })
setMethod("plot",signature(x="flowFrame",y="character"),function(x,y,...) {
	l = length(y)
	if(l==1)
		hist(exprs(x)[,y],...)
	else if(l==2) {
		geneplotter:::smoothScatter(exprs(x)[,y],...)
	}
})
## ==========================================================================


## ==========================================================================
## the $-operator
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"$.flowFrame" <- function(x, val)
    x[,val]
## ==========================================================================


## ==========================================================================
## subsetting method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("[", signature="flowFrame",
          definition=function(x, i, j, ..., drop=FALSE) {
            exprs(x) <-  switch(1+missing(i)+2*missing(j),
                                { exprs(x)[i, j, ..., drop=drop] },
                                { exprs(x)[ , j, ..., drop=drop] },
                                { exprs(x)[i,  , ..., drop=drop] },
                                { exprs(x)[ ,  , ..., drop=drop] } )
            if(!missing(j)) 
              x@parameters = x@parameters[j,]
            x
          })
setMethod("[", signature=signature("flowFrame","filterResult"),
          definition=function(x,i,j,...,drop=FALSE) {
            if(missing(j))
              x[x %in% i,,...,drop=drop]
            else
              x[x %in% i,j,...,drop=drop]
          })
setMethod("[", signature=signature("flowFrame","filter"),
          definition=function(x,i,j,...,drop=FALSE) {
			result = filter(x,i)
            if(missing(j))
              x[result,,...,drop=drop]
            else
              x[result,j,...,drop=drop]
          })
## ==========================================================================


## ==========================================================================
## nrow method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("nrow", signature=signature("flowFrame"),
          definition=function(x)
            return(nrow(x@exprs))
          )
## ==========================================================================


## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ncol", signature=signature("flowFrame"),
          definition=function(x)
            return(ncol(x@exprs))
          )
## ==========================================================================


## ==========================================================================
## show method for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature=signature("flowFrame"),
          definition=function(object) {
            dm <- dim(exprs(object))
            msg <- paste("flowFrame object with ", dm[1], " cells and ", 
                         dm[2], " observables:\n", paste(names(object), 
                         collapse = " "), "\nslot 'description' has ",
                         length(description(object)), " elements\n", sep = "")
            cat(msg)
            return(msg)
          })
## ==========================================================================


## ==========================================================================
## Subset method for flowFrame: Why is this Subset with capital 's' ???
## Why do we need the 'select' parameter? Wouldn't this be equivalent:
## Subset(x[,c(1,3)], subset)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("Subset", signature("flowFrame","filter"),
          definition=function(x,subset,select,...){
            if(!missing(select))
              Subset(x,x %in% subset,select,...)
            else
              Subset(x,x %in% subset,...)
          })
setMethod("Subset",signature("flowFrame","logical"),
          definition=function(x,subset,select,...) {
            if(!missing(select))
              x[subset & !is.na(subset),select]
            else x[subset,]
          })


## ==========================================================================
## split method for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("split", signature("flowFrame","filter"),
          definition=function(x,f,drop=FALSE,...)
          split(x,filter(x,f),drop,...)
          )

#We actually filter on filterResults and multipleFilterResults
setMethod("split", signature("flowFrame","logicalFilterResult"),
          definition=function(x,f,drop=FALSE,population=NULL,...) {
            if(is.null(population)) population=f@filterId
            structure(list(x[f@subSet,],x[!f@subSet,]),
                      names=c(paste(population,"+",sep=""),
                        paste(population,"-",sep="")))
          })

setMethod("split", signature("flowFrame","multipleFilterResult"),function(x,f,drop=FALSE,prefix=NULL,...) {
	nn = if(is.null(prefix)) names(f) else paste(prefix,names(f),sep="")
	structure(lapply(seq(along=f),function(i) x[f[[i]],]),names=nn)
})


setMethod("summary", signature("flowFrame"), 
    function(object, ...) 
        apply(exprs(object), 2, summary)
 )
