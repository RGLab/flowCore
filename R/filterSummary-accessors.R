setMethod("[[",signature("filterSummary","numeric"),function(x,i,j,...) {
	list(name=x@name[i],true=x@true[i],false=x@count[i]-x@true[i],count=x@count[i],p=x@p[i],q=1-x@p[i])
})
setMethod("[[",signature("filterSummary","character"),function(x,i,j,...) {
	i = which(i,names(x))
	x[[i]]
})
setMethod("length",signature("filterSummary"),function(x) length(x@name))
setMethod("names",signature("filterSummary"),function(x) x@name)

setMethod("$",signature("filterSummary","ANY"),function(x,name) {
	switch(name,"n"=x@count,"true"=x@true,"false"=x@count-x@true,"p"=x@p,"q"=1-x@p)
})
setMethod("show",signature("filterSummary"),function(object) {
    if(length(object@name) == 1) {
		cat(sprintf("%s: %d of %d (%.2f%%)\n",object@name,object@true,object@count,object@p*100))
	} else {
		for(i in seq(along=x@name)) {
			cat(sprintf("%s: %d of %d (%.2f%%)\n",object@name[i],object@true[i],object@count[i],object@p[i]*100))
		}
	}
})