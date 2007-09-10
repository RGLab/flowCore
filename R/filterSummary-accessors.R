setMethod("[[",signature("filterSummary","numeric"),function(x,i,j,...) {
	list(name=x@name[i],true=x@true[i],false=x@count[i]-x@true[i],count=x@count[i],p=x@p[i],q=1-x@p[i])
})
setMethod("[[",signature("filterSummary","character"),function(x,i,j,...) {
	i = which(i,names(x))
	x[[i]]
})
setMethod("length",signature("filterSummary"),function(x) length(x@name))
setMethod("names",signature("filterSummary"),function(x) x@name)

setMethod("$",signature("filterSummary","ANY"),function(x,name) "$"(as(x,"data.frame"),name))
