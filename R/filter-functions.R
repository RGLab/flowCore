## ==========================================================================
## function to print filter summary
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print.filterSummary <- function(x,...) {
    
    if(length(x$name) == 1) {
		with(x,cat(sprintf("%s: %d of %d (%.2f%%)\n",name,true,n,100*p)))
	} else {
	for(i in seq(along=x$name))
			with(x,cat(sprintf("%s: %d of %d (%.2f%%)\n",name[i],true[i],n[i],100*p[i])))
	}
}

# quadGate = function(.split,...,ids=names(.split)) {
# 	if(missing(.split)) .split = list(...)
# 	n = names(.split)
# 	makemat = function(i,pos)
# 		matrix(if(pos) c(.split[i],Inf) else c(-Inf,.split[i]),2,dimnames=list(c("min","max"),n[i]))
# 	
# 	gates = list(rectangleGate(makemat(1,FALSE),filterId=paste(ids[1],"-",sep="")),rectangleGate(makemat(1,TRUE),filterId=paste(ids[1],"+",sep="")))
# 	for(i in 2:length(.split)) {
# 		cur = list(rectangleGate(makemat(i,FALSE),filterId=paste(ids[i],"-",sep="")),rectangleGate(makemat(i,TRUE),filterId=paste(ids[i],"+",sep="")))
# 		gates = outer(gates,cur)
# 	}
# 	gates
# }