#### Define population inclusion methods 

setMethod("%in%",signature(x="flowFrame",table="rectangleGate"),function(x,table) {
	e = exprs(x)
	apply(sapply(seq(along=table@parameters),function(i) {
		!is.na(cut(e[,table@parameters[i]],c(table@min[i],table@max[i]),labels=FALSE))
	}),1,all)
})


