## last changed November 2, 2006
## pdh: rectangleGate was still using cytoFrame. I changed it to flowFrame
##      made the filterDetails more uniform among the methods
##      got the right filename for a flowFrame to go in the msg

## ==========================================================================
## Apply polygon filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("applyFilter",
          signature=signature(filter="polygonGate",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-polygonFiltering(filter,flowObject,parent)
              msg <- paste("polygonGate applied on ",
                           deparse(substitute(flowObject)),
                           " (file:",basename(flowObject@description["$FIL"]),
                           ") a ", class(flowObject)," object", sep="")
       			out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          } 
          )
## ==========================================================================



  

## ==========================================================================
##  Apply rectangular filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature(x="flowFrame",table="rectangleGate"),function(x,table) {
	e = if(length(table@parameters)==1) as.matrix(exprs(x)[,table@parameters]) else exprs(x)[,table@parameters]
	apply(sapply(seq(along=table@parameters),function(i) {
		!is.na(cut(e[,i],c(table@min[i],table@max[i]),labels=FALSE))
	}),1,all)
})
setMethod("show","rectangleGate",function(object) {
	cat("Rectangular gate with dimensions:\n")
	for(i in seq(along=object@parameters)) {
		cat("\t")
		cat(object@parameters[i])
		cat(": (")
		cat(paste(object@min[i],object@max[i],sep=","))
		cat(")\n")
	}
})
setMethod("applyFilter",
          signature=signature(filter="rectangleGate",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-rectangleFiltering(filter,flowObject,parent)
               msg <- paste("rectangleGate applied on ",
                            deparse(substitute(flowObject)),
                            " (file:",basename(flowObject@description["$FIL"]),
                            ") a ", class(flowObject)," object", sep="")
               out = list(msg,filter=filter)
              new("filterResult", subSet=selectNew,filterDetails=out)
          })
## ==========================================================================


setMethod("%in%",signature("flowFrame",table="modeGate"),function(x,table) {
	params       = x@parameters
	pp           = structure(rbind(rep(x@bw,length=length(params)),rep(x@n,length=length(params))),dimnames=list(NULL,params))
	x %in% rectGate(apply(pp,2,function(y) find.mode.bounds(x[[names(y)]])))
})

## ==========================================================================
##  Apply norm2Filter on flowFrame Object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature("flowFrame",table="norm2Filter"),function(x,table) {
	if(length(table@parameters) != 2)
		stop("norm2 filters require exactly two parameters.")
	y = exprs(x)[,table@parameters]
	if(length(table@transformation)>0) {
		for(i in seq(along=table@parameters)) {
			f = table@transformation[[table@parameters[i]]]
			if(is.function(f)) y[,i] = f(y[,i])
		}
	}
	if(is.na(match(table@method,c("covMcd","cov.rob"))))
		stop("Method must be either 'covMcd' or 'cov.rob'")
	cov = switch(table@method,
		covMcd = {
			if(nrow(y)>50000) covMcd(y[sample(50000,nrow(y)),]) else covMcd(y)
		},
		cov.rob={cov.rob(y)},
		stop("How did you get here?")
		)
	W  = t(y)-cov$center
	exp(-.5*colSums((solve(cov$cov)%*%W)*W))>exp(-.5*table@scale.factor^2)
})

setMethod("applyFilter",
          signature=signature(filter="norm2Filter",flowObject="flowFrame",parent="ANY"),
          definition=function(filter,flowObject,parent) {
              selectNew <-normFiltering(filter,flowObject,parent)
              msg <- paste("norm2Filter applied on ",
                            deparse(substitute(flowObject)),
                            " (file:",basename(flowObject@description["$FIL"]),
                            ") a ", class(flowObject)," object", sep="")
               out = list(msg,filter=filter,mu=selectNew$mu,S=selectNew$S)
              new("filterResult", subSet=selectNew$sel,filterDetails=out)
          })
## ==========================================================================


