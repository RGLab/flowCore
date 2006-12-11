setMethod("spillover","flowFrame",function(x) {
	# Flow frames have a SPILL keyword that is often
	# used to store the spillover matrix. Attempt to
	# extract it.
	"spillover"
})
setMethod("spillover","flowSet",function(x,unstained=NULL,patt=NULL,fsc="FSC-A",ssc="SSC-A",method="median") {
	if(is.null(unstained)) {
		"is.null(unstained)"
	} else {
		#We often only want spillover for a subset of the columns 
		cols = if(is.null(patt)) colnames(x) else grep(patt,colnames(x),value=TRUE)
		#Ignore the forward and sidescatter channels if they managed to get included
		if(!is.na(match(fsc,cols))) cols = cols[-match(fsc,cols)] 
		if(!is.na(match(ssc,cols))) cols = cols[-match(ssc,cols)]
		#There has got to be a better way of doing this...
		stains = phenoData(x)$name
		stains = stains[-match(unstained,stains)]
		
		#Grab the baseline from the unstained values
		baseline = apply(exprs(x[[unstained]])[x[[unstained]]%in%norm2Filter(fsc,ssc,scale.factor=1.5),cols],2,method)
		
		#Now do the same thing to all the stains and sweep out the baseline to figure out the motion on all
		#of the channels.
		inten = sweep(sapply(stains,function(s)
			apply(exprs(x[[s]])[x[[s]]%in%norm2Filter(fsc,ssc,scale.factor=1.5),cols],2,method)),2,baseline)
		#We assume that the highest intensity channel in each column is the signal channel. If something weird
		#happens here, you probably screwed up your comp controls (or you're using an awful channel like PacO
		#in which case the mean is probably the recommended statistic).
		colnames(inten) = cols[apply(inten,2,which.max)]
		#Now normalize row-wise to figure out the % spillover and ensure that any negative values are set to
		#0 since negative compensation is even more insane than negative fluoresence.
		inten = apply(inten,2,function(x) ifelse(x<0,0,x/max(x)))		
		t(inten[,match(cols,colnames(inten))])
	}
})


setMethod("compensate",signature("flowFrame","matrix"),function(x,spillover) {
	#Make sure we're normalized to [0,1] and then invert
	cols = colnames(spillover)
	e    = exprs(x)
	e[,cols] = e[,cols]%*%solve(spillover/max(spillover))
	exprs(x) = e
	x
})
setMethod("compensate",signature("flowSet","matrix"),function(x,spillover) {
	y = as(structure(lapply(seq(along=x),function(i) compensate(x[[i]],spillover)),names=phenoData(x)$name),"flowSet")
	phenoData(y) = phenoData(x)
	y
})