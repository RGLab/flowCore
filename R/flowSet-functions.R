## ==========================================================================
## read in flowSets, based on the prada function readCytoSet. Most of the actual
## flowSet construction has been moved to the coercion of lists to flowSets. This
## is to facilitate the construction of flowSets from within the R environment (e.g.
## the flowFrames are actually retrieved from the network or a database.)
## ---------------------------------------------------------------------------
read.flowSet = function(files=NULL,path=".",pattern=NULL,phenoData,descriptions,
  sep="\t",...) {
	#A frame of phenoData information
	phenoFrame = NULL
	if(!missing(phenoData)) {
		if(is.character(phenoData) && length(phenoData) == 1)
			phenoFrame = read.phenoData(file.path(path,phenoData),
                          header=TRUE,as.is=TRUE,sep=sep,...)
		else if(is(phenoData,"AnnotatedDataFrame"))
			phenoFrame = phenoData
	}
	
	if(!is.null(phenoFrame)) {
		if(!is.null(files))
			warning("Supplied file names will be ignored, using phenoData names instead.")
		file.names = sampleNames(phenoFrame)
		files      = dir(files,path,full.names=TRUE)
		if(length(files) == 0) 
			stop(paste("Files given by phenoData not found in",path))
	}
	
	if(is.null(files)) {
		files = dir(path,pattern,full.names=TRUE)
		file.names = dir(path,pattern,full.names=FALSE)
		if(length(files)<1)
			stop(paste("No matching files found in ",path))
	} else {
		if(!is.character(files))
			stop("'files' must be a character vector.")
		file.names = files
	}
	#Isn't reading a flowSet easy?
	flowSet = if(!is.null(phenoFrame))
		(as(structure(lapply(files,read.FCS,...),names=file.names),"flowSet") <- phenoFrame)
	else
		as(structure(lapply(files,read.FCS,...),names=file.names),"flowSet")
		
	#If we have data, add it. Otherwise try to use the phenoData list to do something reasonable.
	if(is.null(phenoFrame) && !missing(phenoData)) {
		#Collect the names for each field in the data frame
		field.names = names(phenoData)
		if(is.null(field.names))
			stop("phenoData list must have names")
		field.names = sapply(seq(along=phenoData),function(i) {
			if(length(field.names[i]) == 0) as(phenoData[i],"character") else field.names[i]
		})
		if(!missing(descriptions)) {
			#If the descriptions have names, reorder them as needed.
			if(!is.null(names(descriptions)))
				descriptions = descriptions[field.names]
		} else
			descriptions = field.names
		names(phenoData) = field.names
		phenoData(flowSet) = new("AnnotatedDataFrame",
			data=keyword(flowSet,phenoData),
			varMetadata=data.frame(labelDescriptions=I(descriptions),row.names=field.names))
	}

	flowSet
}
