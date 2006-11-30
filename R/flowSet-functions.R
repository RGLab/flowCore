## read in flowSets, based on the prada function readCytoSet. Most of the actual
## flowSet construction has been moved to the coercion of lists to flowSets. This
## is to facilitate the construction of flowSets from within the R environment (e.g.
## the flowFrames are actually retrieved from the network or a database.)
read.flowSet = function(files=NULL,path=".",pattern=NULL,phenoData,sep="\t",...) {
	if(!missing(phenoData)) {
		if(is.character(phenoData))
			phenoData = read.phenoData(file.path(path,phenoData),header=TRUE,as.is=TRUE,sep=sep,...)
		if(!is(phenoData,"AnnotatedDataFrame"))
			stop("'phenoData' must be of type 'AnnotatedDataFrame'.")
		if(!("name" %in% colnames(pData(phenoData))))
			stop("'phenoData' must contain a column 'name'")
		if(!is.null(files))
			warning("'files' ignored. Using 'name' column from 'phenoData' instead.")
		files = phenoData$name

	}
	if(is.null(files)) {
		files = dir(path,pattern)
		if(length(files)<1)
			stop(paste("No matching files found in ",path))
	}
	if(!is.character(files))
		stop("'files' must be a character vector.")
	flowSet = lapply(files,read.FCS,...)
	names(flowSet) = files
	flowSet = as(flowSet,"flowSet") #Convert from a list to an actual flowSet.
	if(!missing(phenoData)) phenoData(flowSet) = phenoData
	flowSet
}