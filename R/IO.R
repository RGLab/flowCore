## Adapted from F.Hahne 
## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory

## ==========================================================================
## Reading FCS file header and TEXT section only
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
read.FCSheader <- function(files, path=".", keyword=NULL)
{
  stopifnot(is.character(files), length(files)>=1, files!="")

  filenames <- files
  
  if(path != ".")
    files = file.path(path, files)
   
  res <- lapply(files, header)

  if (!is.null(keyword))
    res <- lapply(res, function(x) x[keyword])

  names(res) <- filenames
  res
 
}


header <- function(files){
    con <- file(files, open="rb")
    offsets <- readFCSheader(con)
    txt     <- readFCStext(con, offsets, debug)
    close(con)
    txt
}

## ==========================================================================
## main wrapper to read FCS files
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
read.FCS <- function(filename,
                     transformation="linearize",
                     nlines= NULL,
                     sampling=FALSE,
                     debug=FALSE,
                     alter.names=FALSE,
                     column.pattern=NULL) {
    
    stopifnot(is.character(filename), length(filename)==1, filename!="")
    con <- file(filename, open="rb")
    if(is.logical(transformation) && transformation || !is.null(transformation) && transformation == "linearize") {
        transformation <- TRUE
        scale <- FALSE
    } else if ( !is.null(transformation) && transformation == "scale") {
        transformation <- FALSE
        scale <- TRUE
    }
    else if (is.null(transformation) || is.logical(transformation) && !transformation) {
        transformation <- FALSE 
        scale <- FALSE
    } 
    offsets <- readFCSheader(con)
    txt     <- readFCStext(con, offsets, debug)
    mat     <- readFCSdata(con, offsets, txt, transformation, nlines, sampling, debug, scale, alter.names)
    params  <- makeFCSparameters(colnames(mat),txt)
    close(con)
    
    
    if(!is.null(column.pattern)) {
        n <- colnames(mat)
        i <- grep(column.pattern,n)
        mat <- mat[,i]
        params <- params[i,]
    }
    
    if(is.null(nlines)){
        if(as.integer(readFCSgetPar(txt, "$TOT"))!=nrow(mat))
          stop(paste("file", filename, "seems to corrupted."))
    }
    
    txt[["FILENAME"]] = filename
    if(transformation==TRUE){
        txt[["transformation"]] <-"applied" 
    }
    
    
    description <- strsplit(txt,split=" ")
    names(description) <- names(txt)
    return(new("flowFrame", exprs=mat, description= description, parameters=params))
}


## ==========================================================================
## create AnnotatedDataFrame describing the flow parameters (channels)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
makeFCSparameters <- function(cn,txt) {
	npar = length(cn)
	id   = paste("$P",1:npar,sep="")
	new("AnnotatedDataFrame",
		data=data.frame(row.names=I(id),name=I(cn),desc=I(txt[paste(id,"S",sep="")]),
                  range=as.numeric(txt[paste(id,"R",sep="")])),
		varMetadata=data.frame(row.names=I(c("name","desc","range")),
                  labelDescription=I(c("Name of Parameter","Description of Parameter","Range of Parameter"))))
}


## ==========================================================================
## match FCS parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSgetPar <- function(x, pnam) {
  stopifnot(is.character(x), is.character(pnam)) 
  i <- match(pnam, names(x))
  if(any(is.na(i)))
    stop(paste("Parameter(s)", pnam, "not contained in 'x'"))
  return(x[i])
}


## ==========================================================================
## parse FCS file header
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSheader <- function(con) {
    
    seek(con, 0)
    version <- readChar(con, 6)
    if(!version %in% c("FCS2.0", "FCS3.0"))
      stop("This does not seem to be a valid FCS2.0 or FCS3.0 file")
    
    version <-  substring(version, 4, nchar(version))
    tmp <- readChar(con, 4)
    stopifnot(tmp=="    ")
    
    coffs <- character(6)
    for(i in 1:length(coffs))
      coffs[i] <- readChar(con=con, nchars=8)
    
    ioffs <- c(as.double(version), as.integer(coffs))
    names(ioffs) <- c("FCSversion", "textstart", "textend", "datastart", "dataend", "anastart", "anaend")
    
    if(all(is.na(ioffs[2:5]) || ioffs[2:5]==""))
      stop("Missing header information to start parsing the binary section of the file")
    return(ioffs)
}


## ==========================================================================
## parse FCS file text section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCStext <- function(con, offsets, debug) {
  seek(con, offsets["textstart"])
  txt <- readChar(con, offsets["textend"]-offsets["textstart"]+1)
  txt <- iconv(txt, "", "latin1", sub="byte")
  delimiter <- substr(txt, 1, 1)
  sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter, fixed=TRUE)[[1]]
  
  rv <- c(offsets["FCSversion"], sp[seq(2, length(sp), by=2)])
  names(rv) <- c("FCSversion", sp[seq(1, length(sp)-1, by=2)])
  return(rv)
}


## ==========================================================================
## read FCS file data section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSdata <- function(con, offsets, x, transformation, nlines, sampling, debug, scale, alter.names) {

    endian <- switch(readFCSgetPar(x, "$BYTEORD"),
                     "4,3,2,1" = "big",
                     "2,1" = "big",
                     "1,2" = "little",     
                     "1,2,3,4" = "little",
                     stop(paste("Don't know how to deal with $BYTEORD", readFCSgetPar(x, "$BYTEORD"))))
    
    
    dattype <- switch(readFCSgetPar(x, "$DATATYPE"),
                      "I" = "integer",
                      "F" = "numeric",
                      stop(paste("Don't know how to deal with $DATATYPE", readFCSgetPar(x, "$DATATYPE")))) 
    
    if (readFCSgetPar(x, "$MODE") != "L")
      stop(paste("Don't know how to deal with $MODE", readFCSgetPar(x, "$MODE")))
    
    nrpar    <- as.integer(readFCSgetPar(x, "$PAR"))
    nrowTotal <- as.integer(readFCSgetPar(x, "$TOT"))
    range    <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep="")))
    bitwidth <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
    bitwidth <- unique(bitwidth)
    if(length(bitwidth)!=1)
      stop("Sorry, I am expecting the bitwidth to be the same for all parameters")
    
    ##for DATA segment exceeding 99,999,999 byte.
    if(offsets["FCSversion"] == 3){
        datastart <- as.numeric(readFCSgetPar(x, "$BEGINDATA"))
        dataend <- as.numeric(readFCSgetPar(x, "$ENDDATA"))
        if(offsets["datastart"] != datastart && offsets["datastart"]== 0){
            offsets["datastart"] <-  datastart
        }
        if(offsets["datastart"] != datastart && offsets["datastart"]!= 0){
            print(datastart)
            print(offsets["datastart"])
            stop("The HEADER and the TEXT segment define different starting point to read the data.")
        }
        
        if(offsets["dataend"] != dataend && offsets["dataend"]== 0){
            offsets["dataend"] <-  dataend
        }
        if(offsets["dataend"] != dataend && offsets["dataend"]!= 0){
            stop("The HEADER and the TEXT segment define different ending point to read the data.")
        }
    }
    
    
    size <- bitwidth/8
    if (!size %in% c(1, 2, 4, 8))
       stop(paste("Don't know how to deal with bitwidth", bitwidth))

    ##Read all reports
    if(is.null(nlines)){
        dat <- c()
        seek(con, offsets["datastart"])
        dat <- readBin(con, dattype, n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                   size=size, signed=FALSE, endian=endian)
    }
    ##Read n lines with or without sampling 
    if(length(nlines) == 1){
        if(sampling == FALSE){
            dat <- c()
            offsets["dataend"] <-  offsets["datastart"] + nrpar * size * nlines
            seek(con, offsets["datastart"])
            dat <- readBin(con, dattype, n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                           size=size, signed=FALSE, endian=endian)
        }
        if(sampling == TRUE){
             dat <- c()
             idx <- sort(sample(1: nrowTotal, size=nlines))        
             for (i in 1:length(idx)){
                 startP <- offsets["datastart"] + (idx[i]-1) * nrpar * size
                 endP   <-  startP + nrpar * size
                 seek(con, startP)
                 temp <- readBin(con, dattype, n = (endP - startP+1)/size,
                   size=size, signed=FALSE, endian=endian) 
                 dat <- c(dat, temp)                 
             }
         }
    }
    ##Read a vector of lines without sampling
    if(length(nlines) > 1){
        if(sampling == FALSE){
            dat <- c()
            for (i in 1:length(nlines)){
                startP <- offsets["datastart"] + (nlines[i]-1) * nrpar * size
                endP   <- startP + nrpar * size
                seek(con, startP)
                temp <- readBin(con, dattype, n = (endP - startP + 1)/size,
                                size=size, signed=FALSE, endian=endian)
                dat <- c(dat, temp)
            }
        }
        else{ stop("Error")}
    }
        
    stopifnot(length(dat)%%nrpar==0)
    dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
    cn  <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
    colnames(dat) <- if(alter.names) structure(make.names(cn),names=names(cn)) else cn
    
    if(transformation) {
        ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""))
        ampli <- do.call("rbind",lapply(ampliPar,function(x) as.integer(unlist(strsplit(x,",")))))
        for (i in 1:nrpar){
            if(ampli[i,1] > 0){
                dat[,i] <- 10^((dat[,i]/(range[i]-1))*ampli[i,1])
            }
        }
    }
    if(scale){
        if(transformation) {
            ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""))
            ampli <- do.call("rbind",lapply(ampliPar,function(x) as.integer(unlist(strsplit(x,",")))))		
            for(i in 1:nrpar) {
                dat[,i] = if(ampli[i,1] > 0) dat[,i]/(10^ampli[i,1]) else dat[,i]/(range[i]-1)
            }
        }
        else {
            for(i in 1:nrpar) {
                dat[,i] = dat[,i]/(range[i]-1)
            }
        }
    }
    return(dat) 
}



## ==========================================================================
## read in flowSets, based on the prada function readCytoSet. Most of the actual
## flowSet construction has been moved to the coercion of lists to flowSets. This
## is to facilitate the construction of flowSets from within the R environment (e.g.
## the flowFrames are actually retrieved from the network or a database.)
## ---------------------------------------------------------------------------
read.flowSet = function(files=NULL, path=".", pattern=NULL, phenoData, descriptions, name.keyword,
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
			warning("Supplied file names will be ignored, ",
                                "using phenoData names instead.")
		file.names = sampleNames(phenoFrame)
		files      = dir(path,paste(gsub("\\.","\\\\\\.",file.names),collapse="|"),full.names=TRUE)
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
		file.names = basename(files) ## strip path from names
                if(path != ".")
                  files = file.path(path, files)
               
	}

         
        flowSet = lapply(files, read.FCS,...)
	#Allows us to specify a particular keyword to use as our sampleNames
	#rather than requiring the filename be used. This is handy when something
	#like SAMPLE ID is a more reasonable choice. Sadly reading the flowSet is
	#a lot more insane now.
        
       
        if(!missing(name.keyword))
		names(flowSet) = sapply(flowSet,keyword,name.keyword)
	else
		names(flowSet) = file.names
      
	flowSet = as(flowSet,"flowSet")
      
	if(!is.null(phenoFrame)) { phenoData(flowSet) = phenoFrame } else if(!missing(phenoData)) {
		#Collect the names for each field in the data frame
		field.names = names(phenoData)
                print(field.names)
		if(is.null(field.names))
			stop("phenoData list must have names")
		field.names = sapply(seq(along=phenoData),function(i) {
			if(length(field.names[i]) == 0) as(phenoData[i],"character")
                        else field.names[i]
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
			varMetadata=data.frame(labelDescriptions=I(descriptions),
                          row.names=field.names))
	}
	flowSet
}

