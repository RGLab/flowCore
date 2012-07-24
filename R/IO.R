## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory


## ==========================================================================
## Determine which 'files' are valid FCS files
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isFCSfile <- function(files)
{
    sapply(files, function(f){
        if (file.exists(f)) {
            con <- file(f, open="rb")
            on.exit(close(con))
            version <- readChar(con, 6)
            isTRUE(version %in% c("FCS2.0", "FCS3.0"))
        }
        else FALSE
    })
}


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

header <- function(files,emptyValue=TRUE){
    con <- file(files, open="rb")
    offsets <- readFCSheader(con)
    txt     <- readFCStext(con, offsets,emptyValue=emptyValue)
    close(con)
    txt
}


## ==========================================================================
## main wrapper to read FCS files
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
read.FCS <- function(filename,
                     transformation="linearize",
                     which.lines=NULL,
                     alter.names=FALSE,
                     column.pattern=NULL,
                     decades=0,
                     ncdf=FALSE,
                     min.limit=NULL,
                     dataset=NULL
			 		,emptyValue=TRUE)
{
    ## check file name
    if(!is.character(filename) ||  length(filename)!=1)
        stop("'filename' must be character skalar")
    if(!file.exists(filename))
        stop(paste("'", filename, "' is not a valid file", sep=""))
    con <- file(filename, open="rb")
    on.exit(close(con))

    ## transform or scale data?
    if(is.logical(transformation) && transformation ||
       !is.null(transformation) && transformation == "linearize") {
        transformation <- TRUE
        scale <- FALSE
    } else if ( !is.null(transformation) && transformation == "scale") {
        transformation <- TRUE
        scale <- TRUE
    } else if (is.null(transformation) || is.logical(transformation) &&
               !transformation) {
        transformation <- FALSE 
        scale <- FALSE
    } 

    ## read the file  
    offsets <- findOffsets(con,emptyValue=emptyValue)
    ## check for multiple data sets
    if(is.matrix(offsets))
    {
        nd <- nrow(offsets)
        if(is.null(dataset))
        {
            warning(sprintf("The file contains %d additional data segment%s.\n",
                            nd-1, ifelse(nd>2, "s", "")),
                    "The default is to read the first segment only.\nPlease consider ",
                    "setting the 'dataset' argument.", call.=FALSE)
            offsets <- offsets[1,]
        }
        else    
        {
            if(!is.numeric(dataset) || !dataset %in% seq_len(nd))
                stop(sprintf("Argument 'dataset' must be an integer value in [1,%d].",
                             nd))
            offsets <- offsets[dataset,]
        }
    }
    txt <- readFCStext(con, offsets,emptyValue=emptyValue)
    ## We only transform if the data in the FCS file hasn't already been
    ## transformed before
    if("transformation" %in% names(txt) &&
       txt[["transformation"]] %in% c("applied", "custom"))
       transformation <- FALSE
    mat <- readFCSdata(con, offsets, txt, transformation, which.lines,
                       scale, alter.names, decades, min.limit)
    matRanges <- attr(mat,"ranges")

	
    id <- paste("$P",1:ncol(mat),sep="")
    zeroVals <- as.numeric(sapply(strsplit(txt[paste(id,"E",sep="")], ","),
                                  function(x) x[2]))
    absMin <- apply(mat,2,min,na.rm=TRUE)
    realMin <- pmin(zeroVals,pmax(-111, absMin, na.rm=TRUE), na.rm=TRUE)
    
	if("transformation" %in% names(txt) && txt[["transformation"]] == "custom") {
		for(i in seq_along(colnames(mat))) {
			realMin[i] <- txt[[sprintf("flowCore_$P%sRmin", i)]]
		}
	}

    params <- makeFCSparameters(colnames(mat),txt, transformation, scale,
                                 decades, realMin)
    
    ## only keep certain parameters
    if(!is.null(column.pattern)) {
        n <- colnames(mat)
        i <- grep(column.pattern,n)
        mat <- mat[,i]
        params <- params[i,]
    }

    ## check for validity
    if(is.null(which.lines)){
        if(as.integer(readFCSgetPar(txt, "$TOT"))!=nrow(mat))
            stop(paste("file", filename, "seems to be corrupted."))
    }

	## set transformed flag and fix the PnE and the Datatype keywords
    ## also add our own PnR fields.
    txt[["FILENAME"]] <- filename
    if(transformation==TRUE) { 
       txt[["transformation"]] <-"applied"
       for(p in seq_along(pData(params)$name)) {
          txt[[sprintf("$P%sE", p)]] <- sprintf("0,%g", 0) 
          txt[[sprintf("flowCore_$P%sRmax", p)]] <- matRanges[p] +1
          txt[[sprintf("flowCore_$P%sRmin", p)]] <- realMin[p] 
       }
       txt[["$DATATYPE"]] <- "F"
    }
    ## build description from FCS parameters
	
	if(offsets["FCSversion"]<=2)
	{
		description <- strsplit(txt,split="\n")
	    names(description) <- names(txt)
		
	}else
	{
		description <- strsplit(txt, split=NA) # not really splitting, but converting the data structure}	
	}
	
	
	
    ## the spillover matrix
    spID <- intersect(c("SPILL", "spillover"), names(description))
    if(length(spID)>0){
        sp <- description[[spID]]
        splt <- strsplit(sp, ",")[[1]]
        nrCols <- as.numeric(splt[1])
        cnames <- splt[2:(nrCols+1)]
        vals <- as.numeric(splt[(nrCols+2):length(splt)])
        spmat <- matrix(vals, ncol=nrCols, byrow=TRUE)
        colnames(spmat) <- cnames
        description[[spID]] <- spmat
    }
    tmp <- new("flowFrame", exprs=mat, description= description,
               parameters=params)
    identifier(tmp) <- basename(identifier(tmp))
    if(ncdf && nrow(tmp) > 0){
        if(!require(ncdf))
            stop("You need to have package 'ncdf' installed in order for ",
                 "this feature to work.", call.=FALSE)
        ttmp <- ncdfExpressionMatrix(tmp)
        handler <- new("ncdfHandler", file=ttmp$file, pointer=ttmp$pointer, open=TRUE)
        tmp@exprs <- handler
    }
    return(tmp)
}


## ==========================================================================
## create AnnotatedDataFrame describing the flow parameters (channels)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
makeFCSparameters <- function(cn, txt, transformation, scale, decades,
                              realMin) {

    npar <- length(cn)
    id <- paste("$P",1:npar,sep="")
    rid <- paste("flowCore_", id,"Rmax",sep="")
    original <- is.na(txt[rid[1]])
    range <- origRange <- if(!original) as.numeric(txt[rid]) + 1 else
    as.numeric(txt[paste(id,"R",sep="")])
    range <- rbind(realMin,range-1)

    ## make sure the ranges are transformed along with the data
    if(transformation & !scale){
      
        ampliPar <- txt[paste(id,"E",sep="")]
        noPnE <- is.na(ampliPar)
        if(any(noPnE))
            ampliPar[noPnE] <- "0,0"
        ampli <- do.call(rbind,lapply(ampliPar, function(x)
                                        as.integer(unlist(strsplit(x,",")))))
        for (i in 1:npar)
            if(ampli[i,1] > 0)
                range[,i] <- 10^((range[,i]/(origRange[i]-1))*ampli[i,1])
    }
    else if(scale)
        range[2,] <- rep(10^decades, npar)
       
    new("AnnotatedDataFrame",
        data=data.frame(row.names=I(id),name=I(cn),
        desc=I(txt[paste(id,"S",sep="")]),
        range=as.numeric(txt[paste(id,"R",sep="")]), minRange=range[1,], maxRange=range[2,]),
        varMetadata=data.frame(row.names=I(c("name","desc","range",
                               "minRange", "maxRange")),
        labelDescription=I(c("Name of Parameter","Description of Parameter",
        "Range of Parameter", "Minimum Parameter Value after Transforamtion",
        "Maximum Parameter Value after Transformation"))))
}


## ==========================================================================
## match FCS parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSgetPar <- function(x, pnam, strict=TRUE)
{

    stopifnot(is.character(x), is.character(pnam)) 
    i <- match(pnam, names(x))
    if(any(is.na(i)) && strict)
        stop(paste("Parameter(s)", pnam[is.na(i)], "not contained in 'x'\n"))
    if(!strict)
    {
        if(!all(is.na(i)))
            i[!is.na(i)] <- x[i[!is.na(i)]]
        names(i) <- pnam
        return(i)
    }
    return(x[i])
}


## ==========================================================================
## Find all data sections in a file and record their offsets.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
findOffsets <- function(con,emptyValue=TRUE)
{
    offsets <- readFCSheader(con)
    txt <- readFCStext(con, offsets,emptyValue=emptyValue)
    addOff <- 0
    nd <- as.numeric(txt[["$NEXTDATA"]])
    while(nd != 0)
    {
        addOff <- addOff + nd
        offsets <- rbind(offsets, readFCSheader(con, addOff))
        txt <- readFCStext(con, offsets[nrow(offsets),],emptyValue=emptyValue)
        nd <- as.numeric(txt[["$NEXTDATA"]])
    }
    return(offsets)
}


## ==========================================================================
## parse FCS file header
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSheader <- function(con, start=0)
{
    seek(con, start)
    version <- readChar(con, 6)
    if(!version %in% c("FCS2.0", "FCS3.0"))
        stop("This does not seem to be a valid FCS2.0 or FCS3.0 file")
    
    version <-  substring(version, 4, nchar(version))
    tmp <- readChar(con, 4)
    stopifnot(tmp=="    ")
    
    coffs <- character(6)
    for(i in 1:length(coffs))
        coffs[i] <- readChar(con=con, nchars=8)
    
    ioffs <- c(as.double(version), as.integer(coffs), as.integer(start))
    names(ioffs) <- c("FCSversion", "textstart", "textend", "datastart",
                      "dataend", "anastart", "anaend", "additional")
    ioffs[2:7] <- ioffs[2:7]+ioffs[8]
    
    if(all(is.na(ioffs[2:5]) || ioffs[2:5]==""))
        stop("Missing header information to start parsing the binary ",
             "section of the file")
    return(ioffs)
}



## ==========================================================================
## parse FCS file text section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCStext <- function(con, offsets,emptyValue)
{

    seek(con, offsets["textstart"])
    ## Certain software (e.g. FlowJo 8 on OS X) likes to put characters into
    ## files that readChar can't read, yet readBin, rawToChar and iconv can 
    ## handle just fine.
    txt <- readBin(con,"raw", offsets["textend"]-offsets["textstart"]+1)
    txt <- iconv(rawToChar(txt), "", "latin1", sub="byte")
#	browser()
	if(offsets["FCSversion"]<=2)##
	{
		delimiter <- substr(txt, 1, 1)
		sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter,
				fixed=TRUE)[[1]]
		rv <- c(offsets["FCSversion"], sp[seq(2, length(sp), by=2)])
		names(rv) <- gsub("^ *| *$", "", c("FCSversion", sp[seq(1, length(sp)-1, by=2)]))	
	}else
	{
		#only apply the patch parser to FCS3
		rv = fcs_text_parse(txt,emptyValue=emptyValue)
		rv = c(offsets["FCSversion"], rv)
		names(rv)[1] = "FCSversion" #not sure if this line is necessary	
		names(rv) <- gsub("^ *| *$", "", names(rv))#trim the leading and trailing whitespaces
	}
	
	
	
    return(rv)
}

## ==========================================================================
## a patch to fix the parsing issue when delimiter exists in the keyword value
##,which is allowed only when it is followed by another delimiter immediately
## Note that it is only applies to FCS3.0 because empty value is not valid value.   
##however,this does not conform to FCS2.0,so this patch only applies to FCS3.0 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fcs_text_parse = function(str,emptyValue) {
#	browser()
	pairs = c()
	
	div = substr(str, 1, 1)
	if(nchar(div) != 1) {
		return(c())
	} 
	
	# regexes require double-escaping (*sigh*)
	if(div == "\\") {
		div = "\\\\"
	}
	if(div=="|")
		div<-paste("\\",div,sep="")

	i = 2
	repeat {
		remaining = substr(str, i, nchar(str))
		if(nchar(gsub(" ","",remaining)) == 0) {##filter out white spaces before check the length
			break
		}
		

		key_end_regex <- paste("([^",div,"]",div,")", sep="")
		#############################################################################################
		##1.when empty value is allowed, we assume there is no double delimiters in any values,
		##  otherwise,it would break the parser
		##2.when empty value is not allowed, we safely parse the double delimiters as the valid values
		#############################################################################################
		if(emptyValue) 
		{
			value_end_regex <- div
			final_end_regex <- paste(div,"$", sep="")
		}else
		{
			value_end_regex <- paste("([^",div,"]",div,"[^",div,"])", sep="")
			final_end_regex <- paste("[^",div,"]",div,"$", sep="")			
		}	
			
		
#		browser()
		# find key end
		divider_index = regexpr(key_end_regex, remaining, perl=TRUE)
		if(divider_index < 0) {
			
			if(i==2)
				stop("ERROR: No second divider found\n")
			else
			{
				warning("keyword: ",remaining," is dropped because no value found!")
				break
			}
				
		}
		divider_index = divider_index + 1
		
		value_search = substr(remaining, divider_index + 1, nchar(remaining))
		# find value end
		value_end_index = regexpr(value_end_regex, value_search, perl=TRUE)
		if(emptyValue) 
			value_end_index<-value_end_index-1
			
		if(value_end_index < 0) {
			value_end_index = regexpr(final_end_regex, value_search, perl=TRUE)
			if(value_end_index < 0) {
				if(emptyValue)
					stop("No end found\n There could be double delimiter existing in keyword value.\nPlease set argument 'emptyValue' as FALSE and try again!")
				else
					stop("No end found\n There could be empty keyword value.\nPlease set argument 'emptyValue' as TRUE and try again!")
#				return(c())
#				break
			}
		}
		value_end_index = divider_index + value_end_index
		
		key = substr(remaining, 1, divider_index - 1)
		value = substr(remaining, divider_index + 1, value_end_index)
		
		#replace double delimiters with single one in the final output
		value = gsub(paste(div,div,sep=''), div, value)
		
		#    cat("key: ", key, "\n")
		#    cat("value: ", value, "\n")
		#    cat("-----------------------\n")
		
		pairs = c(pairs, value)
		names(pairs)[length(pairs)] = key
		
		i = i + value_end_index + 1
	}
	
	return(pairs)
}

##read odd bitwidth by reading raw and and operating on raw vector
.readFCSdataRaw<-function(con,dattype,count,size,signed,endian){
#	browser()
	if (size %in% c(1, 2, 4, 8))
	{
		readBin(con=con, what=dattype,n = count,size=size, signed=signed, endian=endian)
	}else
	{
		#read raw byte stream first
		nBytes<-count*size
		oldBytes <- readBin(con=con, what="raw",n = nBytes,size=1)
		#convert to bit vector
		oldBits<-rawToBits(oldBytes)
		#convert the data element to the non-odd  bitwidth
		oldBitWidth<-size*8
		newBitWidth<-2^ceiling(log(oldBitWidth,2))
		newBits<-unlist(lapply(1:count,function(i){
#							browser()
											start<-(i-1)*oldBitWidth+1
											c(raw(newBitWidth-oldBitWidth)
											,oldBits[start:(start+oldBitWidth-1)])
											}
								)
						)		
		#convert raw byte to corresponding type by readBin
		readBin(packBits(newBits,"raw"),what=dattype,n=count,size=newBitWidth/8, signed=signed, endian=endian)
		
		
	}
	
}




## ==========================================================================
## read FCS file data section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSdata <- function(con, offsets, x, transformation, which.lines,
                        scale, alter.names, decades, min.limit=-111) {
    endian <- switch(readFCSgetPar(x, "$BYTEORD"),
                     "4,3,2,1" = "big",
                     "2,1" = "big",
                     "1,2" = "little",     
                     "1,2,3,4" = "little",
                     stop(paste("Don't know how to deal with $BYTEORD",
                                readFCSgetPar(x, "$BYTEORD"))))
    
    dattype <- switch(readFCSgetPar(x, "$DATATYPE"),
                      "I" = "integer",
                      "F" = "numeric",
                      stop(paste("Don't know how to deal with $DATATYPE",
                                 readFCSgetPar(x, "$DATATYPE"))))
    
    if (readFCSgetPar(x, "$MODE") != "L")
        stop(paste("Don't know how to deal with $MODE",
                   readFCSgetPar(x, "$MODE")))
    
    nrpar    <- as.integer(readFCSgetPar(x, "$PAR"))
    nrowTotal <- as.integer(readFCSgetPar(x, "$TOT"))

    if( "transformation" %in% names(x) &&  x[["transformation"]] == "custom"){
       range <- sapply(seq_len(nrpar),function(k){
                x[[sprintf("flowCore_$P%sRmax", k)]]
             })
    } else {
       range <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep="")))
    }
    bitwidth <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
    bitwidth <- unique(bitwidth)
    
    if(length(bitwidth)!=1)
        stop("Sorry, I am expecting the bitwidth to be the same for all ",
             "parameters")
    
    ##for DATA segment exceeding 99,999,999 byte.
    if(offsets["FCSversion"] == 3){
        realOff <- offsets - offsets[8]
        datastart <- as.numeric(readFCSgetPar(x, "$BEGINDATA"))
        dataend <- as.numeric(readFCSgetPar(x, "$ENDDATA"))
        if(realOff["datastart"] != datastart && realOff["datastart"]== 0){
            offsets["datastart"] <-  datastart+offsets[8]
        }
        if(realOff["datastart"] != datastart && realOff["datastart"]!= 0){
            print(datastart)
            print(offsets["datastart"])
            stop("The HEADER and the TEXT segment define different ",
                 "starting point to read the data.")
        }
        
        if(realOff["dataend"] != dataend && realOff["dataend"]== 0){
            offsets["dataend"] <-  dataend+offsets[8]
        }
        if(realOff["dataend"] != dataend && realOff["dataend"]!= 0){
            stop("The HEADER and the TEXT segment define different ending ",
                 "point to read the data.")
        }
    }
    if(bitwidth==10 ){
        if(!gsub(" " ,"", tolower(readFCSgetPar(x, "$SYS"))) ==  "cxp")
            warning("Invalid bitwidth specification.\nThis is a known bug in Beckman ",
                    "Coulter's CPX software.\nThe data might be corrupted if produced ",
                    "by another software.", call.=FALSE)
        else
            warning("Beckma Coulter CPX data.\nCorrected for invalid bitwidth 10.",
                    call.=FALSE)
        bitwidth <- 16
    }
    size <- bitwidth/8
#    if (!size %in% c(1, 2, 4, 8))
#        stop(paste("Don't know how to deal with bitwidth", bitwidth))

    nwhichLines <- length(which.lines)
    ##Read all reports
    if(is.null(which.lines) || (nwhichLines >  nrowTotal)){
        if (nwhichLines >  nrowTotal){
            cat("Warning: the number of lines specified (",nwhichLines,
                ") is greater
                 than the number of collected events (",nrowTotal,
                "). All the events have been read. \n")
        }
        seek(con, offsets["datastart"])
		
		dat <- .readFCSdataRaw(con, dattype,
				count= (offsets["dataend"]-offsets["datastart"]+1)/size,
                       size=size, signed=FALSE, endian=endian)
#        dat <- readBin(con, dattype,
#                       n = (offsets["dataend"]-offsets["datastart"]+1)/size,
#                       size=size, signed=FALSE, endian=endian)
    }else {  ##Read n lines with or without sampling
        if(length(which.lines)==1)
            which.lines <- sample(seq_len(nrowTotal), which.lines)
        which.lines <- sort(which.lines)
        outrange <- length(which(which.lines > nrowTotal))
        if(outrange!=0)
            stop("Some or all the line indices specified are greater that the",
                 "number of collected events.\n")
        dat <- c()
        for (i in 1:length(which.lines)){
            startP <- offsets["datastart"] + (which.lines[i]-1) * nrpar * size
            endP   <-  startP + nrpar * size
            seek(con, startP)
			temp <- .readFCSdataRaw(con, dattype, count= (endP - startP+1)/size,
					size=size, signed=FALSE, endian=endian) 
#            temp <- readBin(con, dattype, n = (endP - startP+1)/size,
#                            size=size, signed=FALSE, endian=endian) 
            dat <- c(dat, temp)                 
        }
    }
    ## stopifnot(length(dat)%%nrpar==0)
    ## Do we want the function to bail out when the above condition is TRUE?
    ## Might be better to assume the data was ok up to this point and
    ## exit gracefully with a warning as done in the following lines...
    ld <- length(dat)
    if(ld %% nrpar != 0){
        dat <- dat[1:(ld - (ld %% nrpar))]
        warning("Error in reading data stream for file '",
                summary(con)$description, "'/nData may be truncated!")
    }

    
    ## apply bitmask for integer data
    if(dattype=="integer"){
        if(length(unique(range))==1)
        {
            usedBits <- log2(range[1])
            if(usedBits<bitwidth)
                dat <- dat %% (2^usedBits)
            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
        }
        else
        {
            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
            for(i in 1:ncol(dat))
            {
                usedBits <- log2(range[i])
                if(usedBits<bitwidth)
                    dat[,i] <- dat[,i] %% (2^usedBits)
            }
        }
    }
    else
    {
        dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
    }
        
    cn  <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
    colnames(dat) <- if(alter.names)  structure(make.names(cn),
                                                names=names(cn))else cn

    ## truncate data at max range
    if(is.na(x["transformation"]))
    {
        for(i in seq_len(ncol(dat)))
            dat[dat[,i]>range[i],i] <- range[i]
        if(!is.null(min.limit))
            dat[dat<min.limit] <- min.limit
    }

    ## transform or scale if necessary
    if(transformation)
    {
       ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""),
             strict=FALSE)
       noPnE <- is.na(ampliPar)
       if(any(noPnE))
       {
          warning("No '$PnE' keyword available for the following channels: ",
                paste(which(noPnE), collapse=", "), "\nUsing '0,0' as default.",
                call.=FALSE)
          ampliPar[noPnE] <- "0,0"
       }
       ampli <- do.call(rbind,lapply(ampliPar, function(x)
                   as.integer(unlist(strsplit(x,",")))))
       for (i in 1:nrpar){
          if(ampli[i,1] > 0){
             dat[,i] <- 10^((dat[,i]/(range[i]-1))*ampli[i,1])
             range[i] <- 10^ampli[i,1]
          }
          else
             range[i] <- range[i]-1

       }
    }
    if(scale){
        d = 10^decades	
        for(i in 1:nrpar)
            if(ampli[i,1] > 0){
               dat[,i] <- d*((dat[,i]-1)/(range[i]-1))
               range[i] <- d*(range[i]/range[i]-1)
            } else{
                dat[,i] <- d*((dat[,i])/(range[i]))
                range[i] <- d
            }
    }
    attr(dat, "ranges") <- range
    return(dat) 
}



## ==========================================================================
## read in flowSets, based on the prada function readCytoSet. Most of the
## actual flowSet construction has been moved to the coercion of lists to
## flowSets. This is to facilitate the construction of flowSets from within
## the R environment (e.g. the flowFrames are actually retrieved from the
## network or a database.)
## ---------------------------------------------------------------------------
read.flowSet <- function(files=NULL, path=".", pattern=NULL, phenoData,
                         descriptions, name.keyword, alter.names=FALSE,
                         transformation="linearize", which.lines=NULL,
                         column.pattern=NULL, decades=0,
                         sep="\t", as.is=TRUE, name, ncdf=FALSE, dataset=NULL,
                         ...)
{
    ## A frame of phenoData information
    phenoFrame <- NULL
    
    ## deal with the case that the phenoData is provided, either as
    ## character vector or as AnnotatedDataFrame.
    if(!missing(phenoData)) {
        if(is.character(phenoData) && length(phenoData) == 1){
            phenoData <- read.AnnotatedDataFrame(file.path(path, phenoData),
                                                 header = TRUE, sep=sep,
                                                 as.is=as.is, ...)
            ## the sampleNames of the Annotated data frame must match the
            ## file names and we try to guess them from the input
            fnams <- grep("name|file|filename", varLabels(phenoData),
                          ignore.case=TRUE)
            if(length(fnams)){
                fn <- as.character(unlist(pData(phenoData[,fnams[1]])))
                if(any(duplicated(fn)))
                    stop("The file names supplied as part of the ",
                         "phenoData are not unique", call.=FALSE)
                sampleNames(phenoData) <- fn
                pd <- pData(phenoData)
                pd[,fnams[1]] <- fn
                pData(phenoData) <- pd
            }
            phenoFrame <- phenoData
        }else if(is(phenoData,"AnnotatedDataFrame")){
            phenoFrame <- phenoData
        }else{if(!is.list(phenoData))
                  stop("Argument 'phenoData' must be of type 'list', ",
                       "'AnnotatedDataFrame' or a filename\n",
                       "of a text file containing the phenotypic information")
          }
    }
    
    ## go on and find the files
    if(!is.null(phenoFrame)) {
        if(!is.null(files))
            warning("Supplied file names will be ignored, ",
                    "using names in the phenoData slot instead.")
        file.names <- sampleNames(phenoFrame)
	files <- file.path(path, file.names)
      	if(!all(file.exists(files)))
            stop(paste("Not all files given by phenoData could be found in",
                       path))
        if(!"name" %in% varLabels(phenoFrame)){
            phenoFrame$name <- files
            varMetadata(phenoFrame)["name",] <- "Filename"
        }
    }else{
        ## if we haven't found files by now try to search according to
        ## 'pattern'
        if(is.null(files)) {
            files <- dir(path,pattern,full.names=TRUE)
            file.names <- dir(path,pattern,full.names=FALSE)
            if(length(files)<1)
                stop(paste("No matching files found in ",path))
        } else {
            if(!is.character(files))
                stop("'files' must be a character vector.")
            file.names <- basename(files) ## strip path from names
            if(path != ".")
                files <- file.path(path, files)    
        }
    }
    
    flowSet <- lapply(files, read.FCS, alter.names=alter.names,
                      transformation=transformation, which.lines=which.lines,
                      column.pattern=column.pattern,
                      decades=decades, ncdf=ncdf)
    ## Allows us to specify a particular keyword to use as our sampleNames
    ## rather than requiring the GUID or the filename be used. This is handy
    ## when something like SAMPLE ID is a more reasonable choice.
    ## Sadly reading the flowSet is a lot more insane now.
    if(!missing(name.keyword)){
        keys <- unlist(sapply(flowSet,keyword,name.keyword))
        if(is.null(keys))
            stop("'", name.keyword, "' is not a valid keyword in any of ",
                 "the available FCS files.", call.=FALSE)
        if(length(keys) != length(flowSet))
            stop("One or several FCS files do not contain a keyword '",
                 name.keyword, "'.", call.=FALSE)
        if(any(duplicated(keys)))
            stop("The values of '", name.keyword, "' are not unique.", call.=FALSE)
        names(flowSet) <- sapply(flowSet,keyword,name.keyword)
    }else{
        names(flowSet) <- make.unique(file.names)
    }
    flowSet <- as(flowSet,"flowSet")
    if(!is.null(phenoFrame))
        phenoData(flowSet) <- phenoFrame
    else if(!missing(phenoData)){
        ##Collect the names for each field in the data frame
        field.names <- names(phenoData)
        if(is.null(field.names))
            stop("phenoData list must have names")
        field.names <- sapply(seq(along=phenoData),function(i) {
            if(length(field.names[i]) == 0) as(phenoData[i],"character")
            else field.names[i]
        })
        if(!missing(descriptions)) {
            ##If the descriptions have names, reorder them as needed.
            if(!is.null(names(descriptions)))
                descriptions = descriptions[field.names]
        } else
        descriptions <- field.names
        names(phenoData) <- field.names
        oldpDat <- pData(flowSet)
        newpDat <- as.data.frame(keyword(flowSet, phenoData))
        sel <- intersect(colnames(newpDat), colnames(oldpDat))
        if(length(sel))
            oldpDat <- oldpDat[, -which(colnames(oldpDat)%in% sel), drop=FALSE]
        newpDat <- cbind(oldpDat, as.data.frame(keyword(flowSet,phenoData)))
        if(any(duplicated(newpDat$name)))
            stop("The character strings in the 'name' variable in the phenoData slot ",
                 "have to be unique.", call.=FALSE)
        phenoData(flowSet) <- new("AnnotatedDataFrame", data=newpDat,
                                  varMetadata=data.frame(labelDescription=I(colnames(newpDat)),
                                                         row.names=colnames(newpDat)))
    }
    ## finally decide on which names to use for the sampleNames, but retain the
    ## original GUIDs in case there are some
    guids <- unlist(fsApply(flowSet, identifier))
    if(any(guids=="anonymous") || !missing(name.keyword))
        guids <- sampleNames(flowSet)
    if(any(duplicated(guids)))
        guids <- make.unique(guids)
    if("GUID" %in% names(description(flowSet[[1]])))
        flowSet <- fsApply(flowSet, function(x){
            keyword(x) <- c(GUID.original=as.character(keyword(x, "GUID")))
            x
        })
    sampleNames(flowSet) <- guids
    if(!missing(name))
        identifier(flowSet) <- name
    flowSet
}



write.AnnotatedDataFrame <- function(frame, file)
{
    con <- file(file, "w")
    on.exit(close(con))
    writeLines(paste(rep("#", ncol(frame)), varLabels(frame),
                     rep(": ", ncol(frame)),
                     varMetadata(frame)$labelDescription, sep=""), con)
    write.table(pData(frame), con, quote=FALSE, sep="\t")  
}


cleanup <- function() if(file.exists(".flowCoreNcdf"))
    unlink(".flowCoreNcdf", recursive=TRUE)





## ==========================================================================
## write new FCS file header
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
writeFCSheader <- function(con, offsets)
{
    seek(con, 0)
    writeChar("FCS3.0    ", con, eos=NULL)
    len <- length(offsets)/2
    for (i in seq_len(len)) {
        indx <- 2*(i-1) +1;
        val1 <- offsets[indx]
        val2 <- offsets[indx+1];
        st1 <- 8 - nchar(val1)
        st2 <- 8 - nchar(val2)
        if( nchar(val1) > 8 || nchar(val2) > 8){
             val1 <- val2 <- 0
        }
        writeChar(paste(paste(rep(" ", 8 - nchar(val1)), collapse=""), val1,
                        collapse="", sep=""), con, eos=NULL)
        writeChar(paste(paste(rep(" ", 8 - nchar(val2)), collapse=""), val2,
                        collapse="", sep=""), con, eos=NULL)
    }
     invisible()
}



## ==========================================================================
## collapse the content of the description slot into a character vector
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
collapseDesc <- function(x)
{
    d <- description(x)
	##make sure there is no empty value for each keyword in order to conform to FCS3.0
#	browser()
	d <- lapply(d, function(y){
				if(length(y)==0)
					return(" ")
				else
					return(sub("^$"," ",y)) 
						
			})
    d <- d[order(names(d))]
	spillName <- intersect(c("SPILL", "spillover"), names(d))
    if(length(spillName) >0){
		mat <-  d[[spillName]]
		rNum <- as.character(nrow(mat))
		clNames <- paste(colnames(mat),sep=",")
		vec <- paste(c(t(mat)),sep=",",collapse=",")
	    d[spillName] <- paste(c(rNum,clNames,vec),sep=",",collapse=",")
	}
    paste("\\", iconv(paste(names(d), "\\", sapply(d, paste, collapse=" "),
                            "\\", collapse="", sep=""), to="latin1",
                      sub=" "), sep="")
   
}



## ==========================================================================
## Write a flowFrame to an FCS file. Unless given explicitely, the filename
## is taken from the idnetifier of the flowFrame. 'what' controls the output
## data type.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
write.FCS <- function(x, filename, what="numeric")
{
    ## Some sanity checking up front
    if(missing(filename))
    {
        filename <- identifier(x)
        if(!length(grep(".", filename, fixed=TRUE)))
            filename <- paste(filename, "fcs", sep=".")
    }
    what <- match.arg(what, c("integer", "numeric", "double"))
    if(!is.character(filename) || length(filename)!=1)
        stop("Argument 'filename' has to be a character of length 1.")
    if(!is(x, "flowFrame"))
        stop("Argument 'x' has to be a 'flowFrame'.") 
    ## Set the mandatory keywords first
    begTxt <- 58
    types <- data.frame(symbol=c("I", "F", "D"),
                        bitwidth=c(2,4,8), stringsAsFactors=FALSE)
    rownames(types) <- c("integer", "numeric", "double")
    orders <- c(little="1,2,3,4", big="4,3,2,1")
    endian <- "big"
    mk <- list("$BEGINANALYSIS"="0",
               "$BEGINDATA"="0",
               "$BEGINSTEXT"=0,
               "$BYTEORD"=orders[endian],
               "$DATATYPE"=types[what, "symbol"],
               "$ENDANALYSIS"="0",
               "$ENDDATA"="0",
               "$ENDSTEXT"="0",
               "$MODE"="L",
               "$NEXTDATA"="0",
               "$PAR"=ncol(x),
               "$TOT"=nrow(x),
               "FCSversion"="3")
    pnb <- as.list(rep(types[what, "bitwidth"]*8, ncol(x)))
    names(pnb) <- sprintf("$P%sB", 1:ncol(x))
    mk <- c(mk, pnb)
#	browser()
    ## FlowJo seems to get confused by empty values in PnS, we fix that here
    pns <- description(x)[sprintf("$P%sS", 1:ncol(x))]
    names(pns) <- sprintf("$P%sS", 1:ncol(x))
#    pns <- lapply(pns, function(y) if(!length(y)) " " else y)
    mk <- c(mk, pns)
    ## We need all PnE keywords and assume "0,0" if they are missing
    pne <- description(x)[sprintf("$P%sE", 1:ncol(x))]
    names(pne) <- sprintf("$P%sE", 1:ncol(x))
    pne[sapply(pne, length)==0] <- "0,0"
    mk <- c(mk, pne)
    ## The same for PnR, "1024" if missing
    pnr <- description(x)[sprintf("$P%sR", 1:ncol(x))]
    names(pnr) <- sprintf("$P%sR", 1:ncol(x))
    pnr[sapply(pnr, length)==0] <- "1024"
    mk <- c(mk, pnr)
    ## Now update the PnN keyword
    pnn <- colnames(x)
    names(pnn) <- sprintf("$P%sN", 1:ncol(x))
    mk <- c(mk, pnn)
    description(x) <- mk
    ## Figure out the offsets based on the size of the initial text section
    ld <-  length(exprs(x)) * types[what, "bitwidth"]
    ctxt <- collapseDesc(x)
    endTxt <- nchar(ctxt) + begTxt -1
    endDat <- ld + endTxt
    endTxt <- endTxt +(nchar(endTxt+1)-1) + (nchar(endDat)-1)
    ## Now we update the header with the new offsets and recalculate
#	browser()
    endDat <- ld + endTxt
    description(x) <- list("$BEGINDATA"=endTxt+1,
                           "$ENDDATA"=endTxt+ld)
    ctxt <- collapseDesc(x)
	
    offsets <- c(begTxt, endTxt, endTxt+1, endTxt+ld, 0,0)
    ## Write out to file
    con <- file(filename, open="wb")
    on.exit(close(con))
    writeFCSheader(con, offsets)
    writeChar(ctxt, con, eos=NULL)
    writeBin(as(t(exprs(x)), what), con, size=types[what, "bitwidth"],
             endian=endian)
	writeChar("00000000", con, eos=NULL)
    filename
}



## ==========================================================================
## Write separate FCS files for each flowFrame in a flowSet. Filenames are
## either taken from the identifiers of the flowSets, or from a character
## vector with the same length as the flowSet, or by prepending 'i_' to the
## value of a character scalar. By default, the files will be written in a
## subdirectry constructed from the flowSet's identifier.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
write.flowSet <- function(x, outdir=identifier(x), filename, ...)
{
    ## Some sanity  checking up front
    if(missing(filename))
    {
        filename <- as.character(fsApply(x, identifier))
    }
    else
    {
        ferr <- paste("Argument 'fileame has to be a character scalar or",
                      "a character vector of the same length as 'x'.")
        if(is.character(filename))
        {
            if(length(filename)==1)
                filename <- paste(seq_len(length(x)), filename, sep="_")
            else if(length(filename) != length(x))
                stop(ferr)
        }
        else
            stop(ferr)
    }
    if(!is.character(outdir) || length(outdir)!=1)
        stop("Argument 'outdir' has to be a character of length 1.")
    if(!file.exists(outdir))
        dir.create(outdir, recursive=TRUE)
    if(!is(x, "flowSet"))
        stop("Argument 'x' has to be a 'flowSet'.")
    for(f in seq_len(length(x)))
    {
        if(!length(grep(".", filename[f], fixed=TRUE)))
            filename[f] <- paste(filename[f], "fcs", sep=".")
        write.FCS(x[[f]], filename=file.path(outdir, filename[f]), ...)
    }
    sampleNames(x) <- filename
    pData(x)$FCS_File <- filename
    write.AnnotatedDataFrame(phenoData(x), file=file.path(outdir, "annotation.txt"))
    outdir
}
