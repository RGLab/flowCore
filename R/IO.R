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
                     which.lines=NULL,
                     debug=FALSE,
                     alter.names=FALSE,
                     column.pattern=NULL,decades=0)
{
    ## check file name
    if(!is.character(filename) ||  length(filename)!=1)
        stop("'filename' must be character skalar")
    if(!file.exists(filename))
        stop(paste("'", filename, "' is not a valid file", sep=""))
    con <- file(filename, open="rb")

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
    offsets <- readFCSheader(con)
    txt     <- readFCStext(con, offsets, debug)
    mat     <- readFCSdata(con, offsets, txt, transformation, which.lines,
                           debug, scale, alter.names, decades)
    params  <- makeFCSparameters(colnames(mat),txt, transformation, scale,
                                 decades)
    close(con)

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

    ## set transformed flag
    txt[["FILENAME"]] = filename
    if(transformation==TRUE){
        txt[["transformation"]] <-"applied" 
    }

    ## build description from FCS parameters
    description <- strsplit(txt,split="\n")
    names(description) <- names(txt)

    ## the spillover matrix
    spID <- intersect(c("SPILL", "spillover"), names(description))
    if(length(spID)>0){
        sp <- description[[spID]]
        nrCols <- as.numeric(substr(sp,1,1))
        sp <- substr(sp,3,nchar(sp))
        cnames <- strsplit(sp, ",")[[1]][1:nrCols]
        vals <- as.numeric(strsplit(sp, ",")[[1]][(nrCols+1):((nrCols*nrCols))])
        spmat <- matrix(vals, ncol=nrCols, byrow=TRUE)
        colnames(spmat) <- cnames
        description[[spID]] <- spmat
    }
    
    return(new("flowFrame", exprs=mat, description= description,
               parameters=params))
}


## ==========================================================================
## create AnnotatedDataFrame describing the flow parameters (channels)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
makeFCSparameters <- function(cn,txt, transformation, scale, decades) {

    npar <- length(cn)
    id <- paste("$P",1:npar,sep="")
    range <- origRange <- as.numeric(txt[paste(id,"R",sep="")])
    range <- rbind(0,range-1)

    ## make sure the ranges are transformed along with the data
    if(transformation & !scale){
        ampliPar <- txt[paste(id,"E",sep="")]
        ampli <- do.call("rbind",lapply(ampliPar, function(x)
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
        range=origRange, minRange=range[1,], maxRange=range[2,]),
        varMetadata=data.frame(row.names=I(c("name","desc","range",
                               "minRange", "maxRange")),
        labelDescription=I(c("Name of Parameter","Description of Parameter",
        "Range of Parameter", "Minimum Parameter Value after Transforamtion",
        "Maximum Parameter Value after Transformation"))))
}


## ==========================================================================
## match FCS parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSgetPar <- function(x, pnam)
{

    stopifnot(is.character(x), is.character(pnam)) 
    i <- match(pnam, names(x))
    if(any(is.na(i)))
        stop(paste("Parameter(s)", pnam, "not contained in 'x'"))
    return(x[i])
}


## ==========================================================================
## parse FCS file header
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSheader <- function(con)
{
    
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
    names(ioffs) <- c("FCSversion", "textstart", "textend", "datastart",
                      "dataend", "anastart", "anaend")
    
    if(all(is.na(ioffs[2:5]) || ioffs[2:5]==""))
        stop("Missing header information to start parsing the binary ",
             "section of the file")
    return(ioffs)
}


## ==========================================================================
## parse FCS file text section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCStext <- function(con, offsets, debug)
{

    seek(con, offsets["textstart"])
    ##Certain software (e.g. FlowJo 8 on OS X) likes to put characters into
    ##files that readChar can't read, yet readBin, rawToChar and iconv can 
    ##handle just fine.
    txt <- readBin(con,"raw", offsets["textend"]-offsets["textstart"]+1)
    txt <- iconv(rawToChar(txt), "", "latin1", sub="byte")
    delimiter <- substr(txt, 1, 1)
    sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter,
                    fixed=TRUE)[[1]]
    rv <- c(offsets["FCSversion"], sp[seq(2, length(sp), by=2)])
    names(rv) <- c("FCSversion", sp[seq(1, length(sp)-1, by=2)])
    return(rv)
}


## ==========================================================================
## read FCS file data section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSdata <- function(con, offsets, x, transformation,  which.lines, debug,
                        scale, alter.names, decades) {
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
    range    <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep="")))
    bitwidth <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
    bitwidth <- unique(bitwidth)
    
    if(length(bitwidth)!=1)
        stop("Sorry, I am expecting the bitwidth to be the same for all ",
             "parameters")
    
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
            stop("The HEADER and the TEXT segment define different ",
                 "starting point to read the data.")
        }
        
        if(offsets["dataend"] != dataend && offsets["dataend"]== 0){
            offsets["dataend"] <-  dataend
        }
        if(offsets["dataend"] != dataend && offsets["dataend"]!= 0){
            stop("The HEADER and the TEXT segment define different ending ",
                 "point to read the data.")
        }
    }
    
    
    size <- bitwidth/8
    if (!size %in% c(1, 2, 4, 8))
        stop(paste("Don't know how to deal with bitwidth", bitwidth))

    
    nwhichLines = length(which.lines)
    ##Read all reports
    if(is.null(which.lines) || (nwhichLines >  nrowTotal)){
        if (nwhichLines >  nrowTotal){
            cat("Warning: the number of lines specified (",nwhichLines,
                ") is greater
                 than the number of collected events (",nrowTotal,
                "). All the events have been read. \n")
        }
        seek(con, offsets["datastart"])
        dat <- readBin(con, dattype,
                       n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                       size=size, signed=FALSE, endian=endian)
    }else {  ##Read n lines with or without sampling 
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
            temp <- readBin(con, dattype, n = (endP - startP+1)/size,
                            size=size, signed=FALSE, endian=endian) 
            dat <- c(dat, temp)                 
        }
    }
    ## stopifnot(length(dat)%%nrpar==0)
    ## Do we want the function to bail out when the above condition is TRUE?
    ## Might be better to assume the data was ok up to this point and
    ## exit gracefully with a warning as done in the following lines...
    if(length(dat) %% nrpar != 0){
        dat <- dat[1:(length(dat) %/% nrpar)]
        warning("Error in reading data stream for file '",
                summary(con)$description, "'/nData may be truncated!")
    }

    ## apply bitmask for integer data
    if(dattype=="integer"){
        usedBits <- unique(log2(range))
        if(usedBits<bitwidth)
            dat <- dat %% (2^usedBits)
    }
    
    dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
    cn  <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
    colnames(dat) <- if(alter.names)  structure(make.names(cn),
                                                names=names(cn))else cn

    ## transform or scale if necessary
    if(transformation) {
        ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""))
        ampli <- do.call("rbind",lapply(ampliPar, function(x)
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
            if(ampli[i,1] > 0)
                dat[,i] <- d*((dat[,i]-1)/(range[i]-1))
            else
                dat[,i] <- d*((dat[,i])/(range[i]))
    }
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
                         descriptions, name.keyword,
                         sep="\t",...)
{
    ## A frame of phenoData information
    phenoFrame = NULL
    
    ## deal with the case that the phenoData is provided, either as
    ## character vector or as AnnotatedDataFrame.
    if(!missing(phenoData)) {
        if(is.character(phenoData) && length(phenoData) == 1){
            phenoData = read.AnnotatedDataFrame(file.path(path, phenoData),
            header = TRUE, as.is = TRUE, sep=sep, ...)
            ## the sampleNames of the Annotated data frame must match the
            ## file names and we try to guess them from the input
            fnams <- grep("name|file|filename", varLabels(phenoData),
                          ignore.case=TRUE)
            if(length(fnams))
                sampleNames(phenoData) <- unlist(pData(phenoData[,fnams[1]]))
        }else if(is(phenoData,"AnnotatedDataFrame")){
            phenoFrame = phenoData
        }else{if(!is.list(phenoData))
                  stop("Argument 'phenoData' must be of type ",
                       "'AnnotatedDataFrame' or a filename\n",
                       "of a text file containing the phenotypic information")
          }
    }
    
    ## go on and find the files
    if(!is.null(phenoFrame)) {
        if(!is.null(files))
            warning("Supplied file names will be ignored, ",
                    "using names in the phenoData slot instead.")
        file.names = sampleNames(phenoFrame)
        files      = dir(path,paste(gsub("\\.","\\\\\\.",file.names),
        collapse="|"),full.names=TRUE)
        if(length(files) != length(file.names)) 
            stop(paste("Not all files given by phenoData could be found in",
                       path))
    }else{
        
        ## if we haven't found files by now try to search according to
        ## 'pattern'
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
    }
    
    flowSet = lapply(files, read.FCS, ...)
    ## Allows us to specify a particular keyword to use as our sampleNames
    ## rather than requiring the filename be used. This is handy when something
    ## like SAMPLE ID is a more reasonable choice. Sadly reading the flowSet is
    ## a lot more insane now.
    if(!missing(name.keyword))
        names(flowSet) <- sapply(flowSet,keyword,name.keyword)
    else
        names(flowSet) <- file.names
    flowSet = as(flowSet,"flowSet")
    if(!is.null(phenoFrame))
        phenoData(flowSet) <- phenoFrame
    else if(!missing(phenoData)) {
        ##Collect the names for each field in the data frame
        field.names <- names(phenoData)
        print(field.names)
        if(is.null(field.names))
            stop("phenoData list must have names")
        field.names = sapply(seq(along=phenoData),function(i) {
            if(length(field.names[i]) == 0) as(phenoData[i],"character")
            else field.names[i]
        })
        if(!missing(descriptions)) {
            ##If the descriptions have names, reorder them as needed.
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

