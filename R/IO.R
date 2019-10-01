## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory


## ==========================================================================
## Determine which 'files' are valid FCS files
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
isFCSfile <- function(files)
{
    sapply(files, function(f){
        if (file.exists(f)) {
            con <- file(f, open="rb")
            on.exit(close(con))
            version <- readChar(con, 6)
            isTRUE(version %in% c("FCS2.0", "FCS3.0", "FCS3.1"))
        }
        else FALSE
    })
}


## ==========================================================================
## Reading FCS file header and TEXT section only
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' Read the TEXT section of a FCS file
#' 
#' Read (part of) the TEXT section of a Data File Standard for Flow Cytometry
#' that contains FACS keywords.
#' 
#' The function \code{read.FCSheader} works with the output of the FACS machine
#' software from a number of vendors (FCS 2.0, FCS 3.0 and List Mode Data LMD).
#' The output of the function is the TEXT section of the FCS files. The user
#' can specify some keywords to limit the output to the information of
#' interest.
#' 
#' @name read.FCSheader
#' 
#' @param files Character vector of filenames.
#' @param path Directory where to look for the files.
#' @param keyword An optional character vector that specifies the FCS keyword
#' to read.
#' @param ... other arguments passed to \code{link[flowCore]{read.FCS}}
#' 
#' @return A list of character vectors. Each element of the list correspond to
#' one FCS file.
#' @author N.Le Meur
#' @seealso \code{link[flowCore]{read.flowSet}},
#' \code{link[flowCore]{read.FCS}}
#' @keywords IO
#' @examples
#' 
#' samp <- read.FCSheader(system.file("extdata",
#'    "0877408774.B08", package="flowCore"))
#' samp
#' 
#' samp <- read.FCSheader(system.file("extdata",
#'    "0877408774.B08", package="flowCore"), keyword=c("$DATE", "$FIL"))
#' samp
#' 
#' @export
read.FCSheader <- function(files, path=".", keyword=NULL, ...)
{
  
    stopifnot(is.character(files), length(files)>=1, files!="")
    filenames <- files
    if(path != ".")
        files = file.path(path, files)
    res <- lapply(files, function(file){
      
                              thisRes <- try(header(file, ...), silent = TRUE)
                              if(class(thisRes) == "try-error"){
                                stop(thisRes, file)
                              }else
                                thisRes
                  })
    if (!is.null(keyword))
        res <- lapply(res, function(x) x[keyword])
    names(res) <- filenames
    res
}

header <- function(files, ...){
    con <- file(files, open="rb")
    offsets <- findOffsets(con, ...)
    txt     <- readFCStext(con, offsets, ...)
    close(con)
    txt
}


## ==========================================================================
## main wrapper to read FCS files
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' Read an FCS file
#' 
#' Check validity and Read Data File Standard for Flow Cytometry
#' 
#' 
#' The function \code{isFCSfile} determines whether its arguments are valid FCS
#' files.
#' 
#' The function \code{read.FCS} works with the output of the FACS machine
#' software from a number of vendors (FCS 2.0, FCS 3.0 and List Mode Data LMD).
#' However, the FCS 3.0 standard includes some options that are not yet
#' implemented in this function. If you need extensions, please let me know.
#' The output of the function is an object of class \code{flowFrame}.
#' 
#' For specifications of FCS 3.0 see \url{http://www.isac-net.org} and the file
#' \url{../doc/fcs3.html} in the \code{doc} directory of the package.
#' 
#' The \code{which.lines} arguments allow you to read a subset of the record as
#' you might not want to read the thousands of events recorded in the FCS file.
#' It is mainly used when there is not enough memory to read one single FCS
#' (which probably will not happen).  It will probably take more time than
#' reading the entire FCS (due to the multiple disk IO).
#' @name read.FCS
#' @aliases read.FCS cleanup isFCSfile
#' @usage 
#' isFCSfile(files)
#' 
#' read.FCS(filename, transformation="linearize", which.lines=NULL,
#'          alter.names=FALSE, column.pattern=NULL, invert.pattern = FALSE,
#'          decades=0, ncdf = FALSE, min.limit=NULL, 
#'          truncate_max_range = TRUE, dataset=NULL, emptyValue=TRUE, 
#'          channel_alias = NULL, ...)
#'          
#' @param files A vector of filenames
#' @param filename Character of length 1: filename
#' @param transformation An character string that defines the type of
#' transformation. Valid values are \code{linearize} (default),
#' \code{linearize-with-PnG-scaling}, or \code{scale}.  The \code{linearize}
#' transformation applies the appropriate power transform to the data. The
#' \code{linearize-with-PnG-scaling} transformation applies the appropriate
#' power transform for parameters stored on log scale, and also a linear
#' scaling transformation based on the 'gain' (FCS \$PnG keywords) for
#' parameters stored on a linear scale. The \code{scale} transformation scales
#' all columns to $[0,10^decades]$.  defaulting to decades=0 as in the FCS4
#' specification.  A logical can also be used: \code{TRUE} is equal to
#' \code{linearize} and \code{FALSE}(or \code{NULL}) corresponds to no
#' transformation.  Also when the transformation keyword of the FCS header is
#' set to "custom" or "applied", no transformation will be used.
#' @param which.lines Numeric vector to specify the indices of the lines to be
#' read. If NULL all the records are read, if of length 1, a random sample of
#' the size indicated by \code{which.lines} is read in. It's used to achieve partial disk IO
#' for the large FCS that can't fit the full data into memory. Be aware the potential slow read
#' (especially for the large size of random sampling) due to the frequent disk seek operations. 
#' @param alter.names boolean indicating whether or not we should rename the
#' columns to valid R names using \code{\link{make.names}}. The default is
#' FALSE.
#' @param column.pattern An optional regular expression defining parameters we
#' should keep when loading the file. The default is NULL.
#' @param invert.pattern logical. By default, \code{FALSE}. If \code{TRUE},
#' inverts the regular expression specified in \code{column.pattern}. This is
#' useful for indicating the channel names that we do not want to read. If
#' \code{column.pattern} is set to \code{NULL}, this argument is ignored.
#' @param decades When scaling is activated, the number of decades to use for
#' the output.
#' @param ncdf Deprecated. Please use 'ncdfFlow' package for cdf based storage.
#' @param min.limit The minimum value in the data range that is allowed. Some
#' instruments produce extreme artifactual values. The positive data range for
#' each parameter is completely defined by the measurement range of the
#' instrument and all larger values are set to this threshold. The lower data
#' boundary is not that well defined, since compensation might shift some
#' values below the original measurement range of the instrument. This can be 
#' set to an arbitrary number or to \code{NULL} (the default value), in which 
#' case the original values are kept.
#' @param truncate_max_range logical type. Default is TRUE. can be optionally
#' turned off to avoid truncating the extreme positive value to the instrument
#' measurement range .i.e.'$PnR'.
#' @param dataset The FCS file specification allows for multiple data segments
#' in a single file. Since the output of \code{read.FCS} is a single
#' \code{flowFrame} we can't automatically read in all available sets. This
#' parameter allows to chose one of the subsets for import. Its value is
#' supposed to be an integer in the range of available data sets. This argument
#' is ignored if there is only a single data segment in the FCS file.
#' @param emptyValue boolean indicating whether or not we allow empty value for
#' keyword values in TEXT segment.  It affects how the double delimiters are
#' treated.  IF TRUE, The double delimiters are parsed as a pair of start and
#' end single delimiter for an empty value.  Otherwise, double delimiters are
#' parsed one part of string as the keyword value.  default is TRUE.
#' @param channel_alias a data.frame used to provide the alias of the channels
#' to standardize and solve the discrepancy across FCS files.It is expected to
#' contain 'alias' and 'channels' column. Each row/entry specifies the common
#' alias name for a collection of channels (comma separated). See examples for
#' details.
#' @param ... ignore.text.offset: whether to ignore the keyword values in TEXT
#' segment when they don't agree with the HEADER.  Default is FALSE, which
#' throws the error when such discrepancy is found.  User can turn it on to
#' ignore TEXT segment when he is sure of the accuracy of HEADER so that the
#' file still can be read.
#' 
#' @return
#' 
#' \code{isFCSfile} returns a logical vector.
#' 
#' \code{read.FCS} returns an object of class
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} that contains the data in
#' the \code{exprs} slot, the parameters monitored in the \code{parameters}
#' slot and the keywords and value saved in the header of the FCS file.
#' 
#' @author F. Hahne, N.Le Meur
#' @seealso \code{\link[flowCore]{read.flowSet}}
#' @keywords IO
#' @examples
#' 
#' ## a sample file
#' fcsFile <- system.file("extdata", "0877408774.B08", package="flowCore")
#' 
#' ## read file and linearize values
#' samp <-  read.FCS(fcsFile, transformation="linearize")
#' exprs(samp[1:3,])
#' description(samp)[3:6]
#' class(samp)
#' 
#' ## Only read in lines 2 to 5
#' subset <- read.FCS(fcsFile, which.lines=2:5, transformation="linearize")
#' exprs(subset)
#' 
#' ## Read in a random sample of 100 lines
#' subset <- read.FCS(fcsFile, which.lines=100, transformation="linearize")
#' nrow(subset)
#' 
#' #manually supply the alias vs channel options mapping as a data.frame
#' map <- data.frame(alias = c("A", "B")
#'                   , channels = c("FL2", "FL4")
#' )
#' fr <- read.FCS(fcsFile, channel_alias = map)
#' fr
#' 
#' @export
read.FCS <- function(filename,
                     transformation="linearize",
                     which.lines=NULL,
                     alter.names=FALSE,
                     column.pattern=NULL,
                     invert.pattern = FALSE,
                     decades=0,
                     ncdf=FALSE,
                     min.limit=NULL,
                     truncate_max_range = TRUE,
                     dataset=NULL,
                     emptyValue=TRUE
		                , channel_alias = NULL
                    , ...)
{
  channel_alias <- check_channel_alias(channel_alias)    
  if(ncdf)
    .Deprecated("'ncdf' argument is deprecated!Please use 'ncdfFlow' package for disk-based data structure.")
    ## check file name
    if(!is.character(filename) ||  length(filename)!=1)
        stop("'filename' must be character scalar")
    if(!file.exists(filename))
        stop(paste("'", filename, "' is not a valid file", sep=""))
    con <- file(filename, open="rb")
    on.exit(close(con))

    ## transform or scale data?
    fcsPnGtransform <- FALSE
    if(is.logical(transformation) && transformation ||
       !is.null(transformation) && transformation == "linearize") {
        transformation <- TRUE
        scale <- FALSE
    } else if ( !is.null(transformation) && transformation == "scale") {
        transformation <- TRUE
        scale <- TRUE
    } else if ( !is.null(transformation) && transformation == "linearize-with-PnG-scaling") {
        transformation <- TRUE
        scale <- FALSE
        fcsPnGtransform <- TRUE
    } else if (is.null(transformation) || is.logical(transformation) &&
               !transformation) {
        transformation <- FALSE
        scale <- FALSE
    }

    ## read the file
    offsets <- findOffsets(con,emptyValue=emptyValue, dataset = dataset, ...)

    txt <- readFCStext(con, offsets,emptyValue=emptyValue, ...)
    ## We only transform if the data in the FCS file hasn't already been
    ## transformed before
    if (fcsPnGtransform) txt[["flowCore_fcsPnGtransform"]] <- "linearize-with-PnG-scaling"
    if("transformation" %in% names(txt) &&
       txt[["transformation"]] %in% c("applied", "custom"))
       transformation <- FALSE
    mat <- readFCSdata(con, offsets, txt, transformation, which.lines,
                       scale, alter.names, decades, min.limit, truncate_max_range, channel_alias)
    matRanges <- attr(mat,"ranges")


    id <- paste("$P",1:ncol(mat),sep="")
    zeroVals <- as.numeric(sapply(strsplit(txt[paste(id,"E",sep="")], ","),
                                  function(x) x[2]))

    absMin <- colMins(mat,,na.rm=TRUE) # replace apply with matrixStats::colMins to speed up
    # absMin <- apply(mat,2,min,na.rm=TRUE)
    realMin <- pmin(zeroVals,pmax(-111, absMin, na.rm=TRUE), na.rm=TRUE)
    
    keep_idx <- seq_along(colnames(mat))
    remove_idx <- NULL
    # Only keep certain parameters
    if(!is.null(column.pattern)) {
      n <- colnames(mat)
      keep_idx <- grep(column.pattern, n, invert = invert.pattern)
      remove_idx <- setdiff(seq_along(colnames(mat)), keep_idx)
      cols <- names(attr(mat, "dimnames")[[2]])
      mat <- mat[,keep_idx,drop=FALSE]
      matRanges <- matRanges[keep_idx]
      names(attr(mat, "dimnames")[[2]]) <- cols[keep_idx]
      attr(mat, "ranges") <- matRanges
      absMin <- absMin[keep_idx]
      realMin <- realMin[keep_idx]
      zeroVals <- zeroVals[keep_idx]
      id <- id[keep_idx]
    }

	if("transformation" %in% names(txt) && txt[["transformation"]] == "custom") {
		for(i in seq_along(colnames(mat))) {
			realMin[i] <- as.numeric(txt[[sprintf("flowCore_$P%sRmin", keep_idx[i])]])
		}
	}

    params <- makeFCSparameters(colnames(mat),txt, transformation, scale,
                                 decades, realMin, id=keep_idx)

    ## check for validity
    if(is.null(which.lines)){
      total_number_of_events <- as.integer(readFCSgetPar(txt, "$TOT"));
        if(total_number_of_events != nrow(mat))
            stop("file", filename, "seems to be corrupted. \n The actual number of cells in data section ("
                 , nrow(mat), ") is not consistent with keyword '$TOT' (", total_number_of_events , ")")
    }

	## set transformed flag and fix the PnE and the Datatype keywords
    ## also add our own PnR fields.
    txt[["FILENAME"]] <- filename
    if(transformation==TRUE) {
       txt[["transformation"]] <-"applied"
       for(p in seq_along(pData(params)$name)) {
          txt[[sprintf("$P%sE", p)]] <- sprintf("0,%g", 0)
          txt[[sprintf("flowCore_$P%sRmax", keep_idx[p])]] <- matRanges[p] +1
          txt[[sprintf("flowCore_$P%sRmin", keep_idx[p])]] <- realMin[p]
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

  # Remove keywords for removed parameters
  if(!is.null(remove_idx)){
    remove_regex <- paste0("(\\$P|^P)", remove_idx, "[A-Z]+")
    remove_keys <- lapply(remove_regex, function(rx) grep(rx, names(description), value=TRUE))
    description <- description[!names(description) %in% do.call(c, remove_keys)]
  }

    ## the spillover matrix
    for(sn in .spillover_pattern){
      sp <- description[[sn]]
      if(!is.null(sp)){
          sp <- txt2spillmatrix(sp)
          if(is.matrix(sp))
          {
            cnames <- colnames(sp)
            cnames <- update_channel_by_alias(cnames, channel_alias, silent = TRUE)
            if(alter.names)
              cnames <- make.names(cnames)
            colnames(sp) <- cnames
            description[[sn]] <- sp  
          }
      }
    }
    tmp <- new("flowFrame", exprs=mat, description= description,
               parameters=params)
    identifier(tmp) <- basename(identifier(tmp))

    return(tmp)
}

txt2spillmatrix <- function(txt){
  splt <- strsplit(txt, ",")[[1]]
  nrCols <- as.numeric(splt[1])
  if(!is.na(nrCols)&&nrCols > 0)
  {
    cnames <- splt[2:(nrCols+1)]
    vals <- as.numeric(splt[(nrCols+2):length(splt)])
    matrix(vals, ncol=nrCols, byrow=TRUE, dimnames = list(NULL, cnames))
  }else
    txt
  
}
## ==========================================================================
## create AnnotatedDataFrame describing the flow parameters (channels)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
makeFCSparameters <- function(cn, txt, transformation, scale, decades,
                              realMin, id) {
    dattype <- switch(readFCSgetPar(txt, "$DATATYPE"),
        "I" = "integer",
        "F" = "numeric",
        "D" = "numeric",
        stop(paste("Don't know how to deal with $DATATYPE",
                readFCSgetPar(txt, "$DATATYPE"))))
    
    npar <- length(cn)
    if(missing(id)){
      id = 1:npar
    }
    id <- paste("$P", id ,sep="")
    
    range <- sapply(id, function(this_id){
      rid <- paste("flowCore_", this_id,"Rmax",sep="")
      original <- is.na(txt[rid])
      if(!original)
        as.numeric(txt[rid]) + 1
      else
        as.numeric(txt[paste(this_id,"R",sep="")])
    })


    origRange <- range
    range <- rbind(realMin,range-1)

    ## make sure the ranges are transformed along with the data
    if(transformation & !scale){

        ampliPar <- txt[paste(id,"E",sep="")]
        noPnE <- is.na(ampliPar)
        if(any(noPnE))
            ampliPar[noPnE] <- "0,0"
        ampli <- do.call(rbind,lapply(ampliPar, function(x)
                                        as.numeric(unlist(strsplit(x,",")))))
        for (i in 1:npar)
            if(ampli[i,1] > 0 && dattype == "integer")
            {
              if(ampli[i,2] == 0)
                ampli[i,2] = 1 #correct f2 value for legacy FCS
              range[,i] <- 10^(range[,i]/(origRange[i]-1)*ampli[i,1])*ampli[i,2]
            }
                
    }
    else if(scale)
        range[2,] <- rep(10^decades, npar)

    desc <- txt[paste(id,"S",sep="")]
    desc <- gsub("^\\s+|\\s+$", "", desc)#trim the leading and tailing whitespaces
    # replace the empty desc with NA
    desc <- sapply(desc, function(thisDesc){
            if(!nzchar(thisDesc))
              NA
            else
              thisDesc
          })
    new("AnnotatedDataFrame",
        data=data.frame(row.names=I(id),name=I(cn),
        desc=I(desc),
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
findOffsets <- function(con,emptyValue=TRUE, dataset = NULL, ...)
{
    offsets <- readFCSheader(con)
    offsets <- matrix(offsets, nrow = 1, dimnames = list(NULL, names(offsets)))
    txt <- readFCStext(con, offsets[1, ],emptyValue=emptyValue, ...)

    addOff <- 0

    if("$NEXTDATA" %in% names(txt)){
      nd <- as.numeric(txt[["$NEXTDATA"]])
    }else
      nd <- 0
    
    txt.list <- list(txt)
    i <- 1
    while(nd != 0)
    {
        i <- i + 1
        addOff <- addOff + nd
        offsets <- rbind(offsets, readFCSheader(con, addOff))
        this.txt <- readFCStext(con, offsets[nrow(offsets),],emptyValue=emptyValue, ...)
        nd <- as.numeric(this.txt[["$NEXTDATA"]])
        txt.list[[i]] <- this.txt
    }

    ## check for multiple data sets
    nDataset <- length(txt.list)
    if(nDataset == 1)
      dataset <- 1
    else
    {
      
      if(is.null(dataset))
      {
        warning(sprintf("The file contains %d additional data segment%s.\n",
                nDataset-1, ifelse(nDataset>2, "s", "")),
            "The default is to read the first segment only.\nPlease consider ",
            "setting the 'dataset' argument.", call.=FALSE)
        dataset <- 1
        
      }
    }
    
    if(!is.numeric(dataset) || !dataset %in% seq_len(nDataset))
      stop(sprintf("Argument 'dataset' must be an integer value in [1,%d].",
              nDataset))
    offsets <- offsets[dataset,]
    txt <- txt.list[[dataset]]
    
#    browser()
    offsets <- checkOffset(offsets, txt, ...)
    return(offsets)
}

#' Fix the offset when its values recorded in header and TEXT don't agree
#' @param offsets the named vector returned by \code{findOffsets}
#' @param x the text segmented returned by \code{readFCStext}
#' @param ignore.text.offset whether to ignore the offset info stored in TEXT segment
#' @param ... not used.
#' @return the updated offsets
checkOffset <- function(offsets, x, ignore.text.offset = FALSE, ...){
  ##for DATA segment exceeding 99,999,999 byte.
  if(offsets["FCSversion"] >= 3){
    realOff <- offsets - offsets[8]

    # Let's not be too strick here as unfortunatelly, some files exported from FlowJo
    # are missing the $BEGINDATA and $ENDDATA keywords and we still need to read those
    datastart <- as.numeric(readFCSgetPar(x, "$BEGINDATA", strict=FALSE))
    if (is.na(datastart)) {
      if (realOff["datastart"] != 0) {
        datastart = realOff["datastart"]
        warning("Missing the required $BEGINDATA keyword! Reading data based on information in the FCS HEADER only.", call.=FALSE)
      } else {
        stop("Don't know where the data segment begins, there was no $BEGINDATA keyword and the FCS HEADER does not say it either.")
      }
    }
    dataend <- as.numeric(readFCSgetPar(x, "$ENDDATA", strict=FALSE))
    if (is.na(dataend)) {
      if (realOff["dataend"] != 0) {
        dataend = realOff["dataend"]
        warning("Missing the required $ENDDATA keyword! Reading data based on information in the FCS HEADER only.", call.=FALSE)
      } else {
        stop("Don't know where the data segment ends, there was no $ENDDATA keyword and the FCS HEADER does not say it either.")
      }
    }

    # when both are present and they don't agree with each other
    if(realOff["datastart"] != datastart && realOff["datastart"]== 0){ #use the TEXT when header is 0
      offsets["datastart"] <-  datastart+offsets[8]
    }
    if(realOff["datastart"] != datastart && realOff["datastart"]!= 0){#trust the header when it is non-zero
      msg <- paste0("The HEADER and the TEXT segment define different starting point ("
                    , offsets["datastart"], ":", datastart
                    , ") to read the data.")
      if(ignore.text.offset)
        warning(msg, " The values in TEXT are ignored!")
      else
        stop(msg)
    }
    #both are present and they don't agree
    if(realOff["dataend"] != dataend && (realOff["dataend"]== 0 || realOff["dataend"]== 99999999)) {#use TEXT when either header is 0 or TEXT is 99999999
      offsets["dataend"] <-  dataend+offsets[8]
    }
    if(realOff["dataend"] != dataend && realOff["dataend"]!= 0 && realOff["dataend"]!= 99999999) {#otherwise trust the header
      msg <- paste0("The HEADER and the TEXT segment define different ending point ("
          , offsets["dataend"], ":", dataend
          , ") to read the data.")
      if(ignore.text.offset)
        warning(msg, " The values in TEXT are ignored!")
      else
        stop(msg)

    }
  }
  offsets
}

## ==========================================================================
## parse FCS file header
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSheader <- function(con, start=0)
{
    seek(con, start)
    version <- readChar(con, 6)
    if(!version %in% c("FCS2.0", "FCS3.0", "FCS3.1"))
        stop("This does not seem to be a valid FCS2.0, FCS3.0 or FCS3.1 file")

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

    if(any(is.na(ioffs[2:5]) | ioffs[2:5]==""))
        stop("Missing header information to start parsing the binary ",
             "section of the file")
    return(ioffs)
}



## ==========================================================================
## parse FCS file text section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCStext <- function(con, offsets,emptyValue = TRUE, cpp = TRUE, ...)
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
        if(cpp)
          rv = fcsTextParse(txt,emptyValue=emptyValue)
        else
		  rv = fcs_text_parse(txt,emptyValue=emptyValue)

        if(!"FCSversion"%in%names(rv))
		  rv <- c(offsets["FCSversion"], rv)

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
# 	if(div == "\\") {
# 		div = "\\\\"
# 	}
# 	if(div=="|")
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

# deprecated by cpp version
sortBytes1 <- function(bytes, byte_order){
  nBytes <- length(byte_order)
  nTotal <- length(bytes)
  nElement <- nTotal / nBytes
  #relative byte order for the entire vector
  byte_orders <- rep(byte_order + 1 , nElement) 
  #absolute byte order for the entire vector
  byte_orders <- byte_orders + rep((seq_len(nElement)-1) * 4, each = nBytes)
  #re-order 
  ind <- order(byte_orders)
  bytes[ind]  
}

# Read raw FCS data. This is a wrapper around .readFCSdataRaw that calls the
# function with at most .Machine$integer.max bytes each time.
.readFCSdataRawMultiple <- function(con, dattype, count, size, signed, endian,
                                    splitInt = FALSE, byte_order) {
  chunk_size <- floor(.Machine$integer.max / size)
  num_chunks <- ceiling(count / chunk_size)

  items <- c()

  for (chunk in seq(num_chunks)) {
    if (chunk == num_chunks) {
      chunk_count <- count - chunk_size * floor(count / chunk_size)
    } else {
      chunk_count <- chunk_size
    }

    items <-
      c(items, .readFCSdataRaw(con, dattype, chunk_count, size, signed, endian,
        splitInt, byte_order))
  }

  items
}

##read odd bitwidth by reading raw and and operating on raw vector
.readFCSdataRaw <- function(con, dattype, count, size, signed, endian,
                            splitInt = FALSE, byte_order) {
  nBytes <- count * size

  if (splitInt && (size != 4 || dattype != "integer")) {
    stop("'splitInt = TRUE' is only valid for uint32")
  }

	if (size %in% c(1, 2, 4, 8)) {
    if (splitInt) {
      # Reorder bytes for mixed endian.
      if (endian == "mixed") {
        byte_order <- as.integer(strsplit(byte_order, ",")[[1]]) - 1
        if(length(byte_order) != size) {
          stop("byte order not consistent with bidwidths")
        }

        bytes <- readBin(con = con, what = "raw", n = nBytes, size = 1)
        newBytes <- sortBytes(bytes, byte_order)
        con <- newBytes
        endian <- "little"
      }
          
      # Read uint32 as two uint16. Coerce count again to make sure that it's
      # within the integer limit.
      splitted <- readBin(con = con, what = dattype, n = as.integer(count * 2),
                          size = size / 2, signed = FALSE, endian = endian)
      uint2double(splitted, endian == "big")
    } else {
      readBin(con = con, what = dattype, n = count, size = size,
              signed = signed, endian=endian)
    }
	} else {
		# Read raw byte stream first.
		oldBytes <- readBin(con = con, what = "raw", n = nBytes, size = 1)
		# Convert to bit vector.
		oldBits <- rawToBits(oldBytes)
		# Convert the data element to the non-odd bitwidth.
		oldBitWidth <- size * 8
		newBitWidth <- 2 ^ ceiling(log(oldBitWidth, 2))
		newBits <-
      unlist(lapply(1:count, function(i) {
				start <- (i - 1) * oldBitWidth + 1
        # Padding zeros.
        c(oldBits[start:(start + oldBitWidth - 1)],
          raw(newBitWidth-oldBitWidth))
      }))
		# Convert raw byte to corresponding type by readBin. packBits is
    # least-significant bit first, so we need to make sure endian is set to
    # "little" instead of the endian used in original FCS.
		readBin(packBits(newBits, "raw"), what = dattype, n = count,
            size = newBitWidth / 8, signed = signed, endian = "little")
	}
}

check_channel_alias <- function(channel_alias){
  if(is.null(channel_alias))
    return (new.env(parent = emptyenv()))
	else if(is(channel_alias, "data.frame"))
	{
	 if(!setequal(c("alias", "channels"), colnames(channel_alias)))
		stop("channel_alias must contain 'alias' and 'channels' columns")
 	 env <- new.env(parent = emptyenv())
	 apply(channel_alias, 1, function(row){
    		channels <- strsplit(split = ",", row["channels"])[[1]]
    		for(c in channels)
    		{
    		  c <- trimws(c)
    		 if(is.null(env[[c]]))
    		   env[[c]] <- trimws(row[["alias"]])
    		 else
    		  stop("multiple entries found in channel_alias for: ", c)
    		}
			})
	return (env)
	}else if(is(channel_alias, "environment"))
		return (channel_alias)
	else
	 stop("channel_alias must be either an environment or a data.frame")
}

update_channel_by_alias <- function(orig_chnl_names, channel_alias, silent = FALSE)
{
  keys <- ls(channel_alias)
  
  new_channels <- unlist(lapply(orig_chnl_names, function(col){
    alias <- channel_alias[[col]]
    if(is.null(alias))
    {
      #try partial match with case insensitive maching
      #escape special characters by enclosing it with \Q and \E
      
      ind <- unlist(lapply(keys, function(key){grepl(paste0("\\Q", key, "\\E"), col, ignore.case = TRUE)}))
      
      if(sum(ind)>1)
        stop(col, " is matched to the multiple entries in the channel_alias: ", paste(keys[ind], " "), "\n Try to modify channel_alias so that channel names are more specific! ")
      else if(sum(ind) == 0)
        return (col)
      else
        return (channel_alias[[keys[ind]]])
    }else
      return (alias)
  }))
  #validity check
  tb <- table(new_channels)
  is.dup <- tb > 1
  if(any(is.dup))
  {
    dup <- names(tb)[is.dup]
    dup.ind <- which(new_channels %in% dup)
    for(chnl in dup){
      chnl.ind <- which(new_channels == chnl)
      new_channels[chnl.ind] <- paste0(chnl, "-", seq(tb[chnl]))
    }
    dt <- data.frame(orig_channel_name = orig_chnl_names[dup.ind], new_channel_name = new_channels[dup.ind])
    if(!silent){
      print(dt)
      warning(paste0("channel_alias: Multiple channels from one FCS are matched to the same alias!\n",
                     "                 Integer suffixes added to disambiguate channels.\n",
                     "                 It is also recommended to verify correct mapping of spillover matrix columns."))
    }
  }
  return (new_channels)
  
}
## ==========================================================================
## read FCS file data section
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
readFCSdata <- function(con, offsets, x, transformation, which.lines,
                        scale, alter.names, decades, min.limit=-111, truncate_max_range = TRUE, channel_alias = NULL) {
    
    byte_order <- readFCSgetPar(x, "$BYTEORD")  
    endian <- switch(byte_order,
                     "4,3,2,1" = "big",
                     "2,1" = "big",
                     "1,2" = "little",
                     "1,2,3,4" = "little",
                     "mixed")
    
    dattype <- switch(readFCSgetPar(x, "$DATATYPE"),
                      "I" = "integer",
                      "F" = "numeric",
                      "D" = "numeric",
                      stop(paste("Don't know how to deal with $DATATYPE",
                                 readFCSgetPar(x, "$DATATYPE"))))

    if (readFCSgetPar(x, "$MODE") != "L")
        stop(paste("Don't know how to deal with $MODE",
                   readFCSgetPar(x, "$MODE")))

    nrpar    <- as.integer(readFCSgetPar(x, "$PAR"))
    nrowTotal <- as.integer(readFCSgetPar(x, "$TOT"))

    if( "transformation" %in% names(x) &&  x[["transformation"]] == "custom"){
      range_str <- sapply(seq_len(nrpar),function(k){
                x[[sprintf("flowCore_$P%sRmax", k)]]
             })
    } else {
        range_str <- readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep=""))
    }

    bitwidth_vec <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
    bitwidth <- unique(bitwidth_vec)
    multiSize <- length(bitwidth) > 1

    if(dattype=="numeric"&&multiSize)
      stop("Sorry, Numeric data type expects the same bitwidth for all parameters!")


    if(dattype=="integer"){
      suppressWarnings(range <- as.integer(range_str))

      #when any channel has range > 2147483647 (i.e. 2^31-1)
      if(any(is.na(range)))
        range <- as.numeric(range_str)

      if(any(is.na(range)))#throws if still fails
        stop("$PnR", range_str[is.na(range)][1], "is larger than R's integer limit:", .Machine$integer.max)
      else if(any(range>2^32)){
        #check if larger than C's uint32 limit ,which should be 2^32-1
        #but strangely(and inaccurately) these flow data uses 2^32 to specifiy the upper bound of 32 uint
        #we try to tolerate this and hopefully there is no such extreme value exsiting in the actual data section
        stop("$PnR", range_str[range>2^32][1], "is larger than C's uint32 limit:", 2^32-1)
      }

      if(multiSize){
        splitInt <- FALSE
      }else
      {
        splitInt <- bitwidth == 32
      }

    }
    else{
      splitInt <- FALSE
      range <- as.numeric(range_str)
      if(any(is.na(range)))
        stop("$PnR", range_str[is.na(range)][1], "is larger than R's numeric limit:", .Machine$double.xmax)
    }



    if(!multiSize){
      if(bitwidth==10){
        if(!gsub(" " ,"", tolower(readFCSgetPar(x, "$SYS"))) ==  "cxp")
          warning("Invalid bitwidth specification.\nThis is a known bug in Beckman ",
                  "Coulter's CPX software.\nThe data might be corrupted if produced ",
                  "by another software.", call.=FALSE)
        else
          warning("Beckma Coulter CPX data.\nCorrected for invalid bitwidth 10.",
                  call.=FALSE)
        bitwidth <- 16
      }
    }
    # multiSize <- T

    if(multiSize){
      size <- bitwidth_vec/8
      signed <- FALSE #dummy. not used in mutliSize logic.
    }else{
      size <- bitwidth/8

      # since signed = FALSE is not supported by readBin when size > 2
      # we set it to TRUE automatically then to avoid warning flooded by readBin
      # It shouldn't cause data clipping since we haven't found any use case where datatype is unsigned integer with size > 16bits
      signed <- !(size%in%c(1,2))
    }


    if(is.null(which.lines)){
      # Read the entire file.
      seek(con, offsets["datastart"])
      nBytes <- offsets["dataend"] - offsets["datastart"] + 1

	    if (multiSize) {
        if (nBytes > .Machine$integer.max) {
          stop(
            paste0("cannot import files with more than ", .Machine$integer.max,
                   " bytes in data segment when file has multiple bitwidths")
          )
        }
        nBytes <- as.integer(nBytes)

#	      if(splitInt&&dattype=="integer")
#	        stop("Mutliple bitwidths with big integer are not supported!")
	      if(endian == "mixed")
	        stop("Cant't handle diverse bitwidths while endian is mixed: ", byte_order)
	      
	      bytes <- readBin(con=con, what="raw",n = nBytes, size = 1)
	      # browser()
	      if(dattype == "numeric" && length(unique(size)) > 1)
	        stop("we don't support different bitwdiths for numeric data type!")
	      dat <- convertRawBytes(bytes, isInt = dattype == "integer", colSize = size, ncol = nrpar, isBigEndian = endian == "big")
	    } else {
	      dat <-
          .readFCSdataRawMultiple(
            con, dattype, count = nBytes / size, size = size, signed = signed,
            endian=endian, splitInt = splitInt, byte_order = byte_order)
	    }
    } else {
      # Read subset of lines, as selected by user.
      if (multiSize) {
        stop("'which.lines' cannot be used with multiple bitwidths")
      }

      # Verify that which.lines is positive and within file limit.
      if (length(which.lines) > 1) {
        if (any(which.lines < 0)) {
          warning("import will skip lines with negative indices")
          which.lines <- which.lines[which.lines > 0]
        }
        if (any(which.lines > nrowTotal)) {
          warning("import will skip lines over number of collected events")
          which.lines <- which.lines[which.lines < nrowTotal]
        }
      }

      if (length(which.lines) == 1) {
        # If a single value is given, sample N lines randomly.
        which.lines <- sample(seq(nrowTotal), which.lines)
      }

      which.lines <- sort(which.lines)
      dat <- c()
      for (i in 1:length(which.lines)){
        startP <- offsets["datastart"] + (which.lines[i] - 1) * nrpar * size
        endP   <- startP + nrpar * size
        seek(con, startP)
        temp <-
          .readFCSdataRawMultiple(
            con, dattype, count = as.integer(endP - startP + 1) / size,
            size = size, signed = signed, endian = endian, splitInt = splitInt,
            byte_order = byte_order)
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

            if(range[1]>0){
              usedBits <- ceiling(log2(range[1]))
              if(usedBits<bitwidth)
                dat <- dat %% (2^usedBits)
            }

            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
        }
        else
        {
            dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
            for(i in 1:ncol(dat))
            {

                if(range[i] > 0){
                  usedBits <- ceiling(log2(range[i]))
                  if(usedBits<bitwidth_vec[i])
                      dat[,i] <- dat[,i] %% (2^usedBits)
                }
            }
        }
    }
    else
    {
        dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
    }

    cn  <- readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))
    cn <- update_channel_by_alias(cn, channel_alias)
    cn <- if(alter.names)  structure(make.names(cn),names=names(cn)) else cn
    dimnames(dat) <- list(NULL, cn)
    ## truncate data at max range
    if(is.na(x["transformation"]))
    {
        if(truncate_max_range){
          for(i in seq_len(ncol(dat)))
            dat[dat[,i]>range[i],i] <- range[i]
        }

        if(!is.null(min.limit))
            dat[dat<min.limit] <- min.limit
    }

    ## Transform or scale if necessary
    # J.Spidlen, Nov 13, 2013: added the flowCore_fcsPnGtransform keyword, which is
    # set to "linearize-with-PnG-scaling" when transformation="linearize-with-PnG-scaling"
    # in read.FCS(). This does linearization for log-stored parameters and also division by
    # gain ($PnG value) for linearly stored parameters. This is how the channel-to-scale
    # transformation should be done according to the FCS specification (and according to
    # Gating-ML 2.0), but lots of software tools are ignoring the $PnG division. I added it
    # so that it is only done when specifically asked for so that read.FCS remains backwards
    # compatible with previous versions.
    fcsPnGtransform <- FALSE
    flowCore_fcsPnGtransform <- readFCSgetPar(x, "flowCore_fcsPnGtransform", strict=FALSE)
    if(!is.na(flowCore_fcsPnGtransform) && flowCore_fcsPnGtransform == "linearize-with-PnG-scaling") fcsPnGtransform <- TRUE
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
                   as.numeric(unlist(strsplit(x,",")))))
       PnGPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "G", sep=""), strict=FALSE)
       noPnG <- is.na(PnGPar)
       if(any(noPnG)) PnGPar[noPnG] <- "1"
       PnGPar = as.numeric(PnGPar)

       for (i in 1:nrpar){
          if(ampli[i,1] > 0 && dattype == "integer"){
             # J.Spidlen, Nov 5, 2013: This was a very minor bug. The linearization transformation
             # for $PnE != "0,0" is defined as:
             # For $PnR/r/, r>0, $PnE/f,0/, f>0: n is a logarithmic parameter with channel values
             # from 0 to r-1. A channel value xc is converted to a scale value xs as xs=10^(f*xc/r).
             # Note the "r" instead of the "r-1" in the formula (which would admitedly make more sense)
			 # However, this is the standard that apparently has been followed by BD and other companies
             # "forever" and it is therefore addoped as such by the ISAC DSTF (see FCS 3.1 specification)
             # To bring this to compliance, I am just changing
             # dat[,i] <- 10^((dat[,i]/(range[i]-1))*ampli[i,1])
             # to
             # dat[,i] <- 10^((dat[,i]/range[i])*ampli[i,1])
             # M.Jiang, Jan 28, 2018. Based on the standard, formula should be  xs = 10^(f1 * xc /(r)) * f2.
             if(ampli[i,2] == 0)
               ampli[i,2] = 1 #correct f2 value for legacy FCS
             dat[,i] <- 10^((dat[,i]/range[i])*ampli[i,1])*ampli[i,2]
             range[i] <- 10^ampli[i,1]*ampli[i,2]
          }
          else if (fcsPnGtransform && PnGPar[i] != 1) {
             dat[,i] <- dat[,i] / PnGPar[i]
             range[i] <- (range[i]-1) / PnGPar[i]
          }
          else
             range[i] <- range[i]-1

       }
    }
    if(scale){
        d = 10^decades
        for(i in 1:nrpar)
            if(ampli[i,1] > 0 && dattype == "integer"){
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
#' Read a set of FCS files
#' 
#' Read one or several FCS files: Data File Standard for Flow Cytometry
#' 
#' There are four different ways to specify the file from which data is to be
#' imported:
#' 
#' First, if the argument \code{phenoData} is present and is of class
#' \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}, then the
#' file names are obtained from its sample names (i.e. row names of the
#' underlying data.frame).  Also column \code{name} will be generated based on
#' sample names if it is not there. This column is mainly used by visualization
#' methods in flowViz.  Alternatively, the argument \code{phenoData} can be of
#' class \code{character}, in which case this function tries to read a
#' \code{AnnotatedDataFrame} object from the file with that name by calling
#' \code{\link[Biobase]{read.AnnotatedDataFrame}(file.path(path,phenoData),\dots{})}.
#' 
#' In some cases the file names are not a reasonable selection criterion and
#' the user might want to import files based on some keywords within the file.
#' One or several keyword value pairs can be given as the phenoData argument in
#' form of a named list.
#' 
#' Third, if the argument \code{phenoData} is not present and the argument
#' \code{files} is not \code{NULL}, then \code{files} is expected to be a
#' character vector with the file names.
#' 
#' Fourth, if neither the argument \code{phenoData} is present nor \code{files}
#' is not \code{NULL}, then the file names are obtained by calling
#' \code{dir(path, pattern)}.
#' 
#' @name read.flowSet
#' 
#' @usage
#' read.flowSet(files=NULL, path=".", pattern=NULL, phenoData,
#'              descriptions,name.keyword, alter.names=FALSE,
#'              transformation = "linearize", which.lines=NULL,
#'              column.pattern = NULL, invert.pattern = FALSE, decades=0, sep="\t",
#'              as.is=TRUE, name, ncdf=FALSE, dataset=NULL, min.limit=NULL,
#'              truncate_max_range = TRUE, emptyValue=TRUE, 
#'              ignore.text.offset = FALSE, channel_alias = NULL, \dots)
#' 
#' @param files Optional character vector with filenames.
#' @param path Directory where to look for the files.
#' @param pattern This argument is passed on to
#' \code{\link[base:list.files]{dir}}, see details.
#' @param phenoData An object of class \code{AnnotatedDataFrame},
#' \code{character} or a list of values to be extracted from the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} object, see details.
#' @param descriptions Character vector to annotate the object of class
#' \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param name.keyword An optional character vector that specifies which FCS
#' keyword to use as the sample names. If this is not set, the GUID of the FCS
#' file is used for sampleNames, and if that is not present (or not unique),
#' then the file names are used.
#' @param alter.names see \code{\link[flowCore]{read.FCS}} for details.
#' @param transformation see \code{\link[flowCore]{read.FCS}} for details.
#' @param which.lines see \code{\link[flowCore]{read.FCS}} for details.
#' @param column.pattern see \code{\link[flowCore]{read.FCS}} for details.
#' @param invert.pattern see \code{\link[flowCore]{read.FCS}} for details.
#' @param decades see \code{\link[flowCore]{read.FCS}} for details.
#' @param sep Separator character that gets passed on to
#' \code{\link[Biobase]{read.AnnotatedDataFrame}}.
#' @param as.is Logical that gets passed on to
#' \code{\link[Biobase]{read.AnnotatedDataFrame}}. This controls the automatic
#' coercion of characters to factors in the \code{phenoData}slot.
#' @param name An optional character scalar used as name of the object.
#' @param ncdf Deprecated. Please refer to 'ncdfFlow' package for cdf based
#' storage.
#' @param dataset see \code{\link[flowCore]{read.FCS}} for details.
#' @param min.limit see \code{\link[flowCore]{read.FCS}} for details.
#' @param truncate.max.range see \code{link[flowCore]{read.FCS}} for details.
#' @param emptyValue see \code{\link[flowCore]{read.FCS}} for details.
#' @param ignore.text.offset see \code{\link[flowCore]{read.FCS}} for details.
#' @param truncate_max_range see \code{\link[flowCore]{read.FCS}} for details.
#' @param channel_alias see \code{\link[flowCore]{read.FCS}} for details.
#' @param \dots Further arguments that get passed on to
#' \code{\link[Biobase]{read.AnnotatedDataFrame}}, see details.
#' 
#' @return An object of class \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @author F. Hahne, N.Le Meur, B. Ellis
#' @keywords IO
#' @examples
#' 
#' fcs.loc <- system.file("extdata",package="flowCore")
#' file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")
#' 
#' samp <- read.flowSet(file.location[1:3])
#' 
#' 
#' @export
read.flowSet <- function(files=NULL, path=".", pattern=NULL, phenoData,
                         descriptions, name.keyword, alter.names=FALSE,
                         transformation="linearize", which.lines=NULL,
                         column.pattern=NULL, invert.pattern = FALSE, decades=0,
                         sep="\t", as.is=TRUE, name, ncdf=FALSE, dataset=NULL,
                         min.limit=NULL, truncate_max_range = TRUE, emptyValue=TRUE
                         , ignore.text.offset = FALSE
                         , channel_alias = NULL
                         , ...)
{
   channel_alias <- check_channel_alias(channel_alias)
    if(ncdf)
      .Deprecated("'ncdf' argument is deprecated!Please use 'ncdfFlow' package for hdf5-based data structure.")
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
            fnams <- grep("file|filename", varLabels(phenoData),
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
            phenoFrame$name <- basename(files)
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
                      column.pattern=column.pattern, invert.pattern = invert.pattern,
                      decades=decades,min.limit=min.limit,emptyValue=emptyValue
                      , ignore.text.offset = ignore.text.offset
                      , dataset=dataset
                      , truncate_max_range = truncate_max_range
                      , channel_alias = channel_alias)
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


#' @export
cleanup <- function() if(file.exists(".flowCoreNcdf"))

     .Deprecated("'ncdf' is deprecated!Please use 'ncdfFlow' package for hdf5-based data structure.")
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
        writeChar(sprintf("%8s", val1), con, eos=NULL)
        writeChar(sprintf("%8s", val2), con, eos=NULL)
    }
     invisible()
}



#' collapse and concatenate the content of the description slot into a single character string
#' @param x flowFrame
#' @examples 
#' fr <- GvHD[[1]]
#' collapseDesc(fr)
#' @noRd 
collapseDesc <- function(x, delimiter = "\\")
{
    d <- description(x)
	##make sure there is no empty value for each keyword in order to conform to FCS3.0
#	browser()
  ##skip flowCore_$PnRMax when trans is not applied 
  trans <- d[["transformation"]]
  if(is.null(trans))
  {
    for(i in 1:ncol(x))
    {
      d[[paste0("flowCore_$P", i, "Rmin")]] <- NULL
      d[[paste0("flowCore_$P", i, "Rmax")]] <- NULL
    }
  }
  d <- collapse_desc(d)
  d <- lapply(d, function(y){
        double_delimiter <- paste0(delimiter, delimiter)
        
        #now  we escape every single delimiter character occurances(include multi-delimiter string) by doubling it
        gsub(delimiter, double_delimiter, y, fixed = TRUE, useBytes = TRUE)
    
  })
  paste(delimiter, iconv(paste(names(d), delimiter, sapply(d, paste, collapse=" "),
						  delimiter, collapse="", sep=""), to="latin1",
				  sub=" "), sep="")
  
}

#' Coerce the list of the keywords into a character 
#' Also flatten spillover matrix into a string
#' 
#' @param d a named list of keywords
#' @param collapse.spill whether to flatten spillover matrix to a string
#' @return a list of strings
#' @examples
#' data(GvHD)
#' fr <- GvHD[[1]]
#' collapse_desc(keyword(fr))
#' @export 
collapse_desc <- function(d, collapse.spill = TRUE)
{
	d <- lapply(d, function(y){
				if(length(y)==0)
					return(" ")
				else
				{
					#make sure spillover matrix doesn't get converted to vector
					if(is.matrix(y))
						return (y)
					else{

					  y <- sub("^$"," ",y)
					}

				}

  
			})
    d <- d[order(names(d))]

  if(collapse.spill)
  {
      
    for(spillName in .spillover_pattern)
    {
       mat <-  d[[spillName]]
       if(!is.null(mat))
       {
    		rNum <- as.character(nrow(mat))
    		clNames <- paste(colnames(mat),sep=",")
    		vec <- paste(c(t(mat)),sep=",",collapse=",")
  	    d[spillName] <- paste(c(rNum,clNames,vec),sep=",",collapse=",")
       }
    }
    
  }
    d
}



## ==========================================================================
## Write a flowFrame to an FCS file. 'what' controls the output
## data type.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' Write an FCS file
#' 
#' Write FCS file from a flowFrame
#' 
#' 
#' The function \code{write.FCS} creates FCS 3.0 standard file from an object
#' of class \code{flowFrame}.
#' 
#' For specifications of FCS 3.0 see \url{http://www.isac-net.org} and the file
#' \url{../doc/fcs3.html} in the \code{doc} directory of the package.
#' 
#' @name write.FCS
#' @aliases write.FCS
#' 
#' @usage 
#' write.FCS(x, filename, what="numeric", delimiter = "|", endian="big")
#' 
#' @param x A \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param filename A character scalar giving the output file name.
#' @param what A character scalar defining the output data type. One in
#' \code{integer}, \code{numeric}, \code{double}. Note that forcing the data
#' type to \code{integer} may result in considerable loss of precision if the
#' data has been transformed. We recommend using the default data type unless
#' disc space is an issue.
#' @param delimiter a single character to separate the FCS keyword/value pairs.
#' Default is : "|"
#' @param endian a character, either "little" or "big" (default), specifying
#' the most significant or least significant byte is stored first in a 32 bit
#' word.
#' 
#' @return
#' 
#' A character vector of the file name.
#' @author F. Hahne
#' @seealso \code{link[flowCore]{write.flowSet}}
#' @keywords IO
#' @examples
#' 
#' ## a sample file
#' inFile <- system.file("extdata", "0877408774.B08", package="flowCore")
#' foo <- read.FCS(inFile, transform=FALSE)
#' outFile <- file.path(tempdir(), "foo.fcs")
#' 
#' ## now write out into a file
#' write.FCS(foo, outFile)
#' bar <- read.FCS(outFile, transform=FALSE)
#' all(exprs(foo) == exprs(bar))
#' 
#' 
#' @export
write.FCS <- function(x, filename, what="numeric", delimiter = "|", endian = "big")
{
  # warning("'write.FCS' is not fully tested and should be considered as experimental.")
    ## Some sanity checking up front
    if(missing(filename))
    {
      stop("output filename needs to be provided!")
#        filename <- identifier(x)
#        if(!length(grep(".", filename, fixed=TRUE)))
#            filename <- paste(filename, "fcs", sep=".")
    }
    what <- match.arg(what, c("integer", "numeric", "double"))
    endian <- match.arg(endian, c("little","big"))
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
    # endian <- "big"
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
    npar <- ncol(x)
    pd <- pData(parameters(x))
    pid <- rownames(pd)
    newid <- seq_len(npar)
    pnb <- as.list(rep(types[what, "bitwidth"]*8, ncol(x)))
    names(pnb) <- paste0("$P", newid, "B")
    mk <- c(mk, pnb)
#	browser()
    orig.kw <- description(x)
    ## We need all PnE keywords and assume "0,0" if they are missing
    pne <- orig.kw[paste0(pid, "E")]
    names(pne) <- sprintf("$P%sE", newid)
    pne[sapply(pne, length)==0] <- "0,0"
    mk <- c(mk, pne)
    
    mat <- exprs(x)
    ## The same for PnR, 
    pnr <- orig.kw[paste0(pid, "R")]
    names(pnr) <- sprintf("$P%sR", newid)
    empty_range_ind <- sapply(pnr, length)==0
    empty_range_chnl <- as.character(pd[pid[empty_range_ind],"name"])
    pnr[empty_range_ind] <- colMaxs(mat[, empty_range_chnl, drop = FALSE])
    mk <- c(mk, pnr)
    
    ## Now update the PnN keyword
    pnn <- colnames(x)
    names(pnn) <- sprintf("$P%sN", newid)
    mk <- c(mk, pnn)
    
    ## Now update the PnS keyword
    pns <- as.vector(pd[["desc"]])
    newid.pns <- newid[!is.na(pns)]
    pns <- pns[!is.na(pns)]
    names(pns) <- sprintf("$P%sS", newid.pns)
    mk <- c(mk, pns)
    
    # Must correct PnX keys for the case that fr was subsetted by bumping them down to be consecutive (like newid)
    names(newid) <-  gsub("\\$", "", pid)
    # bump indices on remaining keys down to their new values
    newnames <- names(x@description)
    for(i in newid){
      newnames <- gsub(paste0(names(newid)[[i]], "([a-zA-z])"), paste0("P", i, "\\1"), newnames)
    }
    names(x@description) <- newnames
    description(x) <- mk
    
    ## Figure out the offsets based on the size of the initial text section
    ld <-  length(mat) * types[what, "bitwidth"]
    ctxt <- collapseDesc(x, delimiter = delimiter)
    
    #try to update the Txt with actual BEGINDATA and  ENDDATA values
    kw.len.old <- 2 #two initial '0' s
    endTxt <- nchar(ctxt, "bytes") + begTxt - 1#must use type = "bytes" as txt could contain special character(e.g. German umlauts) that results in multi-byte write by writeChar
    repeat#keep updating txt until its length stop increasing 
    {
      #compute the offsets based on endTxt
      datastart <- endTxt + 1  
      dataend <-  datastart + ld - 1 
      
      #check if new values cause txt len to increase
      kw.len.new <- nchar(datastart) + nchar(dataend)
      if(kw.len.new > kw.len.old)
      {
        #update the endTxt
        endTxt <- endTxt + kw.len.new - kw.len.old
        kw.len.old <- kw.len.new
        
      }else
        break
    }
  
    description(x) <- list("$BEGINDATA" = datastart,
                           "$ENDDATA" = dataend)
    ctxt <- collapseDesc(x, delimiter = delimiter)

    offsets <- c(begTxt, endTxt, datastart, dataend, 0,0)
    ## Write out to file
    con <- file(filename, open="wb")
    on.exit(close(con))
    writeFCSheader(con, offsets)
    writeChar(ctxt, con, eos=NULL)
    
    # Break writes of data in to 2^30 bytes, to stay under 
    # writeBin limit of 2^31-1
    bitwidth <- types[what, "bitwidth"]
    chunk_size <- (2^30)%/%bitwidth
    mat <- as(t(mat), what)
    n_chunks <- ceiling(length(mat)/chunk_size)
    for (idx in seq_len(n_chunks)){
      writeBin(mat[(((idx-1) * chunk_size) + 1):
                   min(idx * chunk_size, length(mat))],
                   con, size=bitwidth, endian=endian)
    }

    
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
#' Write an FCS file
#' 
#' Write FCS file for each flowFrame in a flowSet
#' 
#' 
#' The function \code{write.flowSet} creates FCS 3.0 standard file for all
#' \code{flowFrames} in an object of class \code{flowSet}. In addition, it will
#' write the content of the \code{phenoData} slot in the ASCII file
#' \code{"annotation.txt"}. This file can subsequently be used to reconstruct
#' the whole \code{flowSet} using the \code{\link[flowCore]{read.flowSet}}
#' function, e.g.:
#' 
#' \code{read.flowSet(path=outdir, phenoData="annotation.txt"}
#' 
#' The function uses \code{\link[flowCore]{write.FCS}} for the actual writing
#' of the FCS files.
#' 
#' @name write.flowSet
#' @aliases write.flowSet
#' 
#' @usage 
#' write.flowSet(x, outdir=identifier(x), filename, \dots)
#' 
#' @param x A \code{\link[flowCore:flowSet-class]{flowSet}}.
#' @param filename A character scalar or vector giving the output file names.
#' By default, the function will use the identifiers of the individual
#' \code{flowFrames} as the file name, potentially adding the \code{.fcs}
#' suffix unless a file extension is already present. Alternatively, one can
#' supply either a character scalar, in which case the prefix \code{i_} is
#' appended (\code{i} being an integer in \code{seq_len(length(x))}), or a
#' character vector of the same length as the \code{flowSet x}.
#' @param outdir A character scalar giving the output directory. As the
#' default, the output of \code{identifier(x)} is used.
#' @param \dots Further arguments that are passed on to
#' \code{\link[flowCore]{write.FCS}}.
#' @return
#' 
#' A character vector of the output directory.
#' @author F. Hahne
#' @seealso \code{link[flowCore]{write.FCS}}
#' @keywords IO
#' @examples
#' 
#' ## sample data
#' data(GvHD)
#' foo <- GvHD[1:5]
#' outDir <- file.path(tempdir(), "foo")
#' 
#' ## now write out into  files
#' write.flowSet(foo, outDir)
#' dir(outDir)
#' 
#' ## and read back in
#' bar <- read.flowSet(path=outDir, phenoData="annotation.txt")
#' 
#' 
#' @export
write.flowSet <- function(x, outdir=identifier(x), filename, ...)
{
    checkClass(x, "flowSet")
    nSample <- length(x)
    ## Some sanity  checking up front
    if(missing(filename))
    {
        filename <- as.character(fsApply(x, identifier))
    }
    else
    {
        nFilename <- length(filename)
        ferr <- paste("Argument filename has to be a character scalar or",
                      "a character vector of the same length as 'x'.")
        if(is.character(filename))
        {
            if(nFilename==1&&nSample>1)
                filename <- paste(seq_len(nSample), filename, sep="_")
            else if(nFilename != nSample)
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
    for(f in seq_len(nSample))
    {
        if(!length(grep(".", filename[f], fixed=TRUE)))
            filename[f] <- paste(filename[f], "fcs", sep=".")
        write.FCS(x[[f]], filename=file.path(outdir, filename[f]), ...)
    }
    sampleNames(x) <- filename#force the sampleNames (thus rownames of Pdata) to be identical to the fcs filenames
    pData(x)$FCS_File <- filename
    write.AnnotatedDataFrame(phenoData(x), file=file.path(outdir, "annotation.txt"))
    outdir
}
