## Adapted from F.Hahne last updated Oct 30 0.6
## For specifications of FACS 3.0 see
## http://www.isac-net.org and the file
## fcs3.html in the doc directory

read.FCS <- function(filename,transformation=NULL,debug=FALSE,alter.names=FALSE)
{
  stopifnot(is.character(filename), length(filename)==1, filename!="")
  con <- file(filename, open="rb")

  if((is.logical(transformation) && transformation) || transformation == "linearize") {
	transformation = TRUE
	scale = FALSE
  } else if(transformation == "scale") {
	transformation = TRUE
	scale = TRUE
  }
  
  offsets <- readFCSheader(con)
  txt     <- readFCStext(con, offsets, debug)
  mat     <- readFCSdata(con, offsets, txt, transformation, debug, scale, alter.names)

  close(con)
  
  if(as.integer(readFCSgetPar(txt, "$TOT"))!=nrow(mat))
    stop(paste("file", filename, "seems to corrupted."))
  if(transformation==TRUE){
      txt[["tranformation"]] <-"applied" 
  }
  return(new("flowFrame", exprs=mat, description=c(txt)))
}

readFCSgetPar <- function(x, pnam) {
  stopifnot(is.character(x), is.character(pnam)) 
  i <- match(pnam, names(x))
  if(any(is.na(i)))
    stop(paste("Parameter(s)", pnam, "not contained in 'x'"))
  return(x[i])
}

readFCSheader <- function(con) {
  seek(con, 0)
  version <- readChar(con, 6)
  if(!version %in% c("FCS2.0", "FCS3.0"))
    stop("This does not seem to be a valid FCS2.0 or FCS3.0 file")
  
  tmp <- readChar(con, 4)
  stopifnot(tmp=="    ")

  coffs <- character(6)
  for(i in 1:length(coffs))
    coffs[i] <- readChar(con=con, nchars=8)
  
  ioffs <- as.integer(coffs)
  names(ioffs) <- c("textstart", "textend", "datastart", "dataend", "anastart", "anaend")
  
  stopifnot(all(!is.na(ioffs) | coffs=="        "), !any(ioffs[1:4]==0))
  return(ioffs)
}

readFCStext <- function(con, offsets, debug) {
  seek(con, offsets["textstart"])
  txt <- readChar(con, offsets["textend"]-offsets["textstart"]+1)
  txt <- iconv(txt, "", "latin1", sub="byte")
  delimiter <- substr(txt, 1, 1)
  sp  <- strsplit(substr(txt, 2, nchar(txt)), split=delimiter, fixed=TRUE)[[1]]
  ## if(length(sp)%%2!=0)
  ##  stop("In readFCStext: unexpected format of the text segment")
  rv <- sp[seq(2, length(sp), by=2)]
  names(rv) <- sp[seq(1, length(sp)-1, by=2)]
  return(rv)
}

readFCSdata <- function(con, offsets, x, transformation, debug, scale, alter.names) {
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
  range    <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "R", sep="")))
  bitwidth <- as.integer(readFCSgetPar(x, paste("$P", 1:nrpar, "B", sep="")))
  bitwidth <- unique(bitwidth)
  if(length(bitwidth)!=1)
    stop("Sorry, I am expecting the bitwidth to be the same for all parameters")

  seek(con, offsets["datastart"])

  size <- bitwidth/8
  if (!size %in% c(1, 2, 4, 8))
    stop(paste("Don't know how to deal with bitwidth", bitwidth))

  dat <- readBin(con, dattype, n = (offsets["dataend"]-offsets["datastart"]+1)/size,
                 size=size, signed=FALSE, endian=endian)

  stopifnot(length(dat)%%nrpar==0)
  dat <- matrix(dat, ncol=nrpar, byrow=TRUE)
  colnames(dat) <- if(alter.names) 
	make.names(readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))) 
  else 
	readFCSgetPar(x, paste("$P", 1:nrpar, "N", sep=""))

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
        ampliPar <- readFCSgetPar(x, paste("$P", 1:nrpar, "E", sep=""))
        ampli <- do.call("rbind",lapply(ampliPar,function(x) as.integer(unlist(strsplit(x,",")))))		
		for(i in 1:nrpar) {
			dat[,i] = if(ampli[i,1] > 0) dat[,i]/(10^ampli[i,1]) else dat[,i]/(range[i]-1)
		}
	}
  return(dat) 
}



