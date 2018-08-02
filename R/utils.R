.spillover_pattern <- c("SPILL", "spillover", "$SPILLOVER")

#' fuzzy match of marker/channel names
#' @param pd pData of parameters of flowFrame
#' @param name \code{character} the string to match
#' @param fix whether to do regexpr match
#' @param parital whether to do the complete word match or parital match
#' @noRd
.flowParamMatch <- function(pd, name, fix = FALSE, partial = FALSE) {
  
  # try to compelete word match by following with a space or the end of string
  if (partial)
    pname <- name
  else
    pname <- paste0(name, "([ ]|$)")
  
  if (fix) {
    ind <- which(toupper(pd$name) %in% toupper(name))
  } else {
    ind <- which(grepl(pname, pd$name, ignore.case = T))
  }
  
  if (length(ind) == 0) {
    # try marker name
    ind <- which(unlist(lapply(pd$desc, function(x) {
                  # split by white space and then match each individual string
                  if (fix) {
                    any(unlist(lapply(strsplit(x, " ", fixed=TRUE), function(y) toupper(y) %in% toupper(name))))
                  } else {
                    grepl(pattern = pname, x, ignore.case = T)
                  }
                })))
  }
  ind
}
#' get channel and marker information from a \code{flowFrame} that matches to the given keyword
#'
#' This function tries best to guess the flow parameter based on the keyword supplied by \code{name}
#' It first does a complete word match(case insensitive) between \code{name} and flow channels and markers.
#' If there are duplcated matches, throw the error. If no matches, it will try the partial match.
#'
#' @return
#' an one-row \code{data.frame} that contains "name"(i.e. channel) and "desc"(i.e. stained marker) columns.
#'
#' @param frm \code{flowFrame} object
#' @param name \code{character} the keyword to match
#' @param ... other arguments: not used.
#' @export
getChannelMarker <- function(frm, name, ...) {
  #escape ( since we see that often times in Cytof data
  
  name <- gsub(")", "\\)", gsub("(", "\\(", name, fixed = TRUE), fixed = TRUE)
  pd <- pData(parameters(frm))
  # try complete match first
  ind <- .flowParamMatch(pd, name, ...)
  
  if (length(ind) > 1) {
    stop("multiple markers matched: ", name)
  }
  
  if (length(ind) == 0) {
    # if no match then give a second try to patial match
    ind <- .flowParamMatch(pd, name, partial = TRUE, ...)
    if (length(ind) == 0)
      stop("can't find ", name) else if (length(ind) > 1)
      stop("multiple markers matched: ", name) else warning(name, " is partially matched with ", pd[ind, c("name", "desc")])
  }
  
  pd <- pd[ind, c("name", "desc")]
  pd[, "name"] <- as.vector(pd[, "name"])
  pd[, "desc"] <- as.vector(pd[, "desc"])
  pd
  
}

LdFlags <- function(){
  libpath <- paste0("lib", Sys.getenv("R_ARCH"), "/libboost_regex.a")
  cat(tools::file_path_as_absolute( base::system.file(libpath, package = "flowCore" )))
}
