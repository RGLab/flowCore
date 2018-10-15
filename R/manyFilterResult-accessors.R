## ==========================================================================
## manyFilterResults create a multiple potentially overlapping populations
## ==========================================================================






## ==========================================================================
## Subsetting to the individual populations. We create a new
## logicalFilterResult.
## ---------------------------------------------------------------------------
#' @export
setMethod("[[",
          signature=signature(x="manyFilterResult"),
          definition=function(x,i,j,drop=FALSE)
      {
          if(is.numeric(i)) i <- names(x)[i]
          if(length(i)!=1)
              stop("Only a single subpopulation can be selected.")
          if(!i %in% names(x))
              stop("Index out of bounds.")
          filterDetails <- x@filterDetails
          names(filterDetails) <- i
          filterDetails$population <- i
          filterDetails$source <- identifier(x)
          dat <- if(i == "rest") !apply(x@subSet, 1, any) else x@subSet[,i]
          new("logicalFilterResult", subSet=dat,
              filterDetails=filterDetails,
              frameId=x@frameId, filterId=i)
      })



## ==========================================================================
## convert to a data frame
## ---------------------------------------------------------------------------
#' @export
as.data.frame.manyFilterResult <- function(x, row.names=NULL,
                                           optional=FALSE,...)
{
    nrows <- length(x)
    nm <- if(nzchar(identifier(x))) identifier(x) else "filter"
    if(is.null(row.names)) {
        if(nrows == 0L)
            row.names <- character(0L)
        else if(length(row.names <- names(x)) == nrows && !any(duplicated(row.names))) {
        } else row.names <- .set_row_names(nrows)
	}
    values <- sapply(names(x),function(i) {
        m <- x[[i]]
        paste(deparse(as(filterDetails(m,identifier(m))$filter,"call")),sep="\n",collapse="\n")
    })
    names(values) <- NULL
    values <- list(values)
    if(!optional)
        names(values) <- nm
    attr(values,"row.names") <- row.names
    class(values) <- "data.frame"
    values
}
