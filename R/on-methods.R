## ==========================================================================
## %on% methods are basically constructors for transformFilter or
## transformList objects. They allow for an inline-syntax for transformation
## operations.
## ==========================================================================





#' Methods for Function \%on\% in Package `flowCore'
#' 
#' This operator is used to construct a \code{transformFilter} that first
#' applies a \code{transformList} to the data before applying the \code{filter}
#' operation. You may also apply the operator to a \code{flowFrame} or
#' \code{flowSet} to obtain transformed values specified in the list.
#' 
#' 
#' @name filter-on-methods
#' @aliases filter-on-methods %on%-methods %on%
#' %on%,filter,transformList-method %on%,filter,transform-method
#' %on%,filter,parameterTransform-method
#' %on%,parameterTransform,flowFrame-method %on%,transform,flowFrame-method
#' %on%,transformList,flowFrame-method %on%,transformList,flowSet-method
#' %on%,ANY,flowSet-method
#' 
#' @param e1 a \code{\linkS4class{filter}}, \code{\linkS4class{transform}},
#' or \code{\linkS4class{transformList}} object
#' @param e2 a \code{\linkS4class{transform}}, \code{\linkS4class{transformList}},
#' \code{\linkS4class{flowFrame}}, or \code{\linkS4class{flowSet}} object
#' 
#' @usage 
#' e1 \%on\% e2
#' 
#' @docType methods
#' @author B. Ellis
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
#' plot(transform("FSC-H"=log, "SSC-H"=log) %on% samp)
#' 
#' 
#' rectangle <- rectangleGate(filterId="rectangleGateI","FSC-H"=c(4.5, 5.5))
#' sampFiltered <- filter(samp, rectangle %on% transform("FSC-H"=log, "SSC-H"=log))
#' res <- Subset(samp, sampFiltered)
#' 
#' plot(transform("FSC-H"=log, "SSC-H"=log) %on% res)
#' 
#' 
## ===========================================================================
## Constructor for a transformFilter for a transform
## ---------------------------------------------------------------------------
#' @export
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="transform"),
          definition=function(e1, e2)
      {
          new("transformFilter",
              filterId=paste(identifier(e1), "on transformed values"),
              transforms=parameterTransform(e2, parameters(e1)),
              filter=e1, parameters=parameters(e1))
      })



## ===========================================================================
## Wrapper for a flowFrame. We transform on the actual raw data
## ---------------------------------------------------------------------------
#' @export
setMethod("%on%",
          signature=signature(e1="transform",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          exprs(e2) <- e1(exprs(e2))
          e2
      })



## ===========================================================================
## Constructor for a transformFilter from a parameterTransform
## ---------------------------------------------------------------------------
#' @export
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="parameterTransform"),
          definition=function(e1, e2)
      {
          new("transformFilter",
              filterId=paste(e1@filterId," on transformed values of ",
              paste(parameters(e2), sep=","),collapse=" "),
              transforms=e2, filter=e1,
              parameters=unique(c(parameters(e1), parameters(e1))))
})



## ===========================================================================
## Wrapper for a flowFrame. We transform on the actual raw data
## ---------------------------------------------------------------------------
#' @export
setMethod("%on%",
          signature=signature(e1="parameterTransform",
                              e2="flowFrame"),
          definition=function(e1, e2)
          {
            params <- if(length(parameters(e1))) colnames(e2) else
            parameters(e1)
            exprs(e2)[,params] <- e1(exprs(e2))[,params]
            e2
          })



## ===========================================================================
## Constructor for a transformFilter from a transformList
## ---------------------------------------------------------------------------
#' @export
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="transformList"),
          definition=function(e1,e2)
      {
          new("transformFilter",
              filterId=paste(e1@filterId, "on transformed values of",
              paste(colnames(e2), collapse=",")),
              transforms=e2,
              filter=e1)
      })



## ==========================================================================
## This is the real workhorse for doing %on% type transformations.
## Eventually, all %on% methods are going to call this.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("%on%",
          signature=signature(e1="transformList",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          x <- exprs(e2)
          pars <- parameters(e2)
          ranges <- range(e2)
          cN <- colnames(x)
          for(y in e1@transforms){
              if( !(y@output %in% cN) )
                  stop(y@output, " is not a variable in the flowFrame")
              x[,y@output] = y@f(x[, y@input])
              mt <- match(y@output, pars$name)
              pData(pars)[mt,c("minRange", "maxRange")] <-
                  y@f(ranges[,y@input])   
			 
          }
          exprs(e2) <- x
          parameters(e2) <- pars
          kw <- updateTransformKeywords(e2)	
		  keyword(e2)[names(kw)] <- kw 
          e2
      })



## ==========================================================================
## General %on% implementation for a flowSet. Basically a wrapper around
## fsApply
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("%on%",
          signature=signature(e1="ANY",
                              e2="flowSet"),
          definition=function(e1, e2) fsApply(e2, "%on%", e1=e1))
