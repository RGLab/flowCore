## ==========================================================================
## fcReferences provide a means to reference into workFlow environments.
## We provide methods to get and assign object from references, as well
## as methods to remove them along with the referenced objects
## ==========================================================================






## ==========================================================================
## Remove an object referenced to by an fcReference object as well as the
## reference object itself (if possible)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Just the object and the reference
setMethod("Rm",
          signature=signature(symbol="fcReference",
                              envir="missing",
                              subSymbol="character"),
          definition=function(symbol, subSymbol, rmRef=TRUE)
      {
          rm(list=identifier(symbol), envir=symbol@env)
          ## This doesn't work because the name space is messing up the
          ## evaluation environment tree.
          ## rm(list=subSymbol, inherits=TRUE)
          ## Instead we are rather bold and remove the reference from the
          ## global env if it exists.
          if(subSymbol %in% ls(globalenv(), all.names=TRUE) &&
             is(get(subSymbol), globalenv()))
              rm(list=subSymbol, envir=globalenv())
      })

## Object and reference and all dependent objects (further down in the
## workFlow tree). This is only recursive for actionItems and views. All
## other objects can be removed without further side effects, but we
## strongly discourage from doing so since the consequences might be
## catastrophic...
setMethod("Rm",
          signature=signature(symbol="fcReference",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol, rmRef=TRUE)
      {
          ##callGeneric(symbol)
          selectMethod("Rm", signature("fcReference", "missing", "character"))(symbol,
                                                                               identifier(symbol))
      })


## For nullReferences we don't need to do anything
setMethod("Rm",
          signature=signature(symbol="fcNullReference",
                              envir="missing",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol, ...){})





## ==========================================================================
## Test for a NULL reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("isNull",
          signature=signature("fcReference"),
          definition=function(f) is(f, "fcNullReference"))



## ==========================================================================



## ==========================================================================
## Resolve a reference, i.e., get the symbol 'ID' from the environment
## 'env'. Resolving NULL references always returns NULL
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## for NULL references
setMethod("get",
          signature=signature(x="fcNullReference", pos = "missing",
          envir = "missing",mode = "missing", inherits = "missing"),
          definition=function(x) NULL)

## for all other references
setMethod("get",
          signature=signature(x="fcReference", pos = "missing",
          envir = "missing",mode = "missing", inherits = "missing"),
          definition=function(x)
          {
              if(!resolved(x, FALSE)){
                  mess <- paste("Unable to resolve reference to object '",
                                x@ID, "'", sep="")
                  if(!is.na(w <- pmatch(tolower(x@ID), 
                                        tolower(ls(x@env))))) 
                      mess <- paste(mess, sprintf("Perhaps you meant %s?",
                                sQuote(ls(x@env)[w])), sep="\n")
                  stop(mess, call.=FALSE)
              }
              get(x@ID, x@env)
          })



## ==========================================================================
## Check whether a reference can be resolved (i.e., whether the referenced
## object exists in the evaluation environment)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resolved <- function(object, verbose=TRUE)
{
    unres <- nchar(object@ID) && !object@ID %in% ls(object@env)
    if(unres && verbose)
        cat("\nUnable to resolve reference to object '",
              object@ID, "'.", sep="")
    return(!unres)
}

