## ==========================================================================
## Accessor to parentView slot. This returns the view object after
## resolving the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parent",
          signature=signature("actionItem"),
          definition=function(object) get(object@parentView))



## ==========================================================================
## Accessor to name slot.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature("actionItem"),
          definition=function(x) x@name)



## ==========================================================================
## Accessor to gate slot. This returns the filter object after resolving
## the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("gate",
          signature=signature("gateActionItem"),
          definition=function(object) get(object@gate))



## ==========================================================================
## Remove an actionItem object from a workFlow. This will traverse down the tree
## and also remove all dependent objects to free their memory.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## For all actionItems we need to remove the associated views.
setMethod("Rm",
          signature=signature(symbol="actionItem",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol)
      {
          curAction <- identifier(symbol)
          aiMatch <- sapply(edgeData(get(envir@tree), attr="actionItem"), identifier)
          depNodes <- gsub(".*\\|", "", names(aiMatch[aiMatch == curAction]))
          lapply(mget(depNodes, envir), Rm, envir)
          return(invisible(NULL))
      })

## For gateActionItems we need to remove the gate and the filterResult
setMethod("Rm",
          signature=signature(symbol="gateActionItem",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol)
      {
          Rm(symbol@gate)
          selectMethod("Rm", signature("actionItem",
                                       "workFlow", "character"))(symbol, envir)
          rm(list=identifier(symbol), envir=envir@env)
          return(invisible(NULL))         
      })

## For transformActionItems we need to remove the transformation object
setMethod("Rm",
          signature=signature(symbol="transformActionItem",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol)
      {
          Rm(symbol@transform)
          selectMethod("Rm", signature("actionItem",
                                       "workFlow", "character"))(symbol, envir)
          rm(list=identifier(symbol), envir=envir@env)
          return(invisible(NULL))
      })

## For compensateActionItems we need to remove the compensation object
setMethod("Rm",
          signature=signature(symbol="compensateActionItem",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol)
      {
          Rm(symbol@compensate)
          selectMethod("Rm", signature("actionItem",
                                       "workFlow", "character"))(symbol, envir)
          rm(list=identifier(symbol), envir=envir@env)
          return(invisible(NULL))
      })
