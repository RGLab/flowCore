## ==========================================================================
## Objects of class view provide wrappers for the results of common flow
## operations in order to organize them in a workflow. There are three
## subclasses: gateView, transformView and compensateView.
## ==========================================================================






## ==========================================================================
## Accessor to action slot. This returns the actionItem object after
## resolving the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("action",
          signature=signature(object="view"),
          definition=function(object) get(object@action))


## ==========================================================================
## Accessor to data slot. This returns the flow data object (flowFrame or
## flowSet) after resolving the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("Data",
          signature=signature(object="view"),
          definition=function(object){
              dat <- get(object@data)
              if(is.null(dat)){
                  parent <- parent(object)
                  ## The parent gate is not be evaluated yet, we do that now
                  if(is(object, "gateView")){
                      parentData <- Data(parent)
                      if(is(parentData, "flowFrame")){
                          newData <- parentData[object@indices[[1]],]
                      }else{
                          newData <- parentData[1:length(parentData)]
                          for(i in seq_along(object@indices))
                              newData[[i]] <- newData[[i]][object@indices[[i]],]
                      }
                      newDataID <- refName(newData)
                      env <- object@env
                      ## We also need to update the appropriate journal entry
                      jid <- grep("journalRef", ls(env), value=TRUE)
                      if(length(jid) != 1)
                          stop("Unable to find journal in this environment.",
                               call.=FALSE)
                      journal <- journal(env)
                      assign(x=newDataID, value=newData, envir=env)
                      actionId <- identifier(action(object))
                      if(actionId %in% names(journal)){
                          journal[[actionId]] <- c(journal[[actionId]],
                                                          newDataID)
                          assign(jid, journal, env)
                      }
                      object@data <- new("fcDataReference", ID=newDataID,
                                         env=env)
                      assign(identifier(object), value=object, envir=env)
                      return(newData)
                  }else{
                      stop("Unable to find data for this view.\n This workFlow object ",
                           "is corrupted beyond repair.", call.=FALSE)
                  }
              }
              return(dat)
          })



## ==========================================================================
## Accessor to parent, linked by an actionItem returns. This the view object
## after resolving the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parent",
          signature=signature(object="view"),
          definition=function(object) parent(action(object)))

setMethod("parent",
          signature=signature(object="NULL"),
          definition=function(object) NULL)


## ==========================================================================
## Get the alias table from a view
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("alias",
          signature=signature(object="view"),
          definition=function(object) get(object@alias))
