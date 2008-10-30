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
          definition=function(object) get(object@data))



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
