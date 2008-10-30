## ==========================================================================
## Objects of class actionItem provide wrappers for the common flow
## operations in order to organize them in a workflow. There are three
## subclasses: gateActionItem, transformActionItem and compensateActionItem.
## ==========================================================================

## ==========================================================================
## Accessor to parentView slot. This returns the view object after
## resolving the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parent",
          signature=signature(object="actionItem"),
          definition=function(object) get(object@parentView))



## ==========================================================================
## Accessor to gate slot. This returns the filter object after resolving
## the reference
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("gate",
          signature=signature(object="gateActionItem"),
          definition=function(object) get(object@gate))



## ==========================================================================
## Get the alias table from an actionItem
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("alias",
          signature=signature(object="actionItem"),
          definition=function(object) get(object@alias))
