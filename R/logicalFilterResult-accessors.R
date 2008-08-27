setMethod("summary",
          signature=signature("logicalFilterResult"),
          definition=function(object,...)
          summary(filterDetails(object,1)$filter,object))


setMethod("names",
          signature=signature("logicalFilterResult"),
          definition=function(x) paste(x@filterId, c("+","-"), sep=""))


## ==========================================================================
## length method, how many populations do we have?
## ---------------------------------------------------------------------------
setMethod("length","logicalFilterResult",function(x) 2)
