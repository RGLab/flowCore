setMethod("summary",signature("logicalFilterResult"),function(object,...) summary(filterDetails(object,1)$filter,object))
