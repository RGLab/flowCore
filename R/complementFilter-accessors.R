## ==========================================================================
## Returns TRUE when the input filter is FALSE.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","complementFilter"),
          function(x,table) !(x%in%table@filters[[1]]))

