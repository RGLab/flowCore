## ==========================================================================
## filter flowFrame object using ellipsoidGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## FIXME: ELLIPOID GATES HAVE TO BE DEFINED BY VECTOR OF MEANS AND
## A COVARIANCE MATRIX

setMethod("%in%",signature(x="flowFrame",table="ellipsoidGate"),function(x,table) {
  e <- if(length(table@parameters)==1) as.matrix(exprs(x)[,table@parameters]) else
          exprs(x)[,table@parameters]
  d1 <- apply(sapply(seq(along=table@parameters), function(i) (e[,i] -
                           table@focus[1,i])^2), 1, function(k) sqrt(sum(k)))
  d2 <- apply(sapply(seq(along=table@parameters), function(i) (e[,i] -
                           table@focus[2,i])^2), 1, function(k) sqrt(sum(k)))
  (d1+d2) <= table@distance
})
