setMethod("%in%",signature("flowFrame","expressionFilter"),function(x,table) {
	data = as.data.frame(exprs(x))
        if(length(table@args)>1){
            res <- eval(table@expr,data,env=table@args)
        }else{
            res <- eval(table@expr,data)
        }
})

          
## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("expressionFilter"),function(object) {
	msg = paste("An expression filter named",identifier(object),":",deparse(object@expr),sep=" ")
	cat(msg)
	cat("\n")
	invisible(msg)
})
