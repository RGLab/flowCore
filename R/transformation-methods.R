## NLM Jan 17
## ==========================================================================
## Transformation function for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...) {
			e <- substitute(list(...))
                        transformed <- as.matrix(transform(as.data.frame(exprs(`_data`)),...))
                        ##Add any new parameter values
			param.names <- colnames(transformed)
                        newParams <- is.na(match(param.names,`_data`@parameters$name))
			#params      = parameters(`_data`)
                        params <- `_data`@parameters$name
                        if(any(newParams)) {
				params <- cbind(params,
                                  data.frame(name=param.names))
                                ##params = cbind(params,
                                ##data.frame(name=param.names[newParams]))
                            }
                        colnames(transformed) <- `_data`@parameters$name
                        new("flowFrame",
				##exprs=transformed[,params$name],
                                exprs=transformed,
				parameters=`_data`@parameters,
				description=description(`_data`))
          })
 
