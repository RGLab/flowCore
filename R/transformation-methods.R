## ==========================================================================
## Transformation function for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...) {
			e = substitute(list(...))
			transformed = as.matrix(transform(as.data.frame(exprs(`_data`)),...))
			#Add any new parameter values
			param.names = colnames(transformed)
			newParams   = is.na(match(param.names,`_data`@parameters$name))
			params      = parameters(`_data`)
			if(any(newParams)) {
				params = rbind(params,
					data.frame(name=param.names[newParams]))
			}
			new("flowFrame",
				exprs=transfomed[,params$name],
				parameters=params,
				description=description(`_data`))
          })
