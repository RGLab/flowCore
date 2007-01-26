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
              if(ncol(transformed) > ncol(`_data`)) {
              	cnames = c(colnames(`_data`),colnames(transformed)[-c(1:ncol(`_data`))])
              }
              else {
              	cnames = colnames(`_data`)
              }
              ##param.names <- colnames(transformed)
              ##newParams <- is.na(match(param.names,`_data`@parameters$name))
              ##params      = parameters(`_data`)$name
              ##if(any(newParams)) {
              ##    params <- cbind(params,
              ##                    data.frame(name=param.names))
                  ##params = cbind(params,
                  ##    data.frame(name=param.names[newParams]))
              ##}
              #colnames(transformed) <- parameters(`_data`)$name
              colnames(transformed) = cnames
              new("flowFrame",
                  exprs=transformed, 
                  parameters=parameters(`_data`),#[,params$name],
                  description=description(`_data`))
          })
## ==========================================================================
## Transform function for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",signature=signature(`_data`="flowSet"),function(`_data`,...) {
    x = `_data`
    y = as(structure(lapply(seq(along=x),
      function(i) transform(`_data`=x[[i]],...)),names=sampleNames(phenoData(x))),"flowSet")
    phenoData(y) = phenoData(x)
    y
})
