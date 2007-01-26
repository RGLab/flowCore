## NLM Jan 17
## ==========================================================================
## Transformation function for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...) {
              e <- substitute(list(...))
              x = `_data`
              transformed <- as.matrix(transform(as.data.frame(exprs(x)),...))
              ##Add any new parameter values
              if(ncol(transformed) > ncol(x)) {
              	cnames = c(colnames(x),colnames(transformed)[-c(1:ncol(x))])
              }
              else {
              	cnames = colnames(x)
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
                  parameters=parameters(x),#[,params$name],
                  description=description(x))
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
