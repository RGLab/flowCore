## The filterSummary class provides a container for the output of summary
## on a filterResult

## We need this to subset filterSummaries of multipleFilterResults
setMethod("[[",
          signature=signature("filterSummary","numeric"),
          definition=function(x,i,j,...)
      {
          list(name=as.vector(x@name[i]), true=as.vector(x@true[i]),
               false=as.vector(x@count[i]-x@true[i]),
               count=x@count[i], p=as.vector(x@p[i]),
               q=as.vector(1-x@p[i]))
      })

setMethod("[[",
          signature=signature("filterSummary","character"),
          definition=function(x,i,j,...)
      {
	i <- which(i==names(x))
	x[[i]]
    })


## Essentially the number of populations (in a multipleFilterResult)
setMethod("length",
          signature=signature("filterSummary"),
          definition=function(x) length(x@name))


## The names of the populations (in a multipleFilterResult)
setMethod("names",
          signature=signature("filterSummary"),
          definition=function(x) x@name)


## This allows for a list-like accessor to the slots (and more...)
## Valid values are 'n', 'count', 'true', 'false', 'name', 'p' and 'q'
setMethod("$",
          signature=signature("filterSummary","ANY"),
          definition=function(x,name)
      {
          switch(name, "n"=x@count, "true"=x@true, "in"=x@true,
                 "false"=x@count-x@true, "out"=x@count-x@true,
                 "p"=x@p, "q"=1-x@p, "count"=x@count, "name"=x@name)
      })


## Provide human-readable output 
setMethod("show",
          signature=signature("filterSummary"),
          definition=function(object) {
              if(length(object@name) == 1) {
                  cat(sprintf("%s: %d of %d events (%.2f%%)\n", object@name ,object@true,
                              object@count, object@p*100))
              } else {
                  for(i in seq(along=object@name)) {
                      cat(sprintf("%s: %d of %d events (%.2f%%)\n", object@name[i],
                                  object@true[i], object@count[i], object@p[i]*100))
                  }
              }
          })


## Return a more machine-readable output in form of a data.frame
setMethod("toTable",
          signature=signature("filterSummary"),
          definition=function(x, ...) {
              data.frame(count=x@count, true=x@true, false=x$false,
                         p=x$p, q=x$q)
          })
