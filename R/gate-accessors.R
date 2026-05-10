## =============================================================================
## Gate objects define boundaries for populations
## =============================================================================

#' @export
setMethod(
  "transform",
  signature = signature(`_data` = "rectangleGate"),
  definition = function(`_data`,
                        translist, 
                        ...) {
    
    gate <- `_data`
    if(!(missing(translist))){
      #check if it is a transformList
      res <- try(class(translist), silent = TRUE)
      if(res != "transformList") {
        err_msg <- attr(res, "condition")[["message"]]
        err_msg <- paste(
          err_msg, 
          "!Please make sure the unnamed argument is a valid 'transformList' object!"
        )
        stop(err_msg)
      } else {
        params <- unname(parameters(gate))
        trans_params <- colnames(translist)
        # transformation required
        for(param in params) {
          if(param %in% trans_params) {
            tfun <- translist@transforms[[param]]@f
            if(!is.infinite(gate@min[param])) {
              gate@min[param] <- tfun(gate@min[param])
            }
            if(!is.infinite(gate@max[param])) {
              gate@max[param] <- tfun(gate@max[param])
            }
          }
        }
        return(gate)
      }
      # apply named transformations of form `FSC-H`=asinhTrans(`FSC-H`)     
    } else {
      coords <- as.matrix(
        transform(
          as.data.frame(
            rbind(
              min = gate@min,
              max = gate@max
            )
          ),
          ...
        )
      )
      gate@min <- unlist(coords[1, , drop = TRUE])
      gate@max <- unlist(coords[2, , drop = TRUE])
      return(gate)
    }
  }
)

#' @export
setMethod(
  "transform",
  signature = signature(`_data` = "polygonGate"),
  definition = function(`_data`,
                        translist, 
                        ...) {
    
    gate <- `_data`
    if(!(missing(translist))){
      #check if it is a transformList
      res <- try(class(translist), silent = TRUE)
      if(res != "transformList") {
        err_msg <- attr(res, "condition")[["message"]]
        err_msg <- paste(
          err_msg, 
          "!Please make sure the unnamed argument is a valid 'transformList' object!"
        )
        stop(err_msg)
      } else {
        params <- unname(parameters(gate))
        trans_params <- colnames(translist)
        # transformation required
        for(param in params) {
          if(param %in% trans_params) {
            tfun <- translist@transforms[[param]]@f
            coords <- gate@boundaries[, param] 
            coords[!is.infinite(coords)] <- tfun(coords[!is.infinite(coords)])
            gate@boundaries[, param] <- coords
          }
        }
        return(gate)
      }
      # apply named transformations of form `FSC-H`=asinhTrans(`FSC-H`)     
    } else {
      gate@boundaries <- as.matrix(
        transform(
          as.data.frame(
            gate@boundaries
          ),
          ...
        )
      )
      return(gate)
    }
  }
)

#' @export
setMethod(
  "transform",
  signature = signature(`_data` = "ellipsoidGate"),
  definition = function(`_data`,
                        translist, 
                        ...) {
    
    # cannot preserve ellipse geometry after transform -> polygonGate
    gate <- as(`_data`, "polygonGate")
    gate <- transform(gate, translist = translist, ...)
    return(gate)
    
  }
)

#' @export
setMethod(
  "transform",
  signature = signature(`_data` = "quadGate"),
  definition = function(`_data`,
                        translist, 
                        ...) {
    
    gate <- `_data`
    if(!(missing(translist))){
      #check if it is a transformList
      res <- try(class(translist), silent = TRUE)
      if(res != "transformList") {
        err_msg <- attr(res, "condition")[["message"]]
        err_msg <- paste(
          err_msg, 
          "!Please make sure the unnamed argument is a valid 'transformList' object!"
        )
        stop(err_msg)
      } else {
        params <- unname(parameters(gate))
        trans_params <- colnames(translist)
        # transformation required
        for(param in params) {
          if(param %in% trans_params) {
            tfun <- translist@transforms[[param]]@f
            if(!is.infinite(gate@boundary[param])) {
              gate@boundary[param] <- tfun(gate@boundary[param])
            }
          }
        }
        return(gate)
      }
      # apply named transformations of form `FSC-H`=asinhTrans(`FSC-H`)     
    } else {
      coords <- as.matrix(
        transform(
          as.data.frame(
            rbind(
              gate@boundary
            )
          ),
          ...
        )
      )
      gate@boundary <- unlist(coords[1, , drop = TRUE])
      return(gate)
    }
  }
)

#' @export
setMethod(
  "transform",
  signature = signature(`_data` = "filters"),
  definition = function(`_data`,
                        translist, 
                        ...) {
    
    # transform each gate in filters list
    gate <- filters(
      lapply(
        unlist(`_data`),
        function(z) {
          transform(
            z,
            translist = translist,
            ...
          )
        }
      )
    )
    return(gate)
  }
)
