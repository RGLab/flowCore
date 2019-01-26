#' @export
transform_gate.default <- function(gate, scale = NULL, deg = NULL, rot_center = NULL, dx = NULL, dy = NULL, center = NULL, ...){
  if(!(class(gate) %in% c("rectangleGate", "quadGate", "ellipsoidGate", "polygonGate"))){
    stop("The object passed to transform_gate is not of a supported gate type")
  }
  params <- list(...)
  if(length(params) > 0){
    keys <- names(params)
    for(i in 1:length(params)){
      # Transfer names appropriately so user doesn't have to
      if(!is.null(names(slot(gate, keys[i]))) && is.null(names(params[[i]]))) names(params[[i]]) <- names(slot(gate, keys[i]))
      if(!is.null(colnames(slot(gate, keys[i]))) && is.null(colnames(params[[i]]))) colnames(params[[i]]) <- colnames(slot(gate, keys[i]))
      if(!is.null(rownames(slot(gate, keys[i]))) && is.null(rownames(params[[i]]))) rownames(params[[i]]) <- rownames(slot(gate, keys[i]))
      slot(gate, keys[i]) <- params[[i]]
    }
  }
  if(!is.null(scale)){
    gate <- scale_gate(gate, scale)
  }
  if(length(c(deg, rot_center)) > 0){
    gate <- rotate_gate(gate, deg, rot_center)
  }
  if(length(c(dx,dy,center) > 0)){
    gate <- shift_gate(gate, dx=dx, dy=dy, center=center)
  }
  return(gate) 
}

#' @export
scale_gate.default <- function(gate, ...){
  stop("scale_gate() does not support this gate type")
}

#' @export
rotate_gate.default <- function(gate, ...){
  stop("rotate_gate() does not support this gate type")
}

#' @export
shift_gate.default <- function(gate, ...){
  stop("shift_gate() does not support this gate type")
}

# TRANSLATION METHODS

#' @export
shift_gate.ellipsoidGate <- function(gate, dx=NULL, dy=NULL, center=NULL){
  if(!is.null(center)){
    if(length(center) == length(gate@mean)){
      names(center) <- names(gate@mean)
      gate@mean <- center
    }else{
      stop("Length of center must match number of parameters in gate")
    }
  }else{
    if(!is.null(dx)){
      if(length(dx) == 1){
        gate@mean[1] <- gate@mean[1] + dx 
      }else if (length(dx) == length(gate@mean)){
        if(!is.null(dy)){
          stop("If dx is non-scalar, no value can be specified for dy")
        }
        gate@mean <- gate@mean + dx
      }else{
        stop("Length of list of changes for mean must be same as number of parameters in gate")
      }
    }
    if(!is.null(dy)){
      if(length(gate@mean) ==  1){
        stop("dy cannot be specified for 1-dimensional gate")
      }
      if(length(dy) == 1){
        gate@mean[2] <- gate@mean[2] + dy
      }else{
        stop("dy must be numeric scalar")
      }
    }
  }
  return(gate)
}

#' @export
shift_gate.rectangleGate <- function(gate, dx=NULL, dy=NULL, center=NULL){
  if(!is.null(center)){
    if(any(is.infinite(gate@min), is.infinite(gate@max))){
      stop("Cannot shift center of rectangleGate with an infinite bound (has no finite center)")
    }else{
      if(length(center) == length(gate@min)){
        old_center <- (gate@max + gate@min) / 2
        diff <- center - old_center
        gate@min <- gate@min + diff
        gate@max <- gate@max + diff
      }else{
        stop("Length of center must match number of parameters in gate")
      }
    }
  }else{
    if(!is.null(dx)){
      if(length(dx) == 1){
        gate@min[1] <- gate@min[1] + dx
        gate@max[1] <- gate@max[1] + dx
      }else if (length(dx) == length(gate@min)){
        if(!is.null(dy)){
          stop("If dx is non-scalar, no value can be specified for dy")
        }
        gate@min <- gate@min + dx
        gate@max <- gate@max + dx
      }else{
        stop("Length of list of changes for mean must be same as number of parameters in gate")
      }
    }
    if(!is.null(dy)){
      if(length(gate@min) == 1){
        stop("dy cannot be specified for 1-dimensional gate")
      }
      if(length(dy) == 1){
        gate@min[2] <- gate@min[2] + dy
        gate@max[2] <- gate@max[2] + dy
      }else{
        stop("dy must be numeric scalar")
      }
    }
  }
  return(gate)
}

#' @export
shift_gate.polygonGate <- function(gate, dx=NULL, dy=NULL, center=NULL){
  if(!is.null(center)){
    if(any(is.infinite(gate@boundaries))){
      stop("Cannot shift center of polygonGate with an infinite bound (has no finite centroid)")
    }else{
      if(length(center) == 2){
        old_center <- poly_centroid(gate@boundaries)
        # Translation vector
        diff <- center - old_center
        gate@boundaries <- sweep(gate@boundaries, 2, -diff)
      }else{
        stop("Invalid number of columns in center argument. Shifting polygon by centroid only valid in 2 dimensions")
      }
      
    }
  }else{
    if(!is.null(dx)){
      if(length(dx) == 1){
        gate@boundaries[, 1] <- gate@boundaries[, 1] + dx
      }else if (length(dx) == ncol(gate@boundaries)){
        if(!is.null(dy)){
          stop("If dx is non-scalar, no value can be specified for dy")
        }
        gate@boundaries <- sweep(gate@boundaries, 2, -dx)
      }else{
        stop("Length of list of changes for mean must be same as number of parameters in gate")
      }
    }
    if(!is.null(dy)){
      if(nrow(gate@boundaries) ==  1){
        stop("dy cannot be specified for 1-dimensional gate")
      }
      if(length(dy) == 1){
        gate@boundaries[, 2] <- gate@boundaries[, 2] + dy
      }else{
        stop("dy must be numeric scalar")
      }
    }
  }
  return(gate)
}

#' @export
shift_gate.quadGate <- function(gate, dx=NULL, dy=NULL, center=NULL){
  if(!is.null(center)){
    if(any(is.infinite(gate@boundary))){
      stop("Cannot shift center of quadGate with an infinite bound (has no finite center)")
    }else{
      if(length(center) == 2){
        names(center) <- names(gate@boundary)
        gate@boundary <- center
      }else{
        stop("Length of center must be precisely 2")
      }
    }
  }else{
    if(!is.null(dx)){
      if(length(dx) == 1){
        gate@boundary[1] <- gate@boundary[1] + dx
      }else if (length(dx) == 2){
        if(!is.null(dy)){
          stop("If dx is non-scalar, no value can be specified for dy")
        }
        gate@boundary[1] <- gate@boundary[1] + dx[1]
        gate@boundary[2] <- gate@boundary[2] + dx[2]
      }else{
        stop("Length of list of changes for mean must be same as number of parameters in gate")
      }
    }
    if(!is.null(dy)){
      if(length(dy) == 1){
        gate@boundary[2] <- gate@boundary[2] + dy
      }else{
        stop("dy must be numeric scalar")
      }
    }
  }
  return(gate)
}

# SCALING METHODS

#' @export
scale_gate.rectangleGate <- function(gate, scale = NULL){
  if(any(is.infinite(gate@min), is.infinite(gate@max))){
    stop("Cannot resize rectangleGate with an infinite bound")
  }
  # Factor can be scalar or vector
  if(!is.null(scale)){
  span <- gate@max - gate@min
  new_span <- span
  center <- (gate@max + gate@min) / 2
    if(length(scale) == 1 || length(scale) == length(gate@min)){
      new_span <- scale*span    
    }else{
      stop("Length of scale must be 1 or same as number of parameters in gate")
    }
  
  half_span <- new_span / 2
  gate@min <- center - half_span
  gate@max <- center + half_span
  }
  return(gate)
}

#' @export
scale_gate.polygonGate <- function(gate, scale = NULL){
  if(any(is.infinite(gate@boundaries))){
    stop("Cannot resize polygonGate with an infinite bound")
  }
  if(!is.null(scale)){
    centroid <- poly_centroid(gate@boundaries)
    diff <- sweep(gate@boundaries, 2, centroid)
    
    # Perform the dilation/contraction
    if(length(scale) %in% c(1,2)){
      # Scale the radial vectors and shift back to centroid
      gate@boundaries <- sweep(sweep(diff, 2, scale, FUN="*"), 2, -centroid)
    }else{
      stop("Too many entries in scale argument. Scaling polygonGate only valid in 2 dimensions")
    }
  }
  return(gate)
}

# Note scale should go c(major axis, minor axis) because of eigen output is descending
#' @export
scale_gate.ellipsoidGate <- function(gate, scale = NULL){
  eigs <- eigen(gate@cov, symmetric = TRUE)
  # Factor can be scalar or vector
  if(!is.null(scale)){
    if(length(scale) == 1 || length(scale) == length(gate@mean)){
      eigs$values <- scale*eigs$values    
    }else{
      stop("Length of scale must be 1 or same as number of parameters in gate")
    }
    # Update cov
    Q <- eigs$vectors
    new_cov <- Q%*%(diag(eigs$values))%*%solve(Q)
    colnames(new_cov) <- colnames(gate@cov)
    rownames(new_cov) <- rownames(gate@cov)
    gate@cov <- new_cov
  }
  return(gate)
}

#' @export
scale_gate.quadGate <- function(gate, scale = NULL){
  # Factor can be scalar or vector of length 2
  if(!is.null(scale)){
    if(length(scale) == 1 || length(scale) == 2){
      gate@boundary <- gate@boundary*scale    
    }else{
      stop("Length of scale must be 1 (uniform scaling) or 2 (an entry for each boundary)")
    }
  }
  return(gate)
}

# ROTATION METHODS

#' @export
rotate_gate.quadGate <- function(gate,  ...){
  stop("rotate_gate is not defined for quadGate objects")
}

#' @export
rotate_gate.rectangleGate <- function(gate, ...){
  stop("rotate_gate is not defined for rectangleGate objects")
}

# For ellipsoidGate, only allowing rotation about ellipse center for now
# Also, degree angular measure (not radian)
#' @export
rotate_gate.ellipsoidGate <- function(gate, deg = NULL, rot_center = NULL){
  if(!is.null(rot_center)){
    stop("rotate_gate only allows rotation about the ellipse center for ellipsoidGate (center argument not allowed)") 
  }
  if(!is.null(deg)){
    rad <- deg*(pi/180)
    rot <- rbind(c(cos(rad), -sin(rad)), c(sin(rad), cos(rad)))
    new_cov <- rot%*%(gate@cov)%*%t(rot)
    colnames(new_cov) <- colnames(gate@cov)
    rownames(new_cov) <- rownames(gate@cov)
    gate@cov <- new_cov
  }
  return(gate)
}

#' @export
rotate_gate.polygonGate <- function(gate, deg = NULL, rot_center = NULL){
  if(!is.null(rot_center)){
    if(length(rot_center) != 2){
      stop("If rot_center is specified, it must be of length 2")
    }
  }else{
    rot_center <- poly_centroid(gate@boundaries)
  }
  # Rotate
  if(!is.null(deg)){
    diff <- sweep(gate@boundaries, 2, rot_center)
    rad <- deg*(pi/180)
    rot <- rbind(c(cos(rad), -sin(rad)), c(sin(rad), cos(rad)))
    diff <- t(rot%*%t(diff))
    colnames(diff) <- colnames(gate@boundaries)
    # Shift back
    gate@boundaries <- sweep(diff, 2, -rot_center)
  }
  return(gate)
}