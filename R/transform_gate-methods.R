## ===========================================================================
## Generics for modifying gates
## ---------------------------------------------------------------------------
#' @export
transform_gate <- function(obj, ...) UseMethod("transform_gate")
#' Simplified geometric transformation of gates
#' 
#' Perform geometric transformations of Gate-type \code{\linkS4class{filter}} objects
#'
#' This method allows changes to the four filter types defined by simple geometric gates (\code{\linkS4class{quadGate}},
#' \code{\linkS4class{rectangleGate}}, \code{\linkS4class{ellipsoidGate}}, and \code{\linkS4class{polygonGate}}) using
#' equally simple geometric transformations (shifting/translation, scaling/dilation, and rotation). The method also
#' allows for directly re-setting the slots of each Gate-type object. Note that these methods are for manually altering
#' the geometric definition of a gate. To easily transform the definition of a gate with an accompanyging scale 
#' transformation applied to its underlying data, see \code{\link[ggcyto]{rescale_gate}}.
#' 
#' First, \code{transform_gate} will apply any direct alterations to the slots of the supplied Gate-type filter object.
#' For example, if "\code{mean = c(1,3)}" is present in the argument list when \code{transform_gate} is called on a
#' \code{ellipsoidGate} object, the first change applied will be to shift the \code{mean} slot to \code{(1,3)}. The method
#' will carry over the dimension names from the gate, so there is no need to provide column or row names with arguments
#' such as \code{mean} or \code{cov} for \code{ellipsoidGate} or \code{boundaries} for \code{polygonGate}.
#' 
#' \code{transform_gate} then passes the geometric arguments (\code{dx}, \code{dy}, \code{deg}, \code{rot_center}, \code{scale}, 
#' and \code{center}) to the methods which perform each respective type of transformation:  
#' \code{\link{shift_gate}}, \code{\link{scale_gate}}, or \code{\link{rotate_gate}}. The order of operations is to first
#' scale, then rotate, then shift. The default behavior of each operation follows that of its corresponding method but for
#' the most part these are what the user would expect. A few quick notes:
#' \itemize{
#' \item \code{rotate_gate} is not defined for \code{rectangleGate} or \code{quadGate} objects, due to their definition as
#' having 1-dimensional boundaries.
#' \item The default center for both rotation and scaling of a \code{polygonGate} is the centroid of the polygon. This
#' results in the sort of scaling most users expect, with a uniform scale factor not distorting the shape of the original polygon.
#' }
#' 
#' 
#' @name transform_gate
#' 
#' @param obj A Gate-type \code{\link{filter}} object (\code{\linkS4class{quadGate}},
#' \code{\linkS4class{rectangleGate}}, \code{\linkS4class{ellipsoidGate}}, or \code{\linkS4class{polygonGate}})
#' 
#' @param scale Either a numeric scalar (for uniform scaling in all dimensions) or numeric vector specifying the factor by 
#' which each dimension of the gate should be expanded (absolute value > 1) or contracted (absolute value < 1). Negative values 
#' will result in a reflection in that dimension. 
#' 
#' For \code{rectangleGate} and \code{quadGate} objects, this amounts to simply
#' scaling the values of the 1-dimensional boundaries. For \code{polygonGate} objects, the values of \code{scale} will be used
#' to determine scale factors in the direction of each of the 2 dimensions of the gate (\code{scale_gate} is not yet defined
#' for higher-dimensional \code{polytopeGate} objects). \strong{Important: } For \code{ellipsoidGate} objects, \code{scale}
#' determines scale factors for the major and minor axes of the ellipse, in that order.
#' 
#' @param deg An angle in degrees by which the gate should be rotated in the counter-clockwise direction.
#' @param rot_center A separate 2-dimensional center of rotation for the gate, if desired. By default, this will
#' be the center for \code{ellipsoidGate} objects or the centroid for \code{polygonGate} objects. The \code{rot_center} argument 
#' is currently only supported for \code{polygonGate} objects. It is also usually simpler to perform a rotation and a translation 
#' individually than to manually specify the composition as a rotation around a shifted center.
#' 
#' @param dx Either a numeric scalar or numeric vector. If it is scalar, this is just the desired shift of the gate in 
#' its first dimension. If it is a vector, it specifies both \code{dx} and \code{dy} as \code{(dx,dy)}.
#' This provides an alternate syntax for shifting gates, as well as allowing shifts of \code{ellipsoidGate} objects
#' in more than 2 dimensions.
#' @param dy A numeric scalar specifying the desired shift of the gate in its second dimension.
#' @param center A numeric vector specifying where the center or centroid should be moved (rather than specifiying \code{dx} 
#' and/or \code{dy})
#' 
#' @param \dots Assignments made to the slots of the particular Gate-type filter object in the form "<slot_name> = <value>"
#' 
#' @return A Gate-type \code{filter} object of the same type as \code{gate}, with the geometric transformations applied
#'
#' @examples
#' \dontrun{
#' # Scale the original gate non-uniformly, rotate it 15 degrees, and shift it
#' transformed_gate <- transform_gate(original_gate, scale = c(2,3), deg = 15, dx = 500, dy = -700)
#' 
#' # Scale the original gate (in this case an ellipsoidGate) after moving its center to (1500, 2000)
#' transformed_gate <- transform_gate(original_gate, scale = c(2,3), mean = c(1500, 2000))
#' }
#'
#' @export
transform_gate.default <- function(obj, scale = NULL, deg = NULL, rot_center = NULL, dx = NULL, dy = NULL, center = NULL, ...){
	gate <- obj
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
scale_gate <- function(obj, ...) UseMethod("scale_gate")

#' Simplified geometric scaling of gates
#' 
#' Scale a Gate-type filter object in one or more dimensions
#' 
#' This method allows uniform or non-uniform geometric scaling of filter types defined by simple geometric gates 
#' (\code{\linkS4class{quadGate}}, \code{\linkS4class{rectangleGate}}, \code{\linkS4class{ellipsoidGate}}, and 
#' \code{\linkS4class{polygonGate}}) Note that these methods are for manually altering
#' the geometric definition of a gate. To easily transform the definition of a gate with an accompanyging scale 
#' transformation applied to its underlying data, see \code{\link[ggcyto]{rescale_gate}}.
#' 
#' The \code{scale} argument passed to \code{scale_gate} should be either a scalar or a vector of the same length
#' as the number of dimensions of the gate. If it is scalar, all dimensions will be multiplicatively scaled uniformly
#' by the scalar factor provided. If it is a vector, each dimension will be scaled by its corresponding entry in the vector.
#' 
#' The scaling behavior of \code{scale_gate} depends on the type of gate passed to it. For \code{rectangleGate} 
#' and \code{quadGate} objects, this amounts to simply scaling the values of the 1-dimensional boundaries. 
#' For \code{polygonGate} objects, the values of \code{scale} will be used to determine scale factors 
#' in the direction of each of the 2 dimensions of the gate (\code{scale_gate} is not yet defined
#' for higher-dimensional \code{polytopeGate} objects). \strong{Important: } For \code{ellipsoidGate} objects, \code{scale}
#' determines scale factors for the major and minor axes of the ellipse, \emph{in that order}. Scaling by a negative factor 
#' will result in a reflection in the corresponding dimension.
#' 
#' @name scale_gate
#' 
#' @param obj A Gate-type \code{\link{filter}} object (\code{\linkS4class{quadGate}},
#' \code{\linkS4class{rectangleGate}}, \code{\linkS4class{ellipsoidGate}}, or \code{\linkS4class{polygonGate}})
#' 
#' @param scale Either a numeric scalar (for uniform scaling in all dimensions) or numeric vector specifying the factor by 
#' which each dimension of the gate should be expanded (absolute value > 1) or contracted (absolute value < 1). Negative values 
#' will result in a reflection in that dimension. 
#' 
#' @param \dots Additional arguments not used
#' 
#' @return A Gate-type \code{filter} object of the same type as \code{gate}, with the scaling applied
#' 
#' @examples
#' \dontrun{
#' # Scales both dimensions by a factor of 5
#' scaled_gate <- scale_gate(original_gate, 5)
#' 
#' # Shrinks the gate in the first dimension by factor of 1/2
#' # and expands it in the other dimension by factor of 3
#' scaled_gate <- scale_gate(original_gate, c(0.5,3))
#' }
#' 
#' @export
scale_gate.default <- function(obj, scale = NULL, ...){
  stop("scale_gate() does not support this gate type")
}

#' @export
rotate_gate <- function(obj, ...) UseMethod("rotate_gate")

#' Simplified geometric rotation of gates
#' 
#' Rotate a Gate-type filter object through a specified angle
#' 
#' This method allows for 2-dimensional geometric rotation of filter types defined by simple geometric gates 
#' (\code{\linkS4class{ellipsoidGate}}, and \code{\linkS4class{polygonGate}}). The method is not defined 
#' for \code{rectangleGate} or \code{quadGate} objects, due to their definition as having 1-dimensional boundaries.
#' Further, keep in mind that the 2-dimensional rotation takes place in the plane where the dimensions
#' of the two variables are evenly linearly scaled. Displaying a rotated ellipse in a plot where the axes are not scaled
#' evenly may make it appear that the ellipse has been distorted even though this is not the case. 
#' 
#' The angle provided in the \code{deg} argument should be in degrees rather than radians. By default, the rotation
#' will be performed around the center of an \code{ellipsoidGate} or the centroid of the area encompassed by
#' a \code{polygonGate}. The \code{rot_center} argument allows for specification of a different center of rotation
#' for \code{polygonGate} objects (it is not yet implemented for \code{ellipsoidGate} objects) but
#' it is usually simpler to perform a rotation and a translation individually than to manually specify 
#' the composition as a rotation around a shifted center.
#' 
#' 
#' @name rotate_gate
#' 
#' @param obj An \code{\linkS4class{ellipsoidGate}} or \code{\linkS4class{polygonGate}}
#' 
#' @param deg An angle in degrees by which the gate should be rotated in the counter-clockwise direction
#' @param rot_center A separate 2-dimensional center of rotation for the gate, if desired. By default, this will
#' be the center for \code{ellipsoidGate} objects or the centroid for \code{polygonGate} objects. The \code{rot_center} argument 
#' is currently only supported for \code{polygonGate} objects.
#' 
#' @param \dots Additional arguments not used
#' 
#' @return A Gate-type \code{filter} object of the same type as \code{gate}, with the rotation applied
#' 
#' @examples
#' \dontrun{
#' #' # Rotates the original gate 15 degrees counter-clockwise
#' rotated_gate <- rotate_gate(original_gate, deg = 15)
#' # Rotates the original gate 270 degrees counter-clockwise
#' rotated_gate <- rotate_gate(original_gate, 270)
#' }
#' 
#' @export
rotate_gate.default <- function(obj, deg = NULL, rot_center = NULL, ...){
  stop("rotate_gate() does not support this gate type")
}

#' @export
shift_gate <- function(obj, ...) UseMethod("shift_gate")

#' Simplified geometric translation of gates
#' 
#' Shift a Gate-type filter object in one or more dimensions
#' 
#' This method allows for geometric translation of filter types defined by simple geometric gates 
#' (\code{rectangleGate}, \code{quadGate}, \code{\linkS4class{ellipsoidGate}}, or \code{\linkS4class{polygonGate}}).
#' The method provides two approaches to specify a translation. For \code{rectangleGate} objects, this will
#' shift the \code{min} and \code{max} bounds by the same amount in each specified dimension. For \code{quadGate}
#' objects, this will simply shift the divinding boundary in each dimension. For \code{ellipsoidGate} objects, this
#' will shift the center (and therefore all points of the ellipse). For \code{polgonGate} objects, this will simply
#' shift all of the points defining the polygon.
#' 
#' The method allows two different approaches to shifting a gate. Through the \code{dx} and/or \code{dy} arguments,
#' a direct shift in each dimension can be provided. Alternatively, through the \code{center} argument, the gate
#' can be directly moved to a new location in relation to the old center of the gate. For \code{quadGate} objects, 
#' this center is the intersection of the two dividing boundaries (so the value of the \code{boundary} slot). For
#' \code{rectangleGate} objects, this is the center of the rectangle defined by the intersections of the centers
#' of each interval. For \code{ellipsoidGate} objects, it is the center of the ellipsoid, given by the \code{mean}
#' slot. For \code{polygonGate} objects, the centroid of the old polygon will be calculated and shifted to the new
#' location provided by \code{center} and all other points on the polygon will be shifted by relation to the centroid.
#' 
#' @name shift_gate
#' 
#' @param obj A Gate-type \code{\link{filter}} object (\code{\linkS4class{quadGate}},
#' \code{\linkS4class{rectangleGate}}, \code{\linkS4class{ellipsoidGate}}, or \code{\linkS4class{polygonGate}})
#' 
#' @param dx Either a numeric scalar or numeric vector. If it is scalar, this is just the desired shift of the gate in 
#' its first dimension. If it is a vector, it specifies both \code{dx} and \code{dy} as \code{(dx,dy)}.
#' This provides an alternate syntax for shifting gates, as well as allowing shifts of \code{ellipsoidGate} objects
#' in more than 2 dimensions.
#' @param dy A numeric scalar specifying the desired shift of the gate in its second dimension.
#' @param center A numeric vector specifying where the center or centroid should be moved (rather than specifiying \code{dx} 
#' and/or \code{dy})
#' @param \dots Additional arguments not used
#' 
#' @return A Gate-type \code{filter} object of the same type as \code{gate}, with the translation applied
#' 
#' @examples
#' \dontrun{
#' # Moves the entire gate +500 in its first dimension and 0 in its second dimension
#' shifted_gate <- shift_gate(original_gate, dx = 500)
#' 
#' #Moves the entire gate +250 in its first dimension and +700 in its second dimension
#' shifted_gate <- shift_gate(original_gate, dx = 500, dy = 700)
#' 
#' # Same as previous
#' shifted_gate <- shift_gate(original_gate, c(500,700))
#' 
#' # Move the gate based on shifting its center to (700, 1000)
#' shifted_gate <- shift_gate(original_gate, center = c(700, 1000))
#' }
#' 
#' @export
shift_gate.default <- function(obj, dx=NULL, dy=NULL, center=NULL, ...){
  stop("shift_gate() does not support this gate type")
}

# TRANSLATION METHODS

#' @noRd
#' @export
shift_gate.ellipsoidGate <- function(obj, dx=NULL, dy=NULL, center=NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
shift_gate.rectangleGate <- function(obj, dx=NULL, dy=NULL, center=NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
shift_gate.polygonGate <- function(obj, dx=NULL, dy=NULL, center=NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
shift_gate.quadGate <- function(obj, dx=NULL, dy=NULL, center=NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
scale_gate.rectangleGate <- function(obj, scale = NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
scale_gate.polygonGate <- function(obj, scale = NULL, ...){
	gate <- obj
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

# Note scale should go c(major axis, minor axis) because eigen output is descending
#' @noRd
#' @export
scale_gate.ellipsoidGate <- function(obj, scale = NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
scale_gate.quadGate <- function(obj, scale = NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
rotate_gate.quadGate <- function(obj, deg = NULL, rot_center = NULL, ...){
  stop("rotate_gate is not defined for quadGate objects")
}

#' @noRd
#' @export
rotate_gate.rectangleGate <- function(obj, deg = NULL, rot_center = NULL, ...){
  stop("rotate_gate is not defined for rectangleGate objects")
}

# For ellipsoidGate, only allowing rotation about ellipse center for now
# Also, degree angular measure (not radian)
#' @noRd
#' @export
rotate_gate.ellipsoidGate <- function(obj, deg = NULL, rot_center = NULL, ...){
	gate <- obj
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

#' @noRd
#' @export
rotate_gate.polygonGate <- function(obj, deg = NULL, rot_center = NULL, ...){
	gate <- obj
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