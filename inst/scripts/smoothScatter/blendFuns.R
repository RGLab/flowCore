## Here are a number of basic image processing and blending
## algorithms. Most of them are based on the documentation to
## the GIMP software but they should be rather general (although
## image processing is a hell of a complex field!!!)

require(RColorBrewer)

## A wrapper around blending functions. Converts color matrices to rgb,
## applies the blending and resets the dimensions of the resulting,
## color vector. This allows to feed bitmaps to the blending functions.
blendWrapper <- function(background, foreground, blendFun, ...){
  d1 <- dim(background)
  d2 <- dim(foreground)
  stopifnot(all(d1==d2))
  bgrgb <- col2rgb(background)
  fgrgb <- col2rgb(foreground)
  blendcol <- blendFun(bg=bgrgb, fg=fgrgb, ...)
  if(!is.null(d1) && !is.null(d2))
    res <- matrix(rgb(blendcol[1,], blendcol[2,], blendcol[3,], max=255),
                  nrow=d1[1], ncol=d1[2])
  else
    res <- rgb(blendcol[1,], blendcol[2,], blendcol[3,], max=255) 
  return(res)
}

##------------------------------------------------------------------------
## Some usefull image processing functions:

## The brightness of a color: a percentage
brightness <- function(col){
  signif((colSums(col)/3)/2.55, 3)
}

## The lightness of a pixel: a percentage
lightness <- function(col){
  signif(((apply(col, 2, max)+apply(col, 2, min))/2)/2.55, 3)
}

## Best reflects the human perception of
## the brighntess of a color: a percentage
luminance <- function(col){
  signif(colSums(col*c(0.3, 0.59, 0.11))/2.55, 3)
}

## The actual color component in hsv colorspace:
## a value from 0 to 360
hue <- function(col){
  col <- col/255
  cmax <- apply(col, 2, max)
  cmin <- apply(col, 2, min)
  delta <- cmax - cmin
  ##h <- rep(0, ncol(col))
  undef <- cmax == 0
  ym <- col[1,]==cmax
  cm <- col[2,]==cmax
  h <- 4+(col[1,]-col[2,])/delta
  h[ym] <- (col[2,ym]-col[3,ym])/delta[ym]
  h[cm] <- 2+(col[3,cm]-col[1,cm])/delta[cm]
  h <- h*60
  sel <- h < 0
  h[sel] <- h[sel]+360
  return(round(h))
}

## The relative colorfulness with respect to brightness:
## a percentage
saturation <- function(col){
  cmax <- apply(col, 2, max)
  cmin <- apply(col, 2, min)
  delta <- cmax - cmin
  s <- rep(0, ncol(col))
  sel <- cmax != 0
  s[sel] = delta/cmax
  return(signif(s*100,3))
}

## The brightness of a color in the hsv space: a percentage
value <- function(col){
  v <- apply(col, 2, max)
  return(signif(v/2.55,3))
}

## invert rgb colors
invert <- function(col){
  return(255-col)
}

## spread rgb colors over the available color range
spread <- function(col){
  d <- dim(col)
  colv <- as.vector(col)
  colv <- round((colv - min(colv))/diff(range(colv)) * 255)
  dim(colv) <- d
  return(colv)
}


############################################################
## The following functions are the actual workhorses for
## blending of colors. They expect matrices of rgb values
## (from 0 to 255) with 3 rows and n>0 columns. Use
## in combination with the blendWrapper function.
###########################################################

## Adds the rgb values of bg and fg color. Blown out values
## are clipped to 255. This lightens an image.
bl_addition <- function(bg, fg){
  d <- dim(bg)
  addcol <- pmin(255, bg+fg)
  dim(addcol) <- d
  return(addcol)
}

## Substracts the fg rgb values from the bg. Negative values
## are clipped to 0. This darkens an image. Note: It is not
## symmetric (bg-fg != fg-bg).
bl_substraction <- function(bg, fg){
  d <- dim(bg)
  addcol <- pmax(0, bg-fg)
  dim(addcol) <- d
  return(addcol)
}

## Absolute difference of rgb values. Colors are darkened but
## don't blow out to black.
bl_difference <- function(bg, fg){
  addcol <- abs(bg-fg)
  return(addcol)
}

## Multiplication of rgb values and normalization back into
## the rgb color cube. This darkens an image
bl_burn <- function(bg, fg){
  addcol <- (bg*fg)/255
  return(addcol)
}

## Division of bg values by foreground values and normalization back into
## the rgb color cube. Division by 0 is avoided by adding 1 and blown
## out colors are clipped. This lightens an image.
bl_dodge <- function(bg, fg){
  d <- dim(bg)
  addcol <- pmin(255, bg/((fg+1)/256))
  dim(addcol) <- d
  return(addcol)
}

## The inverse of 'burn'. This lightens an image.
bl_screen <- function(bg, fg){
  addcol <- 255-((1/255)*((255-bg)*(255-fg)))
  return(addcol)
}

## Combination of 'screen' and 'burn'. This increases the
## contrast of the image.
##FIXME: Produces rubbish, need to check algorithm
bl_overlay <- function(bg, fg){
  addcol <- 1/255*(bg*bl_screen(bg, fg)+(1-bg)*bl_burn(bg, fg))
  return(addcol)
}

## The paralell minima of the rgb values. This darkens an image.
bl_darkenOnly <- function(bg, fg){
  d <- dim(bg)
  addcol <- pmin(bg, fg)
  dim(addcol) <- d
  return(addcol)
}

## The paralell maxima of the rgb values. This lightens an image.
bl_lightenOnly <- function(bg, fg){
  d <- dim(bg)
  addcol <- pmax(bg, fg)
  dim(addcol) <- d
  return(addcol)
}

## Keep saturation and value of bg and hue of fg. 
bl_hue <- function(bg, fg){
  bghsv <- rgb2hsv(bg)
  fghsv <- rgb2hsv(fg)
  addcol <- col2rgb(hsv(fghsv[1,], bghsv[2,], bghsv[3,]))
  return(addcol)
}



## Keep hue and value of bg and saturation of fg. 
bl_saturation <- function(bg, fg){
  bghsv <- rgb2hsv(bg)
  fghsv <- rgb2hsv(fg)
  addcol <- col2rgb(hsv(bghsv[1,], fghsv[2,], bghsv[3,]))
  return(addcol)
}

## Keep hue and saturation of bg and value of fg. 
bl_value <- function(bg, fg){
  bghsv <- rgb2hsv(bg)
  fghsv <- rgb2hsv(fg)
  addcol <- col2rgb(hsv(bghsv[1,], bghsv[2,], fghsv[3,]))
  return(addcol)
}

## Combine fg hue and saturation with bg lightness.
##FIXME: This doesn't look right, need to check algorithm
bl_color<- function(bg, fg){
  bghsv <- rgb2hsv(bg)
  fghsv <- rgb2hsv(fg)
  addcol <- col2rgb(hsv(fghsv[1,], fghsv[2,], lightness(bg)/255))
  return(addcol)
}

## Blending of fg with bg according to fg alpha value. Alpha=0
## means transparent fg, alpha=1 means opaque fg.
bl_alpha <- function(bg,fg, alpha=0.5){
  addcol <- bg*(1-alpha)+fg*alpha
  return(addcol)
}
