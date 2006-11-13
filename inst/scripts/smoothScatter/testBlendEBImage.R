## Use EBImage to test the blending algorithms on 'real' images:
## the convenience fucntion Ebiblend take as arguments two EBImage
## images of exactly the same size (pixel dimensions) and applies
## the blending function blendFun to them. The output of the funcion
## Is the blended image which can be inspected using the display
## function.
library(EBImage)
source("blendFuns.R")

## A wrapper around the actual color blending functions
rgbBlendWrapper <- function(background, foreground, blendFun, ...){
  d1 <- dim(background)
  d2 <- dim(foreground)
  stopifnot(all(d1==d2))
  blendFun(background, foreground, ...)
}

## A function that uses EBImage to handle images for blending
EbiBlend <- function(background, foreground, blendFun, ...){ 
  fgChannels <- t(sapply(channels(foreground), as.vector))*255
  bgChannels <- t(sapply(channels(background), as.vector))*255
  blendedRGB <- rgbBlendWrapper(bgChannels, fgChannels, blendFun, ...)/255
  blendedChannels <- list(r=blendedRGB[1,], g=blendedRGB[2,],
                          b=blendedRGB[3,])
  blendedImg <- toRGB(toRed(blendedRGB[1,])+toGreen(blendedRGB[2,])+
                toBlue(blendedRGB[3,]))
  dim(blendedImg) <- dim(background)
  return(Image(blendedImg, rgb=TRUE))
}
##------------------------------------------------------------------------
## read in a couple of example files to play around with
fg <- read.image("foreground.jpg", rgb=TRUE)
bggray <- read.image("backgray.jpg", rgb=TRUE)
bghue <- read.image("backhue.jpg", rgb=TRUE)
bgsat <- read.image("backsat.jpg", rgb=TRUE)
bgval <- read.image("backval.jpg", rgb=TRUE)
bgcol <- read.image("backcol.jpg", rgb=TRUE)
##-------------------------------------------------------------------------
## apply blending algorithms and show results (note: you need to
## manually close an old image windows in order to display a new one)
display(fg)

display(EbiBlend(bggray, fg, bl_addition))

display(EbiBlend(bggray, fg, bl_substraction))

display(EbiBlend(bggray, fg, bl_difference))

display(EbiBlend(bggray, fg, bl_burn))

display(EbiBlend(bggray, fg, bl_dodge))

display(EbiBlend(bggray, fg, bl_screen))

display(EbiBlend(bggray, fg, bl_overlay)) ##FIXME: algorithm???

display(EbiBlend(bggray, fg, bl_darkenOnly))

display(EbiBlend(bggray, fg, bl_lightenOnly))

display(EbiBlend(bghue, fg, bl_hue))

display(EbiBlend(bgval, fg, bl_value))

display(EbiBlend(bgsat, fg, bl_saturation))

display(EbiBlend(bgcol, fg, bl_color))

display(EbiBlend(bgsat, fg, bl_alpha, alpha=0.5))
display(EbiBlend(bgsat, fg, bl_alpha, alpha=0.7))
display(EbiBlend(bgsat, fg, bl_alpha, alpha=0.3))
display(EbiBlend(bgsat, fg, bl_alpha, alpha=1))


