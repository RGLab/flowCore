source("smoothScatterMult.R")
source("blendFuns.R")

## some not very clever sample data
xy1 <- cbind(rnorm(10000), rnorm(10000, sd=2))
xy2 <- cbind(rnorm(10000,2, sd=1.3), rnorm(10000,4))
xy <- rbind(xy1, xy2)
fact <- factor(rep(0:1, each=nrow(xy1)))

## old smoothScatter plot
smoothScatter(xy)

##########################################################
## smoothScatterMult with different blending algorithms
##########################################################
## alpha blending (alpha of foreground is 0.5
smoothScatterMult(xy, fact=fact, blendFun=bl_alpha, bfArgs=list(alpha=0.5))

## addition of rgb values --> blows out of color space
smoothScatterMult(xy, fact=fact, blendFun=bl_addition)

## substraction of rgb values --> blows out of color space
smoothScatterMult(xy, fact=fact, blendFun=bl_substraction)

## absolute difference of rgb values --> cool!!!
smoothScatterMult(xy, fact=fact, blendFun=bl_difference)

## multiplication of rgb values --> can create black
smoothScatterMult(xy, fact=fact, blendFun=bl_burn)

## division of rgb values --> can create white
smoothScatterMult(xy, fact=fact, blendFun=bl_dodge)

## 'inverse multiplication' of rgb values --> can create white
smoothScatterMult(xy, fact=fact, blendFun=bl_screen)

## combines screen and burn --> algorithm might be wrong
smoothScatterMult(xy, fact=fact, blendFun=bl_overlay)

## paralel min of rgb values --> can create white
smoothScatterMult(xy, fact=fact, blendFun=bl_lightenOnly)

## paralel max of rgb values --> can create white
smoothScatterMult(xy, fact=fact, blendFun=bl_darkenOnly)

## only hue of foreground --> what the hell???
smoothScatterMult(xy, fact=fact, blendFun=bl_hue)

## only value of foreground --> can create white
smoothScatterMult(xy, fact=fact, blendFun=bl_value)

## only saturation of foreground --> outch!!!
smoothScatterMult(xy, fact=fact, blendFun=bl_saturation)

## luminance of background--> algorithm wrong???
smoothScatterMult(xy, fact=fact, blendFun=bl_color)


##########################################################
## Extending and combining blending algorithms:
## One could combine several of the blending modes or
## do even more complex stuff.
##########################################################
## invert the colors of differential blending 
bl_invdiff <- function(fg, bg)
  invert(bl_difference(fg,bg))
smoothScatterMult(xy, fact=fact, blendFun=bl_invdiff)

## spread the values of alpha blending over the whole color range
bl_alphaspread <- function(fg, bg)
  spread(bl_alpha(fg,bg))
smoothScatterMult(xy, fact=fact, blendFun=bl_alphaspread)


##########################################################
## Three different populations: A color mess...
##########################################################
xy3 <- cbind(rnorm(10000, -3, sd=2), rnorm(10000,1, sd=0.8))
xy2 <- rbind(xy1, xy2, xy3)
fact <- factor(rep(0:2, each=nrow(xy1)))
smoothScatter(xy2)
smoothScatterMult(xy2, fac=fact)
