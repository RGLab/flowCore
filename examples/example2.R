## compensation.R
## 
## Read in the plates of compensation tubes, compensate and be sure it is done correctly.
##
## Author: haaland
###############################################################################
library(flowCore)
library(geneplotter)
## lipo188
## get the compensation matrix
comp.mat = as.matrix(read.table("./inst/extdata/compdata/compmatrix",header=TRUE,skip=2))
colnames(comp.mat)
#[1] "FL1.H" "FL2.H" "FL3.H" "FL4.H"
colnames(comp.mat) = c("FL1-H","FL2-H","FL3-H","FL4-H")
## read in the compensation tubes
## I think the default linearlization is the correct thing
comp.fset1 = read.flowSet(path="./inst/extdata/compdata/data/")
ap = phenoData(comp.fset1)
pData(ap)$Tube = 1:5
phenoData(comp.fset1) = ap
pData(ap)
eapply(comp.fset1@frames,function(x) apply(x@exprs,2,range))

## compensate the linearized data
comp.fset1b = compensate(comp.fset1,comp.mat)

eapply(comp.fset1b@frames,function(x) apply(x@exprs,2,range))
truncTrans = truncateTransform("truncate",a=1)
linearTrans = linearTransform(transformationId="linearTransformation",a=2,b=0)
transform(comp.fset1b[[1]],`FL1-H`=linearTrans(`FL1-H`))
apply(comp.fset1b[[1]]@exprs,2,range)
apply(transform(comp.fset1b[[1]],`FL1-H`=linearTrans(`FL1-H`))@exprs,2,range)
apply(transform(comp.fset1b[[1]],`FL1-H`=truncTrans(`FL1-H`))@exprs,2,range)

## truncate values
transform(comp.fset1b[[1]],`FL1-H`=truncTrans(`FL1-H`),`FL2-H`=truncTrans(`FL2-H`),
	`FL3-H`=truncTrans(`FL3-H`),`FL4-H`=truncTrans(`FL4-H`))


comp.fset1c = transform(comp.fset1b,`FL1-H`=truncTrans(`FL1-H`),`FL2-H`=truncTrans(`FL2-H`),
	`FL3-H`=truncTrans(`FL3-H`),`FL4-H`=truncTrans(`FL4-H`))


comp.fset1c = fsApply(comp.fset1b,transform,`FL1-H`=truncTrans(`FL1-H`),`FL2-H`=truncTrans(`FL2-H`),
	`FL3-H`=truncTrans(`FL3-H`),`FL4-H`=truncTrans(`FL4-H`))

comp.fset2b = compensate(comp.fset2,comp088.mat)
eapply(comp.fset2b@frames,function(x) apply(x@exprs,2,range))

plotFour(comp.fsetb[[1]],main="Compensated",logx=FALSE,logy=FALSE)
flowPlot(comp.fsetb


pdf("plots/compensated-scaled.pdf")
eapply(comp.fsetb@frames, function(x) plotFour(x,main="Compensated",logx=FALSE,logy=FALSE))
dev.off()



plotFour(comp1,main="Tube 001: Uncompensated")
comp1.comp = compensate(comp1,comp088.mat)
apply(comp1.comp@exprs,2,range)
plotFour(comp1,main="Tube 001: Compensated on Channel Values")
comp1b = read.FCS("data/lipo088/Comp Tubes/060909.001",transformation="linearize")
plotFour(comp1b,"Tube 001: Linearized Values")
comp1b.comp = compensate(comp1b,comp088.mat)
apply(comp1b.comp@exprs,2,range)



plotFour(comp1b.comp,main="Tube 001: Compensated on Linearized Values")


