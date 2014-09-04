library(flowCore)

comp.mat = as.matrix(read.table(system.file("extdata","compdata","compmatrix",package="flowCore"),header=TRUE,skip=2,check.names=FALSE))

tube.id = function(x,d) as.numeric(gsub("060909.","",d[["$FIL"]]))

comp.fs1 = read.flowSet(path=system.file("extdata","compdata","data",package="flowCore")
#                        ,phenoData=list("Tube"=tube.id)
                        )
#phenoData(comp.fs1 )$Tube 
fsApply(comp.fs1,colMeans,use.exprs=TRUE)
fsApply(comp.fs1,each_col,range,simplify=FALSE)

# Compensate
row.names(comp.mat) = colnames(comp.mat)
comp.mat
comp.fs1 = compensate(comp.fs1,comp.mat)
fsApply(comp.fs1,colMeans,use.exprs=TRUE)

# Do some transformations
truncTrans   = truncateTransform("truncate",a=1)
linearTrans  = linearTransform("linear",a=2,b=0)

#Try a linear transform
fsApply(transform("FL1-H"=linearTrans) %on% comp.fs1,each_col,range,simplify=FALSE)
#Truncate all columns
fsApply(transform("FL1-H"=truncTrans,"FL2-H"=truncTrans,"FL3-H"=truncTrans,"FL4-H"=truncTrans) %on% comp.fs1,
	each_col,range,simplify=FALSE)

#Try to gate on a transformed value. Get ONLY side-scatter values > .3 after norming to [0,1]
normTrans = scaleTransform("norm",a=0,b=1023)
normGate  = rectangleGate("SSC-H"=c(.3,Inf))

fsApply(Subset(comp.fs1,normGate %on% transform("SSC-H"=normTrans)),each_col,range,simplify=FALSE)


# transformList
chnls <- colnames(comp.mat)
transList <- transformList(chnls, logicleTransform())
trans.fs1 <- transform(comp.fs1, transList)
fsApply(trans.fs1,colMeans,use.exprs=TRUE)

#expect the error
chnls <- c(chnls, "dummy")
transList <- transformList(chnls, logicleTransform())
trans.fs1 <- transform(comp.fs1, transList)

