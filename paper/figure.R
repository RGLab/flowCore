## Loading library
library("flowCore")
library("flowQ")
library("flowViz")

## Read a batch of FCS files (find another data set with time channel - Rituximab to change from the GvHD?)
fs = read.flowSet(path = system.file("extdata", package="flowCore"), pattern="0877408774")

fs

## Tranform

## Assess Quality (using fs give me en error) 
qaReport(GvHD[1:3], c("qaProcess.timeline",
                      "qaProcess.timeflow",
                      "qaProcess.cellnumber"))

## Gating 
?polygonGate

## Defining the gate
sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
colnames(sqrcut) <- c("FSC-H","SSC-H")
pg <- polygonGate(filterId="nonDebris", .gate= sqrcut)
pg

fres <- filter(GvHD[1:3], pg)
fres
summary(fres)

## Filtering
## Define filter approach
n2f <- norm2Filter("FSC-H", "SSC-H", filterId="myNorm2Filter")

## Filtering using curv2Filter
fres <- filter(GvHD[1:7], n2f)
fres
summary(fres)

## Visualizing
pdf("n2f-output.pdf")
print(xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD,
            filter=n2f, subset=Patient==5))
dev.off()

## Select subpobulation for futher analysis
nonDebris= Subset(GvHD[1:3], fres)
