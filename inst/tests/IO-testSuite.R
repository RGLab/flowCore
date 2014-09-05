context("Unit Tests Related to Loading/Saving flowFrames from FCS (or other formats")

test_that("test 9-bitwidth FCS", {
      fr3 <- read.FCS("~/rglab/workspace/flowCore/misc/Sample 2.fcs")      
    })

time1<-Sys.time()
fr1<-read.FCS("~/rglab/workspace/QUALIFIER/misc/ITN029ST/20110125240_F06_I025.fcs")
Sys.time()-time1
time1<-Sys.time()
fr2<-read.FCS("~/rglab/workspace/QUALIFIER/misc/ITN029ST/20110125240_F06_I025.fcs")
Sys.time()-time1
identical(fr1,fr2)

fs1<-read.flowSet(list.files(system.file("extdata", package="flowCore"),full=T)[1:3])
save(fs1,file="temp.rda")
fs2<-read.flowSet(list.files(system.file("extdata", package="flowCore"),full=T)[1:3])
load("/loc/no-backup/mike/HIPC/output/temp.rda")
for(i in 1:3)
  print(identical(fs1[[i]],fs2[[i]]))

###test delimiter issue
description(read.FCS("~/rglab/workspace/flowCore/misc/GFP_2Kfold_011911_noKan_QA-1.fcs",emptyValue=F))
read.flowSet("~/rglab/workspace/flowCore/misc/GFP_2Kfold_011911_noKan_QA-1.fcs",emptyValue=F)
description(read.FCS("~/rglab/workspace/QUALIFIER/misc/ITN029ST/20110125240_F06_I025.fcs"))
read.flowSet("~/rglab/workspace/QUALIFIER/misc/ITN029ST/20110125240_F06_I025.fcs")
description(read.FCS(system.file("extdata","0877408774.B08",package="flowCore")))

description(read.FCS("~/rglab/workspace/flowCore/misc/RAINBOW_OK.fcs"))

description(read.FCS("~/rglab/workspace/flowCore/misc/sample_1071.001",emptyValue=T))



read.FCS("~/rglab/workspace/flowCore/misc/Blank.FCS")
read.flowSet("~/rglab/workspace/flowCore/misc/Blank.FCS")

inFile <- system.file("extdata", "0877408774.B08", package="flowCore")
foo <- read.FCS(inFile, transform=FALSE)
outFile <- file.path(tempdir(), "foo.fcs")
## now write out into a file
write.FCS(foo, outFile)
bar <- read.FCS(outFile, transform=FALSE)
description(bar)

data(GvHD)
foo <- GvHD[1:5]
outDir <- file.path(tempdir(), "foo")

# now write out into  files
write.flowSet(foo, outDir)
description(foo[[1]])

library(tools)
setwd("flowCore/inst/doc")
Sweave("HowTo-flowCore.Rnw")
texi2dvi("HowTo-flowCore.tex",pdf=TRUE)



#test Beckman_Coulter_XDP error
library(flowCore)
fs<-read.flowSet(list.files("flowCore/misc/Beckman_Coulter_XDP/",full=T)[1])
fr<-read.FCS(file.path(dirname(list.files("flowCore/misc/Beckman_Coulter_XDP/",full=T)[1]),"tmp.fcs"))

##

a = read.FCS("/loc/no-backup/mike/HIPC/data/Cytotrol/NHLBI/Tcell/CytoTrol_CytoTrol_1.fcs")

# SPILL matrix appears to be correct
print(keyword(a,"SPILL"))

# Write to file
write.FCS(a,"/loc/no-backup/mike/HIPC/output/b.fcs")

# When I read the file back in, the SPILL matrix appears to be malformed.
b = read.FCS("/loc/no-backup/mike/HIPC/output/b.fcs")
print(keyword(b,"SPILL"))

#
fr<-read.FCS("~/rglab/workspace/flowCore/misc/PartecPAI/A0006980.FCS")

##pre-gated data
fr<-read.FCS("/home/wjiang2/rglab/workspace/flowCore/misc/HC002_Col1_P3.fcs")
fr<-read.FCS("/home/wjiang2/rglab/workspace/flowCore/misc/HC002_Col1.fcs")
keyword(fr)$SPIL

