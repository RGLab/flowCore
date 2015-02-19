context("Functional Tests Related to Loading/Saving flowFrames from FCS ")
library(digest)
dataPath <- "~/rglab/workspace/flowCore/misc/"
# 'FILENAME' keyword may change when file path is changed 
# so we hard code it to make the comparsion consistent in case the file is moved


test_that("test odd-bitwidth FCS", {
      fr <- read.FCS(file.path(dataPath, "Sample 2.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy" 
      expect_equal(expectRes[["read.FCS"]][["Sample2"]], digest(fr))
      
      fr <- read.FCS(file.path(dataPath, "oddbitwith/11ColorSmall.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["11ColorSmall"]], digest(fr))
      
      #too slow to test
#      fr <- read.FCS(file.path(dataPath, "oddbitwith/11ColorForOthers.fcs"))
#      keyword(fr)[["FILENAME"]] <- "setToDummy"
#      expect_equal(expectRes[["read.FCS"]][["11ColorForOthers"]], digest(fr))
      
#      fr <- read.FCS("file.path(dataPath, "oddbitwith/11ColorFull.fcs"))
#      keyword(fr)[["FILENAME"]] <- "setToDummy"
#      expect_equal(expectRes[["read.FCS"]][["11ColorFull"]], digest(fr))
    })



test_that("test other FCS", {
    fr <- read.FCS("~/rglab/workspace/QUALIFIER/misc/ITN029ST/20110125240_F06_I025.fcs")
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["ITN029ST"]], digest(fr))
    
    fr <- read.FCS(list.files(system.file("extdata", package="flowCore"),full=T)[1])
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["flowCore"]], digest(fr))
    
    fr <- read.FCS(file.path(dataPath, "Blank.FCS"))
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["Blank"]], digest(fr))
    
    expect_output(fr <- read.FCS(file.path(dataPath, "Bendall et al Cell Sample A_basal.fcs"))
                  , "dropped")
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["Bendall"]], digest(fr))              
})


test_that("test delimiter issue", {
     
      expect_error(read.FCS(file.path(dataPath, "GFP_2Kfold_011911_noKan_QA-1.fcs"))
          , "Empty keyword name detected")
      fr <- read.FCS(file.path(dataPath, "GFP_2Kfold_011911_noKan_QA-1.fcs"), emptyValue = F)
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["GFP2Kfold"]], digest(fr))
      

      fr <- read.FCS(file.path(dataPath, "RAINBOW_OK.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["RAINBOW"]], digest(fr))
      fr <- read.FCS(file.path(dataPath,"RAINBOW_OK.fcs"), emptyValue = F)
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["RAINBOWEmptyValue"]], digest(fr))               

      #\ as delimiter  with empty value
      fr <- read.FCS(file.path(dataPath, "sample_1071.001"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["sample1071"]], digest(fr))
      fr <-read.FCS(file.path(dataPath, "sample_1071.001"),emptyValue=F)
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["sample1071.double"]], digest(fr))
    })



test_that("test Beckman_Coulter_XDP issue", {
      frList <- lapply(list.files(file.path(dataPath, "Beckman_Coulter/Beckman_Coulter_XDP/"),full=T)
                        , function(thisFile){
                          fr <- read.FCS(thisFile)
                          keyword(fr)[["FILENAME"]] <- "setToDummy"
                        })
      expect_equal(expectRes[["read.FCS"]][["BeckmanCoulterXDP"]], digest(frList))
      
    })

test_that("test Beckman_Coulter $SPILLOVER keyword", {
      frList <- lapply(list.files(file.path(dataPath, "Beckman_Coulter"),full=T, pattern = ".fcs")
                    , function(thisFile){
                     fr <- read.FCS(thisFile)
                     keyword(fr)[["FILENAME"]] <- "setToDummy"
                    })
      expect_equal(expectRes[["read.FCS"]][["BeckmanCoulterSPILLOVER"]], digest(frList))
      
    })

test_that("test write.FCS", {
      
    fr <- read.FCS("~/rglab/workspace/flowWorkspace/wsTestSuite/Cytotrol/NHLBI/Tcell/CytoTrol_CytoTrol_1.fcs")
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["NHLBI"]], digest(fr))
    
    # Write to file
    tmp <- tempfile()
    write.FCS(fr,tmp)
    
    # When I read the file back in, the SPILL matrix appears to be malformed.
    fr <- read.FCS(tmp)
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["NHLBIWrite"]], digest(fr))

})


test_that("test U mode", {
  expect_error(read.FCS(file.path(dataPath, "PartecPAI/A0006980.FCS"))
              , "MODE U")
})

test_that("test pre-gated data", {
  fr <- read.FCS(file.path(dataPath, "HC002_Col1_P3.fcs"))
  keyword(fr)[["FILENAME"]] <- "setToDummy"
  expect_equal(expectRes[["read.FCS"]][["HC002_Col1_P3"]], digest(fr))
  
  fr <- read.FCS(file.path(dataPath, "HC002_Col1.fcs"))
  keyword(fr)[["FILENAME"]] <- "setToDummy"
  expect_equal(expectRes[["read.FCS"]][["HC002_Col1"]], digest(fr))

})

test_that("test flowJo exported data with offset = 99999999 and  missing the $BEGINDATA and $ENDDATA keywords ", {
      fr <- read.FCS(file.path(dataPath, "badFlowJoExport.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["badFlowJoExport"]], digest(fr))
    })