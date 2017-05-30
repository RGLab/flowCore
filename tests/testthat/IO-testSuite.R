context("Functional Tests Related to Loading/Saving flowFrames from FCS ")
library(digest)
dataPath <- "~/rglab/workspace/flowCore/misc/"
expectRes <- readRDS("~/rglab/workspace/flowCore/tests/testthat/expectResults.rds")

test_that("DATATYPE:'D'", {
  expect_warning(fr <- read.FCS(file.path(dataPath, "double_precision/wishbone_myleoid_monocyte.fcs")), "Missing the required")
  expect_is(fr, "flowFrame")
  
  expect_warning(fr <- read.FCS(file.path(dataPath, "double_precision/wishbone_thymus_panel1_rep1.fcs")), "Missing the required")
  
  expect_equal(nrow(fr), 250170)
  
})


test_that("multi data segment", {
  expect_warning(fr <- read.FCS(file.path(dataPath, "multi-datasegment.fcs")), "39 additional data")
  expect_is(fr, "flowFrame")
  
  expect_equal(nrow(fr), 1244)
  
  fr <- read.FCS(file.path(dataPath, "multi-datasegment.fcs"), dataset = 1)
  expect_equal(nrow(fr), 1244)

  fr <- read.FCS(file.path(dataPath, "multi-datasegment.fcs"), dataset = 10)
  expect_equal(nrow(fr), 955)
  
})

test_that("FCS with both SPILL and $SPILLOVER present", {

  fr <- read.FCS(file.path(dataPath, "/example-spill-spillover.fcs"))
  expect_equal(keyword(fr)[["SPILL"]], keyword(fr)[["$SPILLOVER"]])
  tmp <- spillover(fr)
  expect_equal(tmp[["SPILL"]], tmp[["$SPILLOVER"]])
  tmp <- tempfile()
  expect_warning(write.FCS(fr, filename = tmp))
  fr <- read.FCS(tmp)  
  expect_equal(keyword(fr)[["SPILL"]], keyword(fr)[["$SPILLOVER"]])
  
})


test_that("test uint_64 + diverse bitwidths + missing $NEXTDATA: '*' ", {
      
      expect_warning(fr <- read.FCS(file.path(dataPath, "uint_64.lxb")))
      
    })


#TODO: investigate why the results is no longer consistent with the archived summary
test_that("mixed endian", {
  fr <- read.FCS(file.path(dataPath, "mixedEndian.fcs"))
  expect_is(fr, "flowFrame")
  # expect_equal(summary(fr), expectRes[["read.FCS"]][["mixedEndian"]])
  
  
})

# 'FILENAME' keyword may change when file path is changed 
# so we hard code it to make the comparsion consistent in case the file is moved
test_that("test special delimiter character: '*' ", {
  
  expect_warning(fr <- read.FCS(file.path(dataPath, "multi_data_segment.LMD")), "1 additional data")
  expect_equal(summary(fr), expectRes[["read.FCS"]][["multi_data_segment"]], tolerance = 0.08)
  
})



test_that("test special delimiter character: '*' ", {
  
    fr <- read.FCS(file.path(dataPath, "specialDelimiter.fcs"))
    expect_equal(summary(fr), expectRes[["read.FCS"]][["specialDelimiter"]], tolerance = 0.001)
  
})


test_that("test flowJo exported data with missing some of PnR keywords ", {
expect_warning(expect_error(fr <- read.FCS(file.path(dataPath, "missing_PnR_flowJoExport.fcs"))
               , "not contained")
               , "Missing the required \\$BEGINDATA keyword")
  
})


test_that("test in consistent datastart between header and TEXT", {
      expect_error(fr <- read.FCS(file.path(dataPath, "Accuri-C6", "Accuri - C6 - A02 Spherotech 8 Peak Beads.fcs"), emptyValue = FALSE)
                   , "HEADER and the TEXT")
     
     expect_warning(fr <- read.FCS(file.path(dataPath, "Accuri-C6", "Accuri - C6 - A02 Spherotech 8 Peak Beads.fcs"), emptyValue = FALSE, ignore.text.offset = TRUE)
                    , "HEADER and the TEXT")
     expect_equal(nrow(fr), 60661)      
     expect_equal(summary(fr), expectRes[["read.FCS"]][["Accuri-C6"]], tolerance = 0.001)
     
     expect_warning(fs <- read.flowSet(file.path(dataPath, "Accuri-C6", "Accuri - C6 - A02 Spherotech 8 Peak Beads.fcs"), emptyValue = FALSE, ignore.text.offset = TRUE)
                    , "HEADER and the TEXT")
     expect_equal(nrow(fs[[1]]), 60661)

   })

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
    fr <- read.FCS(file.path(dataPath, "20110125240_F06_I025.fcs"))
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


# latest R no longer permit overflowed coersion by as.integer
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
    fcsfile <- "~/rglab/workspace/flowWorkspace/wsTestSuite/Cytotrol/NHLBI/Tcell/CytoTrol_CytoTrol_1.fcs"
    fr <- read.FCS(fcsfile)
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["NHLBI"]], digest(fr))
    
    # Write to file
    tmp <- tempfile()
    suppressWarnings(write.FCS(fr,tmp))
    
    # When I read the file back in, the SPILL matrix appears to be malformed.
    fr1 <- read.FCS(tmp)
    keys <- description(fr)
    keys[["$TOT"]] <- trimws(keys[["$TOT"]])
    keys[c("$BEGINDATA", "$ENDDATA")] <- NULL
    keys.new <- description(fr1)
    keys.new[["FILENAME"]] <- "setToDummy"
    #update flowCore range in keyword
    # rng <- range(fr)
    # for(i in 1:ncol(fr))
    # { 
    #   keys[[sprintf("flowCore_$P%sRmin", i)]] <- as.character(rng[1,i])
    #   keys[[sprintf("flowCore_$P%sRmax", i)]] <- as.character(rng[2,i])
    # }
    expect_equal(keys.new[names(keys)], keys)
    expect_equivalent(exprs(fr), exprs(fr1))
    
    # test delimiter(\) escaping 
    description(fr)[["$DATE"]] <- "05\\JUN\\2012"
    suppressWarnings(write.FCS(fr,tmp))
    fr1 <- read.FCS(tmp, emptyValue = F)
    keys.new <- description(fr1)
    keys.new[["FILENAME"]] <- "setToDummy"
    expect_equal(keys.new[["$DATE"]], "05\\\\JUN\\\\2012")
    keys.new[["$DATE"]] <- keys[["$DATE"]]
    expect_equal(keys.new[names(keys)], keys)
    expect_equivalent(exprs(fr), exprs(fr1))
    
    # write it again to see if the existing double delimiter is handled properly
    suppressWarnings(write.FCS(fr1,tmp))
    fr1 <- read.FCS(tmp, emptyValue = F)
    keys.new <- description(fr1)
    keys.new[["FILENAME"]] <- "setToDummy"
    expect_equal(keys.new[["$DATE"]], "05\\\\JUN\\\\2012")
    keys.new[["$DATE"]] <- keys[["$DATE"]]
    expect_equal(keys.new[names(keys)], keys)
    expect_equivalent(exprs(fr), exprs(fr1))

    #test other delimiter
    suppressWarnings(write.FCS(fr,tmp, delimiter = ";"))
    fr1 <- read.FCS(tmp, emptyValue = F)
    keys.new <- description(fr1)
    keys.new[["FILENAME"]] <- "setToDummy"
    expect_equal(keys.new[["$DATE"]], "05\\JUN\\2012")
    keys.new[["$DATE"]] <- keys[["$DATE"]]
    expect_equal(keys.new[names(keys)], keys)
    expect_equivalent(exprs(fr), exprs(fr1))
    
    #when colmn.pattern is used to subset channels in read.FCS
    #make sure the id in $Pn is set properly in write.FCS
    fr_sub <- read.FCS(fcsfile, column.pattern = '-A')
    tmp <- tempfile()
    suppressWarnings(write.FCS(fr_sub , filename = tmp))
    fr1 <- read.FCS(tmp)
    expect_equal(pData(parameters(fr_sub))[["name"]], pData(parameters(fr1))[["name"]], check.attributes = FALSE)
    expect_equal(pData(parameters(fr_sub))[["desc"]], pData(parameters(fr1))[["desc"]], check.attributes = FALSE)
})

test_that("write.flowSet", {
  
  data(GvHD)
  foo <- GvHD[1:2]
  
  
  ## now write out into  files
  outDir <- tempfile()
  suppressWarnings(write.flowSet(foo, outDir))
  expect_equal(dir(outDir), c("annotation.txt", "s5a01.fcs", "s5a02.fcs"))
  
  outDir <- tempfile()
  suppressWarnings(write.flowSet(foo, outDir, filename = c("a")))
  expect_equal(dir(outDir), c("1_a.fcs", "2_a.fcs", "annotation.txt"))
  
  outDir <- tempfile()
  suppressWarnings(write.flowSet(foo, outDir, filename = c("a", "b")))
  expect_equal(dir(outDir), c("a.fcs", "annotation.txt", "b.fcs"))
  
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

test_that("test integer overflow issue", {
      fr <- read.FCS(file.path(dataPath, "intOverFlow.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["intOverFlow"]], digest(fr))
      fr <- read.FCS(file.path(dataPath,"/Beckman_Coulter/MoFlo Astrios EQ 9C bis all.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["MoFlo EQ 9C"]],  digest(fr))
    })

test_that("test diverse Bitwidths", {
      fr <- read.FCS(file.path(dataPath, "diverseBitwidths.fcs"))
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["diverseBitwidths"]], digest(fr))
    })
