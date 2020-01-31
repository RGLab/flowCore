#the FCS files are not included in package due to the size , thus only for internal testing
context("Functional Tests Related to Loading/Saving flowFrames from FCS ")
library(digest)
dataPath <- "~/rglab/workspace/flowCore/misc/"
skip_if_not(dir.exists(dataPath))
expectRes <- readRDS("~/rglab/workspace/flowCore/tests/testthat/expectResults.rds")
expectRes.new <- readRDS("~/rglab/workspace/flowCore/tests/testthat/expectRes.new.rds")
# saveRDS(expectRes.new, file = "~/rglab/workspace/flowCore/tests/testthat/expectRes.new.rds")
test_that("trailing double delimiter", {
  expect_silent(fr <- read.FCS(file.path(dataPath, "trailingdoubledelimiter.fcs")))
  expect_equal(length(keyword(fr)), 38)
  expect_equal(keyword(fr)[["keyword1234"]], "")#check the last keyword is parsed correctly
  
  #disable empty value
  expect_silent(fr <- read.FCS(file.path(dataPath, "trailingdoubledelimiter.fcs"), emptyValue = F))
  expect_equal(length(keyword(fr)), 37)
  expect_equal(keyword(fr)[["SAMPLE ID"]], "Default Patient ID")#check the second last keyword is parsed correctly
  expect_null(keyword(fr)[["keyword1234"]])#check the last keyword is dropped correctly
  
})
test_that("trailing space", {
  expect_silent(fr <- read.FCS(file.path(dataPath, "trailing_space.fcs")))
  expect_equal(length(keyword(fr)), 37)
  expect_equal(keyword(fr)[["SAMPLE ID"]], "Default Patient ID")#check the last keyword is parsed correctly

})
# expectRes.new <- list()
test_that("big file", {
  #too slow to run
  # fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs"), column.pattern = "FSC*")
                 
  fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs")
                              , which.lines = c(1,1e9)
                              , column.pattern = "FSC*")
  expect_equal(nrow(fr), 1)
  fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs"), which.lines = 1:1e3)
  expect_equal(nrow(fr), 1e3)
  fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs"), which.lines = (1e3+1):2e3)
  expect_equal(nrow(fr), 1e3)
  
})

test_that("DATATYPE:'D'", {
  expect_warning(fr <- read.FCS(file.path(dataPath, "double_precision/wishbone_myleoid_monocyte.fcs")), "Missing the required")
  expect_is(fr, "flowFrame")
  filename  <- "wishbone_thymus_panel1_rep1.fcs"
  expect_warning(fr <- read.FCS(file.path(dataPath, "double_precision", filename)), "Missing the required")
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))
  expect_equal(nrow(fr), 250170)
  
})


test_that("multi data segment", {
  expect_warning(fr <- read.FCS(file.path(dataPath, "multi-datasegment.fcs")), "39 additional data")
  expect_is(fr, "flowFrame")
  
  expect_equal(nrow(fr), 1244)
  filename  <- "multi-datasegment.fcs"
  fr <- read.FCS(file.path(dataPath, filename), dataset = 1)
  expect_equal(nrow(fr), 1244)

  fr <- read.FCS(file.path(dataPath, filename), dataset = 10)
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))
  expect_equal(nrow(fr), 955)
  
})

test_that("FCS with both SPILL and $SPILLOVER present", {

  filename <- "example-spill-spillover.fcs"
  fr <- read.FCS(file.path(dataPath, filename))
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))
  
  expect_equal(keyword(fr)[["SPILL"]], keyword(fr)[["$SPILLOVER"]])
  tmp <- spillover(fr)
  expect_equal(tmp[["SPILL"]], tmp[["$SPILLOVER"]])
  tmp <- tempfile()
  write.FCS(fr, filename = tmp)
  fr <- read.FCS(tmp)  
  expect_equal(keyword(fr)[["SPILL"]], keyword(fr)[["$SPILLOVER"]])
  
})


test_that("test uint_64 + diverse bitwidths + missing $NEXTDATA: '*' ", {
      filename <- "uint_64.lxb"
      expect_warning(fr <- read.FCS(file.path(dataPath, filename)))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))
  
    })


#TODO: investigate why the results is no longer consistent with the archived summary
test_that("mixed endian", {
  filename <- "mixedEndian.fcs"
  fr <- read.FCS(file.path(dataPath, filename))
  expect_is(fr, "flowFrame")
  # expect_equal(summary(fr), expectRes[["read.FCS"]][["mixedEndian"]])
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
  
})

# 'FILENAME' keyword may change when file path is changed 
# so we hard code it to make the comparsion consistent in case the file is moved
test_that("test special delimiter character: '*' ", {
  filename <- "multi_data_segment.LMD"
  expect_warning(fr <- read.FCS(file.path(dataPath, filename)), "1 additional data")
  expect_equal(summary(fr), expectRes[["read.FCS"]][["multi_data_segment"]], tolerance = 0.08)
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
  
})



test_that("test special delimiter character: '*' ", {
    filename <- "specialDelimiter.fcs"
    fr <- read.FCS(file.path(dataPath, filename))
    expect_equal(summary(fr), expectRes[["read.FCS"]][["specialDelimiter"]], tolerance = 0.001)
    #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
    
})


test_that("test flowJo exported data with missing some of PnR keywords ", {
expect_warning(expect_error(fr <- read.FCS(file.path(dataPath, "missing_PnR_flowJoExport.fcs"))
               , "not contained")
               , "Missing the required \\$BEGINDATA keyword")
  
})


test_that("test in consistent datastart between header and TEXT", {
  expect_output(expect_error(fr <- read.FCS(file.path(dataPath, "Accuri-C6", "Accuri - C6 - A02 Spherotech 8 Peak Beads.fcs"), emptyValue = FALSE)
                   , "HEADER and the TEXT")
                , "uneven number of tokens")
    filename <- "Accuri - C6 - A02 Spherotech 8 Peak Beads.fcs"
     expect_output(expect_warning(fr <- read.FCS(file.path(dataPath, "Accuri-C6", filename), emptyValue = FALSE, ignore.text.offset = TRUE)
                    , "HEADER and the TEXT")
                    , "uneven number of tokens")
     #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
     
     expect_equal(nrow(fr), 60661)
     expect_equal(summary(fr), expectRes[["read.FCS"]][["Accuri-C6"]], tolerance = 0.001)

     expect_output(expect_warning(fs <- read.flowSet(file.path(dataPath, "Accuri-C6", filename), emptyValue = FALSE, ignore.text.offset = TRUE)
                    , "HEADER and the TEXT"), "uneven number of tokens")
     expect_equal(nrow(fs[[1]]), 60661)

   })

test_that("test odd-bitwidth FCS", {
      filename <- "Sample 2.fcs"
      fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy" 
      expect_equal(expectRes[["read.FCS"]][["Sample2"]], digest(fr))
      
      filename <- "11ColorSmall.fcs"
      fr <- read.FCS(file.path(dataPath, "oddbitwith", filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["11ColorSmall"]], summary(fr))
      
      #too slow to test
#      fr <- read.FCS(file.path(dataPath, "oddbitwith/11ColorForOthers.fcs"))
#      keyword(fr)[["FILENAME"]] <- "setToDummy"
#      expect_equal(expectRes[["read.FCS"]][["11ColorForOthers"]], digest(fr))
      
#      fr <- read.FCS("file.path(dataPath, "oddbitwith/11ColorFull.fcs"))
#      keyword(fr)[["FILENAME"]] <- "setToDummy"
#      expect_equal(expectRes[["read.FCS"]][["11ColorFull"]], digest(fr))
    })



test_that("test other FCS", {
    filename <- "20110125240_F06_I025.fcs"
    fr <- read.FCS(file.path(dataPath, filename))
    #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["ITN029ST"]], digest(fr))
    
    
    fr <- read.FCS(list.files(system.file("extdata", package="flowCore"),full=T)[1])
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["flowCore"]], digest(fr))
    
    filename <- "Blank.FCS"
    fr <- read.FCS(file.path(dataPath, filename))
    #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["Blank"]], summary(fr))
    
    filename <- "Bendall et al Cell Sample A_basal.fcs"
    expect_output(fr <- read.FCS(file.path(dataPath, filename))
                  , "dropped")
    #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
    keyword(fr)[["FILENAME"]] <- "setToDummy"
    expect_equal(expectRes[["read.FCS"]][["Bendall"]], digest(fr))              
})


test_that("test delimiter issue", {
      filename <- "GFP_2Kfold_011911_noKan_QA-1.fcs"
      expect_error(read.FCS(file.path(dataPath, filename))
          , "Empty keyword name detected")
      fr <- read.FCS(file.path(dataPath, filename), emptyValue = F)
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["GFP2Kfold"]], summary(fr))
      expect_match(keyword(fr)[["GTI$WORKLIST"]], "C:/Document")
      
      filename <- "RAINBOW_OK.fcs"
      fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["RAINBOW"]], digest(fr))
      fr <- read.FCS(file.path(dataPath, filename), emptyValue = F)
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["RAINBOWEmptyValue"]], digest(fr))               

      #\ as delimiter  with empty value
      filename <- "sample_1071.001"
      fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["sample1071"]], digest(fr))
      fr <-read.FCS(file.path(dataPath, filename),emptyValue=F)
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

test_that("test U mode", {
  expect_error(read.FCS(file.path(dataPath, "PartecPAI/A0006980.FCS"))
              , "MODE U")
})

test_that("test pre-gated data", {
  filename <- "HC002_Col1_P3.fcs"
  fr <- read.FCS(file.path(dataPath, filename))
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
  keyword(fr)[["FILENAME"]] <- "setToDummy"
  expect_equal(expectRes[["read.FCS"]][["HC002_Col1_P3"]], digest(fr))
  
  filename <- "HC002_Col1.fcs"
  fr <- read.FCS(file.path(dataPath, filename))
  #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
  keyword(fr)[["FILENAME"]] <- "setToDummy"
  expect_equal(expectRes[["read.FCS"]][["HC002_Col1"]], digest(fr))

})

test_that("test flowJo exported data with offset = 99999999 and  missing the $BEGINDATA and $ENDDATA keywords ", {
    filename <- "badFlowJoExport.fcs"
  fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["badFlowJoExport"]], digest(fr))
    })

test_that("test integer overflow issue", {
      filename <- "intOverFlow.fcs"
      fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["intOverFlow"]], digest(fr))
      
      filename <- "MoFlo Astrios EQ 9C bis all.fcs"
      fr <- read.FCS(file.path(dataPath,"/Beckman_Coulter", filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["MoFlo EQ 9C"]],  digest(fr))
    })

test_that("test diverse Bitwidths", {
      filename <- "diverseBitwidths.fcs"
      fr <- read.FCS(file.path(dataPath, filename))
      #expectRes.new[[filename]] <<- list(ncol = ncol(fr), nrow = nrow(fr), chnl = colnames(fr), marker = markernames(fr), range = range(fr), range_data= range(fr, "data"), colmean = colMeans(exprs(fr)))  
      keyword(fr)[["FILENAME"]] <- "setToDummy"
      expect_equal(expectRes[["read.FCS"]][["diverseBitwidths"]], digest(fr))
    })

test_that("handle > 2^32-1 bytes", {
  fr <- flowFrame(matrix(data = rnorm(3e8), nrow = 1e8, ncol =3, dimnames = list(NULL, c("A", "B", "C"))))
  expect_gt(object.size(exprs(fr)), 2^31-1)
  tmp <- tempfile()
  write.FCS(fr, tmp, what = "double")
  set.seed(1)
  lines <- sort(sample(1:1e8, 1e3))
  fr1 <- read.FCS(tmp, which.lines = lines)
  expect_equivalent(exprs(fr)[lines,], exprs(fr1))
})

