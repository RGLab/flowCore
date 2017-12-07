context("Functional Tests Related to Loading/Saving flowFrames from FCS ")
library(digest)
dataPath <- "~/rglab/workspace/flowCore/misc/"
expectRes <- readRDS("~/rglab/workspace/flowCore/tests/testthat/expectResults.rds")
expectRes.new <- readRDS("~/rglab/workspace/flowCore/tests/testthat/expectRes.new.rds")
# saveRDS(expectRes.new, file = "~/rglab/workspace/flowCore/tests/testthat/expectRes.new.rds")

# test_that("Miltenyi's Macsquantify", {
#   expect_warning(
#     fr <- read.FCS(file.path(dataPath, "Miltenyi/1696_12017-04-30.0001_compatible.fcs"))
#     , "Missing the required")
#   expect_is(fr, "flowFrame")
#   
#   range(fr, type = "data")
#   range(fr)
#   expect_warning(fr <- read.FCS(file.path(dataPath, "double_precision/wishbone_thymus_panel1_rep1.fcs")), "Missing the required")
#   
#   expect_equal(nrow(fr), 250170)
#   
# })
# expectRes.new <- list()
test_that("big file", {
  expect_error(fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs")), "integer limits")
  expect_error(fr <- read.FCS(file.path(dataPath, "gigantic_file.fcs"), which.lines = c(1,1e9)), "number of collected events")
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
  expect_warning(write.FCS(fr, filename = tmp))
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
    expect_equal(expectRes[["read.FCS"]][["Blank"]], digest(fr))
    
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
      expect_equal(expectRes[["read.FCS"]][["GFP2Kfold"]], digest(fr))
      
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

test_that("test write.FCS", {
    fcsfile <- "~/rglab/workspace/flowWorkspace/wsTestSuite/Cytotrol/NHLBI/Tcell/CytoTrol_CytoTrol_1.fcs"
    fr <- read.FCS(fcsfile)
    expect_equal(keyword(fr)[["transformation"]], "applied")
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
    expect_equal(keys.new[names(keys)], keys)
    expect_equivalent(exprs(fr), exprs(fr1))
    
    #disable default linearize trans
    fr_notrans <- read.FCS(fcsfile, transformation = FALSE)
    expect_null(keyword(fr_notrans)[["transformation"]])
    #flowCore_$PnR and transformation keywords should be absent now
    #and there are should be no other difference in keywords between the two read
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr_notrans)))]
    expect_equal(length(missing.keys), 25)
    expect_true(all(grepl("(flowCore_\\$P)|(transformation)", missing.keys)))
    #any the resulted write will produce no trans related keyword r
    suppressWarnings(write.FCS(fr_notrans,tmp))
    fr1 <- read.FCS(tmp, transformation = FALSE)
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr1)))]
    expect_equal(length(missing.keys), 25)
    expect_true(all(grepl("(flowCore_\\$P)|(transformation)", missing.keys)))
    # when default linearize is enabled
    fr1 <- read.FCS(tmp)
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr1)))]
    expect_equal(length(missing.keys), 0)
    
    #transform fr
    fr.trans <- transform(fr_notrans, estimateLogicle(fr_notrans, markernames(fr_notrans)))
    expect_equal(keyword(fr.trans)[["transformation"]], "custom")
    #new keywords flowCore_$P* has been inserted
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr.trans)))]
    expect_equal(length(missing.keys), 0)
    suppressWarnings(write.FCS(fr.trans,tmp))
    #these keywords remains even disable trans when read  it back
    fr1 <- read.FCS(tmp, transformation = FALSE)
    expect_equal(keyword(fr1)[["transformation"]], "custom")
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr1)))]
    expect_equal(length(missing.keys), 0)
    #and transformation flag has no effect on read when it is already custom
    fr1 <- read.FCS(tmp)
    expect_equal(keyword(fr1)[["transformation"]], "custom")
    missing.keys <- names(keys)[which(!names(keys) %in% names(description(fr1)))]
    expect_equal(length(missing.keys), 0)
    
    
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
