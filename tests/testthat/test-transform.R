comp.fs1 <- read.flowSet(path = system.file("extdata","compdata","data",package="flowCore"))
fr <- comp.fs1[[1]]
test_that("compensate & transform", {
      
      comp.mat <- as.matrix(read.table(system.file("extdata","compdata","compmatrix",package="flowCore"),header=TRUE,skip=2,check.names=FALSE))
      
      tube.id <- function(x,d) as.numeric(gsub("060909.","",d[["$FIL"]]))
      
      
      
      
    # Compensate with single comp
      row.names(comp.mat) <- colnames(comp.mat)
      
      comp.fs2 <- compensate(comp.fs1,comp.mat)
      
      expect_equal(fsApply(comp.fs2, colMeans, use.exprs = TRUE), expectRes[["trans.comp"]])
              
      #extra chanls
      comp.mat1 <- comp.mat[c(1:4,1), c(1:4,1)]
      colnames(comp.mat1)[5] <- "dd"
      rownames(comp.mat1)[5] <- "dd"
      expect_error(compensate(fr, comp.mat1), "not present in the flowFrame")

      
      # list
      comp <- sapply(sampleNames(comp.fs1), function(sn)comp.mat, simplify = FALSE)
      comp.fs3 <- compensate(comp.fs1,comp)
      expect_equal(fsApply(comp.fs3, colMeans, use.exprs = TRUE), expectRes[["trans.comp"]])
      
      # unmatched names
      names(comp)[1] <- "dd"
      expect_error(compensate(comp.fs1,comp), regexp = "must match")
      
      #unmatched length
      comp <- comp[1:3]
      expect_error(compensate(comp.fs1,comp), regexp = "must match")
      
      #modify comp[5]
      comp <- sapply(sampleNames(comp.fs1), function(sn)comp.mat, simplify = FALSE)
      comp[[5]][2] <- 0.001
      comp.fs4 <- compensate(comp.fs1,comp)
      expect_failure(expect_equal(fsApply(comp.fs4, colMeans, use.exprs = TRUE)
                  , expectRes[["trans.comp"]]), regexp = "8.399298e-06")
      
      #extra comp element
      comp <- sapply(sampleNames(comp.fs1), function(sn)comp.mat, simplify = FALSE)
      comp[["dd"]] <- 1:10
      comp.fs <- compensate(comp.fs1, comp)
      expect_equal(fsApply(comp.fs, colMeans, use.exprs = TRUE), expectRes[["trans.comp"]])
      
      # Do some transformations
      truncTrans <- truncateTransform("truncate",a=1)
      linearTrans <- linearTransform("linear",a=2,b=0)
      
      #linearTrans is not visible from transform method
      #since non-standard evaluation  will only evaluate its direct parent environment
      #Try a linear transform
#      expect_equal(fsApply(transform(comp.fs1, `FL1-H` = linearTrans(`FL1-H`)), each_col,range)
#                  , expectRes[["trans.linear"]])
      
      #the new lazy evaluation is still not working perfectly within test R session
      transList <- transformList("FL1-H", linearTrans)
      thisRes <- fsApply(transform(comp.fs, transList), each_col,range)
      names(dimnames(thisRes)[[2]]) <- NULL
      expect_equal(thisRes, expectRes[["trans.linear"]])
      
      
      #Truncate all columns
      transList <- transformList(c("FL1-H", "FL2-H", "FL3-H", "FL4-H"), truncTrans)
      thisRes <- fsApply(transform(comp.fs, transList),each_col,range)
      names(dimnames(thisRes)[[2]]) <- NULL
      expect_equal(thisRes, expectRes[["trans.trunc"]])
      
      expectRes[["trans.trunc"]] <- thisRes
      #Try to gate on a transformed value. Get ONLY side-scatter values > .3 after norming to [0,1]
      normTrans <- scaleTransform("norm",a=0,b=1023)
      normGate  <- rectangleGate("SSC-H"=c(.3,Inf))
      transList <- transformList("SSC-H", normTrans)
      expect_equal(fsApply(Subset(transform(comp.fs,transList), normGate),each_col,range)
                  , expectRes[["trans.normTrans"]])
      
              
      # transformList
      chnls <- colnames(comp.mat)
      transList <- transformList(chnls, logicleTransform())
      trans.fs1 <- transform(comp.fs, transList)
      expect_equal(fsApply(trans.fs1,colMeans,use.exprs=TRUE)
                  , expectRes[["trans.transformList"]])
      
      #list of transformList
      trans.list <- sapply(sampleNames(comp.fs), function(sn)transList)
      trans.fs1 <- transform(comp.fs, trans.list)
      expect_equal(fsApply(trans.fs1,colMeans,use.exprs=TRUE)
                   , expectRes[["trans.transformList"]])
      
      trans.list[[1]] <- logicleTransform()
      expect_error(trans.fs1 <- transform(comp.fs, trans.list), "a valid 'transformList'")
      
      trans.list[[1]] <- trans.list[[2]]
      names(trans.list)[1] <- "d"
      expect_error(trans.fs1 <- transform(comp.fs, trans.list), "consistent with flow data")
      
      #expect the error by giving bad channel name
      chnls <- c(chnls, "dummy")
      transList <- transformList(chnls, logicleTransform())
      expect_error(transform(comp.fs, transList), "dummy is not a variable in the flowFrame")
     
      
      
    })

test_that("hyperlogGml2", {
  
  trans <- hyperlogtGml2("FL1-H")
  trans <- eval(trans)
  res <- trans(fr)
  expect_equal(summary(res), expectRes[["hyperlogGml2"]], tol = 2e-4)
})

test_that("logicle", {
  
  trans <- logicletGml2("FL1-H", A = 2)
  trans <- eval(trans)
  res <- trans(fr)
  expect_equal(summary(res), expectRes[["logicleGml2"]], tol = 2e-4)
  
  trans <- logicleTransform(a = 2)
  raw <- exprs(fr)[,"FL1-H"]
  res <- trans(raw)
  expect_equal(summary(res), expectRes[["logicle"]], tol = 2e-4)
  
  inv <- inverseLogicleTransform(trans, "")
  expect_equal(inv(res), raw)
  
  #invalid params
  trans <- logicleTransform(a = 2, w = 2)
  expect_error(res <- trans(raw), "W = 2")
})

test_that("biexponential", {

  trans <- biexponentialTransform()
  raw <- exprs(fr)[,"FL1-H"]
  res <- trans(raw)
  expect_equal(summary(res), expectRes[["biexponential"]], tol = 2e-4)
})
