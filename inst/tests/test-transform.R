comp.fs1 <- read.flowSet(path = system.file("extdata","compdata","data",package="flowCore")
                         #                        ,phenoData=list("Tube"=tube.id)
                              )
test_that("transform", {
      
      comp.mat <- as.matrix(read.table(system.file("extdata","compdata","compmatrix",package="flowCore"),header=TRUE,skip=2,check.names=FALSE))
      
      tube.id <- function(x,d) as.numeric(gsub("060909.","",d[["$FIL"]]))
      
      
      
      
    # Compensate
      row.names(comp.mat) <- colnames(comp.mat)
      
      comp.fs1 <- compensate(comp.fs1,comp.mat)
      
      expect_equal(fsApply(comp.fs1, colMeans, use.exprs = TRUE)
                  , expectRes[["trans.comp"]])
      
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
      thisRes <- fsApply(transform(comp.fs1, transList), each_col,range)
      names(dimnames(thisRes)[[2]]) <- NULL
      expect_equal(thisRes, expectRes[["trans.linear"]])
      
      
      #Truncate all columns
      transList <- transformList(c("FL1-H", "FL2-H", "FL3-H", "FL4-H"), truncTrans)
      thisRes <- fsApply(transform(comp.fs1, transList),each_col,range)
      names(dimnames(thisRes)[[2]]) <- NULL
      expect_equal(thisRes, expectRes[["trans.trunc"]])
      
      expectRes[["trans.trunc"]] <- thisRes
      #Try to gate on a transformed value. Get ONLY side-scatter values > .3 after norming to [0,1]
      normTrans <- scaleTransform("norm",a=0,b=1023)
      normGate  <- rectangleGate("SSC-H"=c(.3,Inf))
      transList <- transformList("SSC-H", normTrans)
      expect_equal(fsApply(Subset(transform(comp.fs1,transList), normGate),each_col,range)
                  , expectRes[["trans.normTrans"]])
      
              
      # transformList
      chnls <- colnames(comp.mat)
      transList <- transformList(chnls, logicleTransform())
      trans.fs1 <- transform(comp.fs1, transList)
      expect_equal(fsApply(trans.fs1,colMeans,use.exprs=TRUE)
                  , expectRes[["trans.transformList"]])
      
      #expect the error by giving bad channel name
      chnls <- c(chnls, "dummy")
      transList <- transformList(chnls, logicleTransform())
      expect_error(transform(comp.fs1, transList), "dummy is not a variable in the flowFrame")
      
      
    })

test_that("hyperlogGml2", {
  fr <- comp.fs1[[1]]
  trans <- hyperlogtGml2("FL1-H")
  trans <- eval(trans)
  res <- trans(fr)
  expect_equal(summary(res), expectRes[["hyperlogGml2"]])
})

test_that("logicle", {
  fr <- comp.fs1[[1]]
  trans <- logicletGml2("FL1-H", A = 2)
  trans <- eval(trans)
  res <- trans(fr)
  expect_equal(summary(res), expectRes[["logicleGml2"]])
  
  trans <- logicleTransform(a = 2)
  raw <- exprs(fr)[,"FL1-H"]
  res <- trans(raw)
  expect_equal(summary(res), expectRes[["logicle"]])
  
  inv <- inverseLogicleTransform("", trans)
  expect_equal(inv(res), raw)
})
