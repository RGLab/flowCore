test_that("transform", {
      
      comp.mat <- as.matrix(read.table(system.file("extdata","compdata","compmatrix",package="flowCore"),header=TRUE,skip=2,check.names=FALSE))
      
      tube.id <- function(x,d) as.numeric(gsub("060909.","",d[["$FIL"]]))
      
      comp.fs1 <- read.flowSet(path = system.file("extdata","compdata","data",package="flowCore")
#                        ,phenoData=list("Tube"=tube.id)
                  )
      
      
    # Compensate
      row.names(comp.mat) <- colnames(comp.mat)
      
      comp.fs1 <- compensate(comp.fs1,comp.mat)
      
      expect_equal(fsApply(comp.fs1, colMeans, use.exprs = TRUE)
                  , expectRes[["trans.comp"]])
      
      # Do some transformations
      truncTrans <- truncateTransform("truncate",a=1)
      linearTrans <- linearTransform("linear",a=2,b=0)
      
      #Try a linear transform
      
      expect_equal(fsApply(transform(comp.fs1, `FL1-H` = linearTrans(`FL1-H`)), each_col,range)
                  , expectRes[["trans.linear"]])
      #Truncate all columns
      expect_equal(fsApply(transform(comp.fs1, `FL1-H`=truncTrans(`FL1-H`)
                                                ,`FL2-H`=truncTrans(`FL2-H`)
                                                ,`FL3-H`=truncTrans(`FL3-H`)
                                                ,`FL4-H`=truncTrans(`FL4-H`))
                            ,each_col,range)
                    , expectRes[["trans.trunc"]])
      
 
      #Try to gate on a transformed value. Get ONLY side-scatter values > .3 after norming to [0,1]
      normTrans <- scaleTransform("norm",a=0,b=1023)
      normGate  <- rectangleGate("SSC-H"=c(.3,Inf))
      
      expect_equal(fsApply(Subset(transform(comp.fs1, `SSC-H` = normTrans(`SSC-H`)), normGate),each_col,range)
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

