test_that("filter", {
      set.seed(123)
      ## create the data to be used
      ## these are three interesting wells from a BD FACS CAP(TM) plate
      ## with PBMS (perpheral blood monocytes) on the plate
      b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"),transformation="scale")
      e07 = read.FCS(system.file("extdata","0877408774.E07",package="flowCore"),transformation="scale")
      f06 = read.FCS(system.file("extdata","0877408774.F06",package="flowCore"),transformation="scale")
      
      filter1 = rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))

      b08.result1 = filter(b08,filter1)
      expect_is(b08.result1, "filterResult")
      expect_equal(sum(b08 %in% filter1), 8291)
      
      e07.result1 = filter(e07,filter1)
      expect_equal(sum(e07 %in% filter1), 8514)
      
      f06.result1 = filter(f06,filter1)
      expect_equal(sum(f06 %in% filter1), 8765)
      
      
    })

