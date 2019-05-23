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
      
      
          
      filter2 = flowStats::norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")
      b08.result2 = filter(b08,filter2 %subset% b08.result1)
      expect_equal(filterDetails(b08.result2), expectRes[["filter.norm2Filter"]])
      expect_equal(sum(b08 %in% (filter2 %subset% b08.result1)), 6496)
      
      e07.result2 = filter(e07,filter2 %subset% e07.result1)
      expect_equal(sum(e07 %in% (filter2 %subset% e07.result1)), 6416)
      
      f06.result2 = filter(f06,filter2 %subset% f06.result1)
      expect_equal(sum(f06 %in% (filter2 %subset% f06.result1)), 6959, tol = 2e-3)
      
      ## the third-fifth gates get the positive cells for the marker in FL1-H
      ## this is a really interesting example because it illustrates that there
      ## are two subpopulations. Naturally we would like to automatically find them
      ## In this case we want to now what percent the positive population in FL1-H is of the
      ## total population

      filter3 = rectangleGate("FL1-H"=c(.4,Inf),filterId="FL1-H+")
      b08.result3 = filter(b08,filter3 %subset% b08.result2)
      expect_equal(summary(b08.result3),expectRes[["filter.b08.result3"]]) 
      # An example of doing manipulations with filterResult objects 
      # since they are now also filter objects. The %subset% argument
      # has a special summary meaning because it does it's calculations
      # relative to the RHS argument, whereas filterResults are ALWAYS
      # wrt to the full dataset. So, the value below should differ from 
      # the value above.
      expect_equal(summary(b08.result3 %subset% b08.result2), expectRes[["filter.b08.result3.subset"]])


      filter4 = flowStats::norm2Filter("FSC-H","FL1-H",scale.factor=2,filterId="FL1-H+")
      b08.result4 = filter(b08,filter4 %subset% b08.result2)
      expect_equal(summary(b08.result4),expectRes[["filter.b08.result3.subset2"]]) 

      expect_equal(summary(b08.result4 %subset% b08.result2)@p, 0.4802)
      
      #We now have NOT gates so we can proceed!
      b08.result5 = filter(b08,filter4 %subset% (b08.result2 & !b08.result4))
      expect_equal(summary(b08.result5)@true, 1407)

      expect_equal(summary(b08.result4 %subset% (b08.result4 | b08.result5))$p, 0.4802)
      
    })

