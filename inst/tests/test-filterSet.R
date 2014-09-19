test_that("filterSet", {
    b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"),transformation="scale")
    fs <- new("filterSet")
    fs[["filter1"]] <- rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))
    
    expect_is(fs[["filter1"]], "rectangleGate")
    expect_is(fs["filter1"], "filterReference")
    
    expect_equal(sum(b08 %in% fs["filter1"]), 8291)
    
    b08.result1 <- filter(b08,fs["filter1"])
    expect_is(b08.result1, "filterResult")
    expect_equal(summary(b08.result1)@p, 0.8291)
    
    fs[[""]] <- norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")
    
    filter1 <- rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))
    filter2 <- norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")
    
    #these can't be evaluated within local session
#    expect_is(as("filter1","filter"), "filterReference")
#    expect_is(as(as.name("filter1"),"filter"), "filterReference")
    
#    expect_is(as(~ filter1 %subset% filter2,"filter"), "subsetFilter")
    
    fs[["Combined"]] <- ~ `Live Cells` %subset% filter1
    expect_is(as(fs,"list"), "list")
    
    expect_is(sort(fs,dependencies=TRUE), "character")
    f <- filter(b08,fs)
    expect_is(f, "manyFilterResult")
    expect_is(as.data.frame(f), "data.frame")
    expect_is(f@subSet[1:10,], "matrix")
    expect_equal(colSums(f@subSet), c(filter1 = 8291, `Live Cells` =  6717, "Combined" = 6496))

    expect_is(f[[1]], "filterResult")
    expect_is(f[[2]], "filterResult")
    expect_is(f[[3]], "filterResult")
    
    expect_equal(rownames(f@dependency)[rowSums(f@dependency)==0], "Combined")
    
    expect_is(split(b08,f,flowSet=TRUE), "flowSet")
    expect_is(split(b08,f,drop=TRUE,flowSet=TRUE), "flowSet")
    

})

