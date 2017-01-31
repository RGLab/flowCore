context("flowFrame methods...")
data("GvHD")
fr <- GvHD[[1]]
test_that("range", {
  expect_equal(range(fr)  ,data.frame("FSC-H" = c(0,1023)
                                      ,"SSC-H" = c(0,1023)
                                      ,"FL1-H" = c(1,10000)
                                      ,"FL2-H" = c(1,10000)
                                      ,"FL3-H" = c(1,10000)
                                      ,"FL2-A" = c(0,1023)
                                      ,"FL4-H" = c(1,10000)
                                      ,"Time" = c(0,1023)
                                      , row.names = c("min", "max")
                                      , check.names = FALSE
                                      )
                )
  
  expect_equal(range(fr, type = "data")  ,data.frame("FSC-H" = c(59,1023)
                                                    ,"SSC-H" = c(6,1023)
                                                    ,"FL1-H" = c(1,10000)
                                                    ,"FL2-H" = c(1.000,9221.666)
                                                    ,"FL3-H" = c(1.000,1131.784)
                                                    ,"FL2-A" = c(0,1023)
                                                    ,"FL4-H" = c(1,1162.77)
                                                    ,"Time" = c(1, 755)
                                                    , row.names = c("min", "max")
                                                    , check.names = FALSE
                                                  )
              , tolerance = 4e-7)
               
})
