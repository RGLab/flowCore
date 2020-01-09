fr <- expectRes[["emptyFrame"]]

test_that("flowParamMatch ",{
      pd <- pData(parameters(fr))
      #complete word fixed match 
      #by channel name
      expect_equivalent(.flowParamMatch(pd, "<B710-A>", fix = TRUE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "<b710-a>", fix = TRUE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "B710-A", fix = TRUE, partial = FALSE), integer(0))
      #by marker name
      expect_equivalent(.flowParamMatch(pd, "CD4 PcpCy55", fix = TRUE, partial = FALSE), integer(0))#TOFIX
      expect_equivalent(.flowParamMatch(pd, "CD4", fix = TRUE, partial = FALSE), 5)#TOFIX
      expect_equivalent(.flowParamMatch(pd, "cd4", fix = TRUE, partial = FALSE), 5)#TOFIX
      expect_equivalent(.flowParamMatch(pd, "CD4 Pcp", fix = TRUE, partial = FALSE), integer(0))
      
      #complete word regexp match 
      #by channel name
      expect_equivalent(.flowParamMatch(pd, "<B710-A>", fix = FALSE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "<b710-a>", fix = FALSE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "B710-A", fix = FALSE, partial = FALSE), integer(0))
      #by marker name
      expect_equivalent(.flowParamMatch(pd, "CD4 PcpCy55", fix = FALSE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "CD4", fix = FALSE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "cd4", fix = FALSE, partial = FALSE), 5)
      expect_equivalent(.flowParamMatch(pd, "cd", fix = FALSE, partial = FALSE), integer(0))
      
      #partial regexp match 
      #by channel name
      expect_equivalent(.flowParamMatch(pd, "<B710-A>", fix = FALSE, partial = TRUE), 5)
      expect_equivalent(.flowParamMatch(pd, "<b710-a>", fix = FALSE, partial = TRUE), 5)
      expect_equivalent(.flowParamMatch(pd, "B710-A", fix = FALSE, partial = TRUE), 5)
      expect_equivalent(.flowParamMatch(pd, "b710", fix = FALSE, partial = TRUE), 5)
      expect_equivalent(.flowParamMatch(pd, "T710", fix = FALSE, partial = TRUE), integer(0))
      #by channel name
      expect_equivalent(.flowParamMatch(pd, "CD4 PcpCy55", fix = FALSE, partial = TRUE), 5)
      expect_equivalent(.flowParamMatch(pd, "CD4", fix = FALSE, partial = TRUE), c(5,11))
      expect_equivalent(.flowParamMatch(pd, "cd4", fix = FALSE, partial = TRUE), c(5,11))
      expect_equivalent(.flowParamMatch(pd, "cd45", fix = FALSE, partial = TRUE), 11)
      expect_equivalent(.flowParamMatch(pd, "cd46", fix = FALSE, partial = TRUE), integer(0))
      
    })

test_that("getChannelMarker ",{
      
      res <- data.frame(name = "<B710-A>", desc = "CD4 PcpCy55", row.names = "$P5", stringsAsFactors = FALSE)
      # match by partial/full marker name     
      expect_equivalent(getChannelMarker(fr, "CD4"), res)
      expect_equal(getChannelMarker(fr, "CD4 PcpCy55"), res)
      expect_equal(getChannelMarker(fr, "cd4"), res)
      expect_error(getChannelMarker(fr, "cd41"), "can't find cd41")
      expect_error(getChannelMarker(fr, "cd"), "multiple markers matched")
      
      # match by channel name
      expect_equal(getChannelMarker(fr, "<B710-A>"), res)
      suppressWarnings(expect_equal(getChannelMarker(fr, "B710-A"), res))
      expect_warning(getChannelMarker(fr, "B710-A"), "partially matched")
      
      res <- data.frame(name = "<G780-A>", desc = "CD45RA PECy7", row.names = "$P11", stringsAsFactors = FALSE)
      suppressWarnings(expect_equivalent(getChannelMarker(fr, "CD45"), res))
      expect_warning(getChannelMarker(fr, "CD45"), "partially matched")
      expect_equivalent(getChannelMarker(fr, "CD45RA"), res)
    
      #succeed when exact match found  
      fr <- flowFrame(matrix(1:100, ncol = 3, dimnames = list(NULL, c("cd18", "cd186", "cd18 apc"))))
      expect_error(getChannelMarker(fr, "cd1")[["name"]], "multiple")
      expect_equal(getChannelMarker(fr, "cd18")[["name"]], "cd18")
      expect_equal(getChannelMarker(fr, "cd186")[["name"]], "cd186")
      
      fr <- flowFrame(matrix(1:100, ncol = 3, dimnames = list(NULL, c("cd4", "cd186", "cd18 apc"))))
      expect_equal(getChannelMarker(fr, "cd18")[["name"]], "cd18 apc")
    })
