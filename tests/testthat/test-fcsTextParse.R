context("FCS header parsing...")

test_that("fcsTextParse", {
      
      txt <- "/k1/v1/k2/v2/"
      res <- c(k1 = "v1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      
      #single characer value
      txt <- "/k1/1/k2/v2/"
      res <- c(k1 = "1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      txt <- "/k1/1/k/v/"
      res <- c(k1 = "1", k = "v")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      
      #without ending delimiter
      txt <- "/k1/v1/k2/v2"
      res <- c(k1 = "v1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      txt <- "/k1/v1/k2/v2/k3"
      expect_output(thisRes <- fcsTextParse(txt, empty = TRUE), "dropped")
      expect_equal(thisRes, res)
      
      #unpaired kw
      txt <- "/k1/v1/k2/v2/k3/"
      expect_output(fcsTextParse(txt, empty = TRUE), "uneven number")
            
      #with double delimiters
      txt <- "/k1/v//1/k2/v2/"
      expect_error(fcsTextParse(txt, empty = TRUE), "Empty keyword name")
      res <- c(k1 = "v/1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = FALSE), res)
      
      
      #empty value
      txt <- "/k1//k2/v2/"
      res <- c(k1 = "", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      #It will treat empty value as the double delimiter case and there is no way to detect this kind of logic error
      res <- c(`k1/k2` =  'v2')
      expect_equal(fcsTextParse(txt, empty = FALSE), res)  

      #special delimiter: \
      txt <- "\\k1\\v1\\k2\\v2\\"
      res <- c(k1 = "v1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      #with empty val
      txt <- "\\k1\\\\k2\\v2\\"
      res <- c(k1 = "", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), c(`k1\\k2` = "v2"))
      
      #special delimiter: |
      txt <- "|k1|v1|k2|v2|"
      res <- c(k1 = "v1", k2 = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      
      #special character: $
      txt <- "/$k1/v1/$k2/v2/"
      res <- c(`$k1` = "v1", `$k2` = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
      
      txt <- "\\k1\\1\\k2\\v2\\"
      res <- c(`k1` = "1", `k2` = "v2")
      expect_equal(fcsTextParse(txt, empty = TRUE), res)
      expect_equal(fcsTextParse(txt, empty = F), res)
    })
