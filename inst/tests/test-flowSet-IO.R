context("flowSet IO...")

data(GvHD)
fs <- GvHD[1:2]
expectPD <- pData(fs)
expectPD[["Patient"]] <- as.integer(as.vector(expectPD[["Patient"]]))
expectPD[["Visit"]] <- as.integer(expectPD[["Visit"]])
expectPD[["name"]] <- I(paste0(expectPD[["name"]], ".fcs"))
rownames(expectPD) <- paste0(rownames(expectPD), ".fcs")

tmpdir <- tempfile()

test_that("write.flowSet", {
      suppressWarnings(write.flowSet(fs, tmpdir))
})

test_that("read.flowSet", {
      
      files <- list.files(tmpdir, pattern = "fcs")
      #no phenoData supplied
      fs1 <- read.flowSet(files, path = tmpdir)
      expect_equal(pData(fs1), expectPD[, "name", drop = F])
      
      
      anno <- list.files(tmpdir, pattern = "txt", full = T)
      pd <- Biobase::read.AnnotatedDataFrame(anno)
      pd[["name"]] <- I(paste0(pd[["name"]], ".fcs"))
      #with phenoData supplied
      suppressWarnings(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd))
      pData(fs1)["FCS_File"] <- NULL
      expect_equal(pData(fs1), expectPD)
      
      #pd without name
      pData(pd)[["name"]] <- NULL
      suppressWarnings(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd))
      pData(fs1)[["name"]] <- I(pData(fs1)[["name"]])
      pData(fs1)["FCS_File"] <- NULL
      expect_equal(pData(fs1), expectPD)
      
      #pd with wrong name
      pd <- Biobase::read.AnnotatedDataFrame(anno)
      pData(pd)[["name"]] <- paste0(pData(pd)[["name"]], "dummy")
      suppressWarnings(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd))
      expect_equal(pData(fs1), pData(pd))
      
      #create duplicated folder
      tmpdir1 <- tempfile()
      suppressWarnings(write.flowSet(fs, tmpdir1))
      #try to read both folders in
      files <- list.files(tmpdir, pattern = "fcs", full = T)
      files1 <- list.files(tmpdir1, pattern = "fcs", full = T)
      fs2 <- read.flowSet(c(files, files1))
      #check duplicates
      sn <- basename(files)
      sn1 <- paste0(sn, ".1")
      expect_equal(sampleNames(fs2), c(sn, sn1)) 
      
    })

test_that("phenoData<-", {
      pd <- expectPD
      
      # without name column
      pd[["name"]] <- NULL
      
      sampleNames(fs) <- paste0(sampleNames(fs), ".fcs") 
      pData(fs) <- pd
      #name column is added
      expect_equal(expectPD[["name"]], I(pData(fs)[["name"]]))
      
      # name to be different from rownames
      pd[["name"]] <- letters[1:2] 
      pData(fs) <- pd
      expect_equal(pData(fs)[["name"]], pd[["name"]])
      expect_equal(rownames(pData(fs)), rownames(expectPD))
      
      
    })