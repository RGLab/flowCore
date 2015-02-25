context("flowSet IO...")

data(GvHD)
fs <- GvHD[1:2]
expectPD <- pData(fs)
expectPD[["Patient"]] <- as.integer(as.vector(expectPD[["Patient"]]))
expectPD[["Visit"]] <- as.integer(expectPD[["Visit"]])

tmpdir <- tempfile()

test_that("write.flowSet", {
      suppressWarnings(write.flowSet(fs, tmpdir))
})

test_that("read.flowSet", {
      
      files <- list.files(tmpdir, pattern = "fcs")
      #no phenoData supplied
      fs1 <- read.flowSet(files, path = tmpdir)
      sampleNames(fs1) <- gsub(".fcs", "", sampleNames(fs1)) 
      expect_equal(pData(fs1), expectPD[, "name", drop = F])
      
      
      anno <- list.files(tmpdir, pattern = "txt", full = T)
      pd <- Biobase::read.AnnotatedDataFrame(anno)
      
      #with phenoData supplied
      suppressWarnings(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd))
      sampleNames(fs1) <- gsub(".fcs", "", sampleNames(fs1))
      pData(fs1)["FCS_File"] <- NULL
      expect_equal(pData(fs1), expectPD)
      
      #pd without name
      pData(pd)[["name"]] <- NULL
      suppressWarnings(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd))
      sampleNames(fs1) <- gsub(".fcs", "", sampleNames(fs1))
      pData(fs1)["FCS_File"] <- NULL
      expect_equal(pData(fs1), expectPD)
      
      #pd with wrong name
      pd <- Biobase::read.AnnotatedDataFrame(anno)
      pData(pd)[["name"]] <- paste0(pData(pd)[["name"]], "dummy")
      suppressWarnings(expect_error(fs1 <- read.flowSet(files, path = tmpdir, phenoData = pd), "'name' column is not consistent with rownames "))
      
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
      pData(fs) <- pd
      #name column is added
      expect_equal(expectPD[["name"]], pData(fs)[["name"]])
      
      # wrong name
      pd[["name"]] <- letters[1:2] 
      expect_error(pData(fs) <- pd, "'name' column is not consistent with rownames ")
      
      
    })