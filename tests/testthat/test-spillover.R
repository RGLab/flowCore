context("Test spillover methods")

control_path <- system.file("extdata", "compdata", "data",
                            package="flowCore")
matchfile_path <- system.file("extdata", "compdata", "comp_match",
                              package="flowCore")
controls <- lapply(dir(control_path, full.names=TRUE), 
                   read.FCS)
filenames <- sapply(controls, keyword, "$FIL")
names(controls) <- c("UNSTAINED", "FL1-H", "FL2-H", "FL4-H", "FL3-H")
controls <- as(controls, "flowSet")
# Scramble the columns (particularly fsc, ssc) to make sure
# the methods still handle it appropriately
controls <- controls[,c(3,1,4,5,2,6,7)]


test_that("spillover: Match columns using regexpr", {
  ref_file <- system.file("extdata", "compdata", "compref1", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  comp <- spillover(controls, unstained = "UNSTAINED", fsc="FSC-H", ssc="SSC-H", patt = "-H", stain_match="regexpr")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})

test_that("spillover: Match columns using ordered", {
  ref_file <- system.file("extdata", "compdata", "compref2", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  # The spillover matrix here should (appropriately) be incorrect bc the rows are
  # not in the order of the columns.
  comp <- spillover(controls, unstained = "UNSTAINED", fsc="FSC-H", ssc="SSC-H", patt = "-H", stain_match="ordered")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})

test_that("spillover: Match columns using intensity", {
  ref_file <- system.file("extdata", "compdata", "compref2", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  comp <- spillover(controls, unstained = "UNSTAINED", fsc="FSC-H", ssc="SSC-H", patt = "-H", stain_match="intensity")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})


test_that("spillover_match: Using path to dir with files", {
  ref_file <- system.file("extdata", "compdata", "compref3", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  matched <- spillover_match(path=control_path, fsc="FSC-H", ssc="SSC-H", matchfile = matchfile_path)
  comp <- spillover(matched, prematched = TRUE, fsc="FSC-H", ssc="SSC-H", patt = "-H",
                    useNormFilt = TRUE, method = "mode")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})

test_that("spillover_match: Using preconstructed flowSet with filenames as sample names", {
  sampleNames(controls) <- filenames
  ref_file <- system.file("extdata", "compdata", "compref3", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  matched <- spillover_match(controls, fsc="FSC-H", ssc="SSC-H", matchfile = matchfile_path)
  comp <- spillover(matched, prematched = TRUE, fsc="FSC-H", ssc="SSC-H", patt = "-H",
                    useNormFilt = TRUE, method = "mode")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})



