controls <- lapply(dir(system.file("extdata", "compdata", "data",
                                   package="flowCore"), full.names=TRUE), 
                   read.FCS,
                   column.pattern="-H")

test_that("spillover: Match columns using regexpr", {
  ref_file <- system.file("extdata", "compdata", "compref", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  names(controls) <- c("UNSTAINED", "FL1-H", "FL2-H", "FL4-H", "FL3-H")
  controls <- as(controls, "flowSet")
  comp <- spillover(controls, unstained = "UNSTAINED", fsc="FSC-H", ssc="SSC-H", stain_match="regexpr")
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})
