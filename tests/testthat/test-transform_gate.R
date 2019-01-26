context("Test gate transform methods")

test_polygon <-  rbind(
  c(2336.005, 2000.0125),
  c(2144.009, 1520.0400),
  c(2224.007, 912.0965),
  c(2336.005, 752.0873),
  c(2576.003, 688.0765),
  c(3311.534, 1051.7277),
  c(3584.000, 1280.0648),
  c(3600.000, 2016.0120),
  c(3568.000, 2240.0067),
  c(3360.000, 2304.0057),
  c(3200.001, 2224.0070),
  c(3040.001, 2192.0076),
  c(2704.002, 2048.0111),
  c(2560.003, 2032.0115))
colnames(test_polygon) <-  c("x", "y")
test_polygonGate <- polygonGate(filterId = "TestPolyGate", .gate = test_polygon)

test_cov <- rbind(c(230000, 90000), c(90000, 360000))
colnames(test_cov) <- c("x", "y")
test_mean <- c(700, 1600)
test_ellipsoidGate <- ellipsoidGate(filterId = "TestEllipseGate", cov = test_cov, mean = test_mean)

test_minmax <- list("x"=c(300, 800), "y"=c(1300,1900))
test_rectangleGate <- rectangleGate(filterId="TestRectGate", test_minmax)

test_boundary <- c(743, 1894)
names(test_boundary) <-  c("x", "y")
test_quadGate <- quadGate(filterId= "TestQuadGate", .gate=test_boundary)

# Just some factors used for testing
test_scale_uniform <- 4
test_scale_split <- c(2,7)
tol <- 1e-10
dx1 <- 42
dy1 <- 105
dx2 <- -17
dy2 <- 74
new_center <- c(932, 1437)
new_center2 <- c(841, -35)
test_angle <- 117
test_angle2 <- 63

test_that("scale polygonGate", {
  centroid <- poly_centroid(test_polygonGate@boundaries)
  scaled_uniform <- scale_gate(test_polygonGate, test_scale_uniform)
  scaled_indiv <- scale_gate(test_polygonGate, test_scale_split)
  diff_original <- sweep(test_polygonGate@boundaries, 2, centroid)
  diff_uniform <- sweep(scaled_uniform@boundaries, 2, centroid)
  diff_indiv <- sweep(scaled_indiv@boundaries, 2, centroid)
  match1 <- matrix(rep(test_scale_uniform, 28), c(14,2))
  colnames(match1) <-  c("x", "y")
  match2 <- matrix(c(rep(test_scale_split[1],14), rep(test_scale_split[2],14)), c(14,2))
  colnames(match2) <-  c("x", "y")
  expect_equal(diff_uniform / diff_original, match1)
  expect_equal(diff_indiv / diff_original, match2)
})

test_that("shift polygonGate", {
  diff <- new_center - poly_centroid(test_polygonGate@boundaries)
  shifted <- shift_gate(test_polygonGate, dx1)
  expect_equal(test_polygonGate@boundaries[,1] + dx1, shifted@boundaries[,1])
  expect_equal(test_polygonGate@boundaries[,2], shifted@boundaries[,2])
  shifted <- shift_gate(test_polygonGate, dx2, dy1)
  expect_equal(test_polygonGate@boundaries[,1] + dx2, shifted@boundaries[,1])
  expect_equal(test_polygonGate@boundaries[,2] + dy1, shifted@boundaries[,2])
  shifted <- shift_gate(test_polygonGate, c(dx2, dy1))
  expect_equal(test_polygonGate@boundaries[,1] + dx2, shifted@boundaries[,1])
  expect_equal(test_polygonGate@boundaries[,2] + dy1, shifted@boundaries[,2])
  shifted <- shift_gate(test_polygonGate, center = new_center)
  expect_equal(sweep(test_polygonGate@boundaries, 2, -diff), shifted@boundaries)
})

test_that("rotate polygonGate", {
  centroid <- poly_centroid(test_polygonGate@boundaries)
  rotated <- rotate_gate(test_polygonGate, test_angle)
  diff_original <- sweep(test_polygonGate@boundaries, 2, centroid)
  diff_rotated <- sweep(rotated@boundaries, 2, centroid)
  # Any given vector from the centroid should be rotated same amount
  a <- diff_original[7,]
  b <- diff_rotated[7,]
  # This will just give the positive angle, hence the abs
  theta <- (acos(sum(a*b)/(sqrt(sum(a * a)) * sqrt(sum(b * b)))))*180/pi
  expect_equal(theta, abs(test_angle))
})

test_that("scale ellipsoidGate", {
  scaled_uniform <- scale_gate(test_ellipsoidGate, test_scale_uniform)
  scaled_indiv <- scale_gate(test_ellipsoidGate, test_scale_split)
  evs_original <- eigen(test_ellipsoidGate@cov, symmetric = TRUE)
  evs_uniform <- eigen(scaled_uniform@cov, symmetric = TRUE)
  evs_indiv <- eigen(scaled_indiv@cov, symmetric = TRUE)
  expect_equal(test_ellipsoidGate@mean, scaled_uniform@mean)
  expect_equal(evs_original$values * test_scale_uniform, evs_uniform$values)
  expect_true(all.equal(evs_original$vectors[,1], evs_uniform$vectors[,1], tolerance = tol) || 
                all.equal(evs_original$vectors[,1], -evs_uniform$vectors[,1], tolerance = tol))
  expect_true(all.equal(evs_original$vectors[,2], evs_uniform$vectors[,2], tolerance = tol) || 
                all.equal(evs_original$vectors[,2], -evs_uniform$vectors[,2], tolerance = tol))
  # Accounts for re-sort that occurs if major axis becomes minor
  new_evs <- evs_original$values * test_scale_split
  sort_order <- if(new_evs[1] < new_evs[2]) c(2,1) else c(1,2)
  expect_equal(test_ellipsoidGate@mean, scaled_indiv@mean)
  expect_equal(sort(evs_original$values * test_scale_split, decreasing = TRUE), evs_indiv$values)
  # Allow for -1 scaling of eigenvectors
  expect_true(isTRUE(all.equal(evs_original$vectors[,1], evs_indiv$vectors[,sort_order[1]], tolerance = tol)) || 
                isTRUE(all.equal(evs_original$vectors[,1], -evs_indiv$vectors[,sort_order[1]], tolerance = tol)))
  expect_true(isTRUE(all.equal(evs_original$vectors[,2], evs_indiv$vectors[,sort_order[2]], tolerance = tol)) || 
                isTRUE(all.equal(evs_original$vectors[,2], -evs_indiv$vectors[,sort_order[2]], tolerance = tol)))
})


test_that("shift ellipsoidGate", {
  diff <- new_center - test_ellipsoidGate@mean
  shifted <- shift_gate(test_ellipsoidGate, dx1)
  expect_equal(test_ellipsoidGate@mean[1] + dx1, shifted@mean[1])
  expect_equal(test_ellipsoidGate@mean[2], shifted@mean[2])
  expect_equal(test_ellipsoidGate@cov, shifted@cov)
  shifted <- shift_gate(test_ellipsoidGate, dx2, dy1)
  expect_equal(test_ellipsoidGate@mean[1] + dx2, shifted@mean[1])
  expect_equal(test_ellipsoidGate@mean[2] + dy1, shifted@mean[2])
  expect_equal(test_ellipsoidGate@cov, shifted@cov)
  shifted <- shift_gate(test_ellipsoidGate, c(dx2, dy1))
  expect_equal(test_ellipsoidGate@mean[1] + dx2, shifted@mean[1])
  expect_equal(test_ellipsoidGate@mean[2] + dy1, shifted@mean[2])
  expect_equal(test_ellipsoidGate@cov, shifted@cov)
  shifted <- shift_gate(test_ellipsoidGate, center = new_center)
  expect_equal(test_ellipsoidGate@mean + diff, shifted@mean)
  expect_equal(test_ellipsoidGate@cov, shifted@cov)
})

test_that("rotate ellipsoidGate", {
  rad <- test_angle*(pi/180)
  rot <- rbind(c(cos(rad), -sin(rad)), c(sin(rad), cos(rad)))
  rotated_cov <- rot%*%(test_ellipsoidGate@cov)%*%t(rot)
  colnames(rotated_cov) <- colnames(test_ellipsoidGate@cov)
  rotated <- rotate_gate(test_ellipsoidGate, test_angle)
  expect_equal(rotated_cov, rotated@cov)
  a <- eigen(test_ellipsoidGate@cov, symmetric = TRUE)$vectors[,1]
  b <- eigen(rotated@cov, symmetric = TRUE)$vectors[,1]
  theta <- (acos(sum(a*b)/(sqrt(sum(a * a)) * sqrt(sum(b * b)))))*180/pi
  expect_equal(theta, abs(test_angle))
})

test_that("scale rectangleGate", {
  scaled_uniform <- scale_gate(test_rectangleGate, test_scale_uniform)
  scaled_indiv <- scale_gate(test_rectangleGate, test_scale_split)
  center_original <- (test_rectangleGate@max + test_rectangleGate@min) / 2
  center_uniform <- (scaled_uniform@max + scaled_uniform@min) / 2
  center_indiv <- (scaled_indiv@max + scaled_indiv@min) / 2
  span_original <- test_rectangleGate@max - test_rectangleGate@min
  span_uniform <- scaled_uniform@max - scaled_uniform@min
  span_indiv <- scaled_indiv@max - scaled_indiv@min
  expect_equal(center_original, center_uniform)
  expect_equal(center_original, center_indiv)
  expect_equal(span_original*test_scale_uniform, span_uniform)
  expect_equal(span_original*test_scale_split, span_indiv)
})

test_that("shift rectangleGate", {
  shifted <- shift_gate(test_rectangleGate, dx1)
  expect_equal(test_rectangleGate@min[1] + dx1, shifted@min[1])
  expect_equal(test_rectangleGate@max[1] + dx1, shifted@max[1])
  expect_equal(test_rectangleGate@min[2], shifted@min[2])
  expect_equal(test_rectangleGate@max[2], shifted@max[2])
  shifted <- shift_gate(test_rectangleGate, dx2, dy1)
  expect_equal(test_rectangleGate@min[1] + dx2, shifted@min[1])
  expect_equal(test_rectangleGate@max[1] + dx2, shifted@max[1])
  expect_equal(test_rectangleGate@min[2] + dy1, shifted@min[2])
  expect_equal(test_rectangleGate@max[2] + dy1, shifted@max[2])
  shifted <- shift_gate(test_rectangleGate, c(dx2, dy1))
  expect_equal(test_rectangleGate@min[1] + dx2, shifted@min[1])
  expect_equal(test_rectangleGate@max[1] + dx2, shifted@max[1])
  expect_equal(test_rectangleGate@min[2] + dy1, shifted@min[2])
  expect_equal(test_rectangleGate@max[2] + dy1, shifted@max[2])
  shifted <- shift_gate(test_rectangleGate, center = new_center)
  center_shifted <- (shifted@max + shifted@min) / 2
  named_new_center <- new_center
  names(named_new_center) <- names(test_rectangleGate@min)
  expect_equal(named_new_center, center_shifted)
})

test_that("rotate rectangleGate", {
  expect_error(rotate_gate(test_rectangleGate, test_angle))
})

test_that("scale quadGate", {
  scaled_uniform <- scale_gate(test_quadGate, test_scale_uniform)
  scaled_indiv <- scale_gate(test_quadGate, test_scale_split)
  expect_equal(test_quadGate@boundary*test_scale_uniform, scaled_uniform@boundary)
  expect_equal(test_quadGate@boundary*test_scale_split, scaled_indiv@boundary)
})

test_that("shift quadGate", {
  shifted <- shift_gate(test_quadGate, dx1)
  expect_equal(test_quadGate@boundary[1] + dx1, shifted@boundary[1])
  expect_equal(test_quadGate@boundary[2], shifted@boundary[2])
  shifted <- shift_gate(test_quadGate, dx2, dy1)
  expect_equal(test_quadGate@boundary[1] + dx2, shifted@boundary[1])
  expect_equal(test_quadGate@boundary[2] + dy1, shifted@boundary[2])
  shifted <- shift_gate(test_quadGate, c(dx2, dy1))
  expect_equal(test_quadGate@boundary[1] + dx2, shifted@boundary[1])
  expect_equal(test_quadGate@boundary[2] + dy1, shifted@boundary[2])
  shifted <- shift_gate(test_quadGate, center = new_center)
  named_new_center <- new_center
  names(named_new_center) <- names(test_quadGate@boundary)
  expect_equal(named_new_center, shifted@boundary)
})

test_that("rotate quadGate", {
  expect_error(rotate_gate(test_quadGate, test_angle))
})

test_that("transform polygonGate", {
  # Mix dilation, rotation, and translation. Just make sure it is appropriately applying all 3 transformations in that order.
  # Logic of each of the transformations is mostly covered by earlier individual tests.
  transformed <- transform_gate(test_polygonGate, scale = test_scale_split, dx = c(dx1, dy1), deg = test_angle)
  composition <- shift_gate(rotate_gate(scale_gate(test_polygonGate, scale = test_scale_split), deg = test_angle), dx = c(dx1, dy1))
  expect_equal(transformed@boundaries, composition@boundaries)
  
  # Direct parameter modification followed by another composite transformation. Make sure it behaves in that order.
  new_polygon <- sweep(test_polygon, 2, -c(dx1, dy1))
  new_id <- "AnotherPolyGate"
  transformed <- transform_gate(test_polygonGate, scale = test_scale_split, boundaries = new_polygon, 
                                dx = c(dx2, dy2), filterId = new_id, deg = test_angle)
  pre_shifted <- shift_gate(test_polygonGate, c(dx1, dy1))
  composition <- shift_gate(rotate_gate(scale_gate(pre_shifted, scale = test_scale_split), deg = test_angle), dx = c(dx2, dy2))
  expect_equal(transformed@boundaries, composition@boundaries)
  expect_equal(transformed@filterId, new_id)
})


test_that("transform ellipsoidGate", {
  # Mix dilation, rotation, and translation
  transformed <- transform_gate(test_ellipsoidGate, scale = test_scale_split, dx = c(dx1, dy1), deg = test_angle)
  composition <- shift_gate(rotate_gate(scale_gate(test_ellipsoidGate, scale = test_scale_split), deg = test_angle), dx = c(dx1, dy1))
  expect_equal(transformed@mean, composition@mean)
  expect_equal(transformed@cov, composition@cov)
  
  # Direct parameter modification followed by another composite transformation. Make sure it behaves in that order.
  rad <- test_angle2*(pi/180)
  rot <- rbind(c(cos(rad), -sin(rad)), c(sin(rad), cos(rad)))
  new_cov <- rot%*%(test_ellipsoidGate@cov)%*%t(rot)
  new_id <- "AnotherEllipsoidGate"
  new_mean <- test_ellipsoidGate@mean + c(dx1, dy1)
  transformed <- transform_gate(test_ellipsoidGate, scale = test_scale_split, cov = new_cov, 
                                dx = c(dx2, dy2), mean = new_mean, filterId = new_id, deg = test_angle)
  pre_processed <- shift_gate(rotate_gate(test_ellipsoidGate, test_angle2), c(dx1, dy1))
  composition <- shift_gate(rotate_gate(scale_gate(pre_processed, scale = test_scale_split), deg = test_angle), dx = c(dx2, dy2))
  expect_equal(transformed@mean, composition@mean)
  expect_equal(transformed@cov, composition@cov)
  expect_equal(transformed@filterId, new_id)
})


test_that("transform rectangleGate", {
  # Rotation is invalid, so error
  expect_error(transform_gate(test_rectangleGate, scale = test_scale_split, dx = c(dx1, dy1), deg = test_angle))
  # Just dilation and translation
  transformed <- transform_gate(test_rectangleGate, scale = test_scale_split, dx = c(dx1, dy1))
  composition <- shift_gate(scale_gate(test_rectangleGate, scale = test_scale_split), dx = c(dx1, dy1))
  expect_equal(transformed@min, composition@min)
  expect_equal(transformed@max, composition@max)
  
  # Direct parameter modification followed by another composite transformation. Make sure it behaves in that order.
  new_id <- "AnotherRectangleGate"
  new_min <- test_rectangleGate@min + c(dx1, dy1)
  new_max <- test_rectangleGate@max + c(dx1, dy1)
  transformed <- transform_gate(test_rectangleGate, scale = test_scale_split, min = new_min, 
                                dx = c(dx2, dy2), max = new_max, filterId = new_id)
  pre_processed <- shift_gate(test_rectangleGate, c(dx1, dy1))
  composition <- shift_gate(scale_gate(pre_processed, scale = test_scale_split), dx = c(dx2, dy2))
  expect_equal(transformed@min, composition@min)
  expect_equal(transformed@max, composition@max)
  expect_equal(transformed@filterId, new_id)
})

test_that("transform quadGate", {
  # Rotation is invalid, so error
  expect_error(transform_gate(test_quadGate, scale = test_scale_split, dx = c(dx1, dy1), deg = test_angle))
  # Just dilation and translation
  transformed <- transform_gate(test_quadGate, scale = test_scale_split, dx = c(dx1, dy1))
  composition <- shift_gate(scale_gate(test_quadGate, scale = test_scale_split), dx = c(dx1, dy1))
  expect_equal(transformed@boundary, composition@boundary)
  
  # Direct parameter modification followed by another composite transformation. Make sure it behaves in that order.
  new_id <- "AnotherQuadGate"
  new_boundary <- test_quadGate@boundary + c(dx1, dy1)
  transformed <- transform_gate(test_quadGate, scale = test_scale_split, boundary = new_boundary, 
                                dx = c(dx2, dy2), filterId = new_id)
  pre_processed <- shift_gate(test_quadGate, c(dx1, dy1))
  composition <- shift_gate(scale_gate(pre_processed, scale = test_scale_split), dx = c(dx2, dy2))
  expect_equal(transformed@boundary, composition@boundary)
  expect_equal(transformed@filterId, new_id)
})


