context("transformation of gate co-ordinates")

# PE inverse transformer with cofactor = 50
pe_sinh_trans <- function(x, cofactor = 50) {
  sinh(x) * cofactor
}

# APC inverse transformer with cofactor = 100
apc_sinh_trans <- function(x, cofactor = 100) {
  sinh(x) * cofactor
}

# extract inverse transformers to transformList
inv_trans <- transformList(
  c("PE-A", "APC-A"),
  list(
    "PE-A" = pe_sinh_trans,
    "APC-A" = apc_sinh_trans
  )
)

# NOTE: here we create gates on the transformed scale 
# PE-A range = c(0, 9.26) - asinh(x/50)
# APC-A range = c(0, 8.56) - asinh(x/100)

# rectangleGate
rg <- rectangleGate(
  "PE-A" = c(4, 6),
  "APC-A" = c(3.5, 7),
  filterId = "rect"
)

# polygonGate
pg <- as(rg, "polygonGate")

# ellipsoidGate
eg <- ellipsoidGate(
  .gate = matrix(
    c(
      6879, 
      3612, 
      3612, 
      5215
    ), 
    ncol=2,
    dimnames = list(
      c("PE-A", "APC-A"), 
      c("PE-A", "APC-A")
    )
  ),
  mean = c(
    "PE-A" = 4.5,
    "APC-A" = 4
  ),
  filterId = "ellipse"
)

# quadGate
qg <- quadGate(
  "PE-A" = 4.5,
  "APC-A" = 4,
  filterId = "quad"
)

# filters - all supported gate types
gates <- filters(
  list(
    "rect" = rg,
    "poly" = pg,
    "ellipse" = eg,
    "quad" = qg
  )
)

# test filters method to cover all gate type tests
test_that(
  "transform gate co-ordinates", {
    # transform to get linear gates
    gates_inv <- transform(
      gates,
      inv_trans
    )
    # rectangleGate
    expect_equal(
      rbind(
        gates_inv[["rect"]]@min,
        gates_inv[["rect"]]@max
      ),
      matrix(
        c(
          1364.5,
          1654.3,
          10085.7,
          54831.6
        ),
        nrow = 2,
        ncol = 2,
        byrow = TRUE,
        dimnames = list(
          NULL,
          c("PE-A", "APC-A")
        )
      ),
      tolerance = 0.1
    )
    
    # polygonGate
    expect_equal(
      gates_inv[["poly"]]@boundaries,
      matrix(
        c(
          1364.5,  1654.3,
          10085.7,  1654.3,
          10085.7, 54831.6,
          1364.5, 54831.6
        ),
        ncol = 2,
        nrow = 4,
        byrow = TRUE,
        dimnames = list(
          NULL,
          c("PE-A", "APC-A")
        )
      ), 
      tolerance = 0.1
    )
    
    # ellipsoidGate -> calls polygonGate method -> no need to check coords here
    expect_is(
      gates_inv[["ellipse"]],
      "polygonGate"
    )
    
    # quadGate
    expect_equal(
      gates_inv[["quad"]]@boundary,
      c(
        "PE-A" = 2250.2,
        "APC-A" = 2729
      ),
      tolerance = 0.1
    )
  }
)