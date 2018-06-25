## tests for functions in MP_neighbors.R

cat("\ntest_neighbors.R\n")


###############################################################################
## Prep data objects

## create small dataset just for testing
## this one has points along the diagonal (easy to see which points are closest)
dtestA = cbind(A=c(1,2,6,7), B=c(1,2,6,7))
rownames(dtestA) = paste0("S", 1:nrow(dtestA))
ddistA = dist(dtestA)


## another small dataset with equal distances between items
dtestB = cbind(A=c(0.9, 1, 2, 1, 2), B=c(0.9, 1, 1, 2, 2))
rownames(dtestB) = paste0("S", 0:(nrow(dtestB)-1))
ddistB = dist(dtestB)




###############################################################################
## Tests for normalizing distance matrices


test_that("normalization of distances performs error checks", {
  ## input must be a matrix-like object
  expect_error(MPrankNeighbors(letters))
  expect_error(MPrankNeighbors(NA))
  expect_error(MPrankNeighbors(NULL))
  ## non-square matrices are not god
  expect_error(MPrankNeighbors(dtestA))
  ## matrices with empty rows/cols are not good
  expect_error(MPrankNeighbors(dtestA[, c()]))
  expect_error(MPrankNeighbors(dtestA[c(), 1]))
})


test_that("normalization of dstiances gives error on non-numerics", {
  ## matrices must have numeric data (not letters of factors)
  badmatrix = data.frame(A=1:4, B=1:4, C=1:4, D=1:4)
  rownames(badmatrix) = LETTERS[1:4]
  expect_silent(MPrankNeighbors(badmatrix))
  badmatrix["A", "C"] = "zero"
  expect_error(MPrankNeighbors(badmatrix))
})


test_that("normalized distances are always matrices", {
  normA = MPrankNeighbors(ddistA)
  expect_equal(class(normA), "matrix")
})


test_that("normalized distances are 0-1", {
  unit_range = c(0, 1)
  normA = MPrankNeighbors(ddistA)
  expect_equal(range(normA), unit_range)
  normB = MPrankNeighbors(ddistB)
  expect_equal(range(normB), unit_range)
})


test_that("normalized distances have zero diagonals", {
  normA = MPrankNeighbors(ddistA)
  normB = MPrankNeighbors(ddistB)
  expect_equal(sum(diag(normA))+sum(diag(normB)), 0)  
})

