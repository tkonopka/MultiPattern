## tests for functions in MP_neighbors.R

cat("\ntest_neighbors.R ")


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


test_that("normalized distances error checks", {
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




###############################################################################
## Tests for finding neighbors

if (FALSE) {
test_that("finding nearest neighbors from S1", {
  ## item one is closest to item 2
  output_one= MPgetNeighborSet(ddistA, "S1", maxrank=1)
  output_full= MPgetNeighborSet(ddistA, "S1", maxrank=0)
  expected_one= c(S2=1)
  expected_full= c(S2=1, S3=2, S4=3)
  expect_equal(output_one, expected_one)
  expect_equal(output_full, expected_full)
})

test_that("finding nearest neighbors from S3", {
  ## item two is close to item 1 and 3
  output_two= MPgetNeighborSet(ddistA, "S3", maxrank=2)
  output_full= MPgetNeighborSet(ddistA, "S3", maxrank=0)
  expected_two= c(S4=1, S2=2)
  expected_full = c(S4=1, S2=2, S1=3)
  expect_equal(output_two, expected_two)
  expect_equal(output_full, expected_full)
})

test_that("finding nearest neighbors with ties", {
  ## check some hand-picked items
  output_full= MPgetNeighborSet(ddistB, "S1", maxrank=0)
  expected_full = c(S0=1, S2=2.5, S3=2.5, S4=4)
  expect_equal(output_full, expected_full)
  ## up to rank 2 should only give one hit (others are tied at higher ranks)
  output_two= MPgetNeighborSet(ddistB, "S1", maxrank=2)
  expected_two = c(S0=1)
  expect_equal(output_two, expected_two)
})
}

