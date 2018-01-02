## tests for functions in MP_distances.R

cat("\ntest_distances.R ")


###############################################################################
## Tests for imputing values in matrices


test_that("remove NAs from matrices or data frames", {
  ## create a test object
  xx = cbind(A=c(1,2,NA), B=c(NA, 4, 5))

  ## impute using mean from entire matrix
  output = MPimputeNAs(xx)
  expected = cbind(A=c(1,2,3), B=c(3, 4, 5))
  expect_equal(output, expected)
  
  ## impute by columns
  output2 = cbind(apply(xx, 2, MPimputeNAs))
  expected2 = cbind(A=c(1,2,1.5), B=c(4.5, 4, 5))
  expect_equal(output2, expected2)
})

test_that("remove infinite values from matrices or data frames", {
  ## create a test object
  xx = cbind(A=c(1,2,Inf), B=c(-Inf, 4, 5))

  ## impute using mean from entire matrix
  output = MPimputeNAs(xx)
  expected = cbind(A=c(1,2,5), B=c(1, 4, 5))
  expect_equal(output, expected)
})

test_that("remove infinite values from matrices or data frames (edge cases)", {
  ## create a test object with all bad values
  xx = cbind(A=c(Inf,Inf), B=c(-Inf, Inf))

  ## imputing here normalizes the values into 0,1 range
  expect_warning(MPimputeNAs(xx))
  output = suppressWarnings(MPimputeNAs(xx))
  expected = cbind(A=c(1, 1), B=c(0, 1))
  expect_equal(output, expected)
})



###############################################################################
## Tests for individual distance functions

test_that("classic euclidean/canberra/manhattan distances (easy cases)", {
  xx = cbind(A=c(1,2,3), B=c(2,3,4), C=c(3,4,5))
  rownames(xx) = letters[1:3]
  ## compute standard distances
  output_euclidean = dist.euclidean(xx)
  expected_euclidean = dist(xx, method="euclidean")
  output_manhattan = dist.manhattan(xx)
  expected_manhattan = dist(xx, method="manhattan")
  output_canberra = dist.canberra(xx)
  expected_canberra = dist(xx, method="canberra")
  ## specialized functions and stats:: functions should give same results
  expect_equal(output_euclidean, expected_euclidean, check.attributes=FALSE)
  expect_equal(output_manhattan, expected_manhattan, check.attributes=FALSE)
  expect_equal(output_canberra, expected_canberra, check.attributes=FALSE)
})


test_that("classic euclidean/canberra/manhattan distances (cases with NAs)", {
  xx = cbind(A=c(1,2,NA), B=c(2,Inf,4), C=c(3,4,5))
  rownames(xx) = letters[1:3]
  ## compute distances
  output_euclidean = dist.euclidean(xx)
  output_manhattan = dist.manhattan(xx)
  output_canberra = dist.canberra(xx)
  ## outputs should be all numeric without infinite values
  ## for euclidean
  expect_equal(sum(is.na(output_euclidean)), 0)
  expect_equal(sum(!is.finite(output_euclidean)), 0)
  ## for manhattan
  expect_equal(sum(is.na(output_manhattan)), 0)
  expect_equal(sum(!is.finite(output_manhattan)), 0)
  ## for canberra
  expect_equal(sum(is.na(output_canberra)), 0)
  expect_equal(sum(!is.finite(output_canberra)), 0)
})


test_that("additional dist functions spearman/pc1", {
  ## make a larger dataset so that spearman does not complain
  xx = cbind(A=c(1,2,3), B=c(2,3,4), C=c(3,4,5), D=11:13,
             E=c(6, 2, 8), F=c(1, 1, 0), G=c(1,4,1), H=c(0, 3, 8))
  rownames(xx) = letters[1:3]
  ## compute distances
  output_euclidean = dist.pc1(xx)
  output_spearman = dist.spearman(xx)
  output_log2euc = dist.log2euclidean(xx)
  ## tests are naive, should have non-zero values
  expect_gt(sum(output_euclidean), 0)
  expect_gt(sum(output_spearman), 0)
  expect_gt(sum(output_log2euc), 0)
})

test_that("tests for distances based on clustering", {
  
  xx = cbind(A=c(1,2,3,4,5), B=c(2,3,4,5,3), C=c(3,4,5,1,0), D=11:15)
  rownames(xx) = letters[1:5]
  
  ## compute distances
  output_clustreg = dist.clust(xx, clust.alt=FALSE)
  output_clustalt = dist.clust(xx, clust.alt=TRUE)
  output_pamreg = dist.clust(xx, clust.method="pam", clust.alt=FALSE)
  output_pamalt = dist.clust(xx, clust.method="pam", clust.alt=TRUE)
  
  ## tests are naive, should have non-zero values
  expect_gt(sum(output_clustreg), 0)
  expect_gt(sum(output_clustalt), 0)
  expect_gt(sum(output_pamreg), 0)
  expect_gt(sum(output_pamalt), 0)  

  ## the reg and alt should have different distances
  expect_gt(sum(abs(output_clustreg-output_clustalt)), 0)
  expect_gt(sum(abs(output_pamreg-output_pamalt)), 0)

  ## cannot have too many clusters
  expect_error(dist.clust(xx, clust.k=10))
})

test_that("tests for distances based on clustering without rownames", {
  xx = cbind(A=c(1,2,3,4,5), B=c(2,3,4,5,3), C=c(3,4,5,1,0), D=11:15)
  output = dist.clust(xx, clust.alt=FALSE)
  output_matrix = as.matrix(output)
  ## the dist function will create sample names using integers
  expect_equal(colnames(output_matrix), rownames(output_matrix))
  expect_equal(colnames(output_matrix), as.character(1:5))
})
