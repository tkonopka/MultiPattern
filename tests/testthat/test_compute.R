## tests for computing based on a MultiPattern configuration

cat("\ntest_compute.R ")


###############################################################################
## Prep data objects

## samples
num.samples = 10
snames = paste0("S", 1:num.samples)


## object used while testing adding configurations
prep.list = list(one="D1", both=c("D1", "D2"))
dist.list = list(euc=dist.euclidean, man=dist.manhattan)


## object with multiple datasets and configurations
mptest = MPnew(snames, data=list(A=MPdata4S[, 1:2], B=MPdata6S[,1:2]))
mptest$settings$subsample.N = 10
mptest = MPeasyConfig(mptest, 
                      type=list(A="hclust", B=c("euclidean", "manhattan")),
                      )
mptest = MPremove(mptest, config=grep("3", names(mptest$configs), value=T))

## a large dataset (for triggering large object messages)
largeN = 1e4
largedata = cbind(A=1:largeN, B=1:largeN, C=1:largeN, D=1:largeN)
rownames(largedata) = paste0("S", 1:largeN)
mplarge = MPnew(rownames(largedata), data=list(large=largedata))
mplarge = MPeasyConfig(mplarge, type=c("euclidean", "manhattan"))




###############################################################################
## Tests computing distances


test_that("compute distances gives error with non-MP object", {
  expect_error(MPgetDistances(1:10))
})


test_that("compute all distances in a MP analysis", {
  result = MPgetDistances(mptest, verbose=FALSE)
  expect_equal(sort(names(result)), sort(names(mptest$configs)))
})


test_that("compute some distances in an MP analysis", {
  myconf = names(mptest$configs)
  myconf = myconf[grep("clust.A", myconf)]
  result = MPgetDistances(mptest, configs=myconf)
  expect_equal(sort(names(result)), sort(myconf))
})


test_that("compute distances gives error when configs don't match object", {
  myconf = c(names(mptest$configs)[1], "bad_name")
  expect_error(MPgetDistances(mptest, configs=myconf))
})


test_that("MPgetDistances is silent for small datasets", {
  myconf = names(mptest$configs)[1]
  expect_silent(MPgetDistances(mptest, configs=myconf, verbose=FALSE))
  ## this should be silent even when verbose is on because dataset is small
  expect_silent(MPgetDistances(mptest, configs=myconf, verbose=TRUE))
})


test_that("MPgetDistances prints time warnings for large datasets", {
  myconf = names(mplarge$configs)[1]
  expect_silent(MPgetDistances(mplarge, configs=myconf, verbose=FALSE))
  expect_message(MPgetDistances(mplarge, configs=myconf, verbose=TRUE))
})




###############################################################################
## Tests computing average meta-distances

test_that("compute avg meta-distances gives error with non-MP object", {
  expect_error(MPgetAverageMetaDistance(1:10))
})


test_that("compute avg meta-distances gives messages", {
  expect_silent(MPgetAverageMetaDistance(mptest, verbose=FALSE))
  expect_message(MPgetAverageMetaDistance(mptest, verbose=TRUE))
})


test_that("compute averge meta-dstances in an MP analysis (using subsampling)", {
  result1 = MPgetAverageMetaDistance(mptest, verbose=FALSE)
  expect_equal(colnames(result1), names(mptest$configs))
  expect_equal(rownames(result1), names(mptest$configs))
  ## should not give exactly same results on multiple repeats
  result2 = MPgetAverageMetaDistance(mptest, verbose=FALSE)
  expect_gt(sum(abs(result2-result1)), 0)
})


test_that("compute averge meta-dstances gives error when subsampling mis-specified", {
  expect_error(MPgetAverageMetaDistance(mptest, verbose=FALSE, subsample.N=1))
  expect_error(MPgetAverageMetaDistance(mptest, verbose=FALSE, subsample.N=1.5))
})


test_that("compute averge meta-dstances uses fractional or integer subsample size", {
  mptest$settings$subsample.R = 4
  set.seed(1234)
  result1 = MPgetAverageMetaDistance(mptest, subsample.N=0.8, verbose=FALSE)
  set.seed(1234)
  result2 = MPgetAverageMetaDistance(mptest, subsample.N=0.8*length(mptest$items),
                                     verbose=FALSE)
  expect_equal(result1, result2)
})




###############################################################################
## Tests computing complete meta-distances (no sub-sampling)

test_that("compute meta-distances gives error with wrong input class", {
  expect_error(MPgetMetaDistances(mptest))
})


test_that("compute meta-distances using all configs/samples", {
  mpsims = MPgetDistances(mptest, verbose=FALSE)
  result.1 = MPgetMetaDistance(mpsims)
  result.2 = MPgetMetaDistance(mpsims)
  expect_equal(result.1, result.2)
  result.a2 = MPgetMetaDistance(mpsims, alpha=2)
  result.b1 = MPgetMetaDistance(mpsims, beta=1)
  result.aneg = MPgetMetaDistance(mpsims, alpha=-0.5)
  expect_false(identical(result.1, result.a2))
  expect_false(identical(result.1, result.b1))
  expect_false(identical(result.1, result.aneg))
})




###############################################################################
## Tests for getting representatives

test_that("find representatives gives error with unusual inputs", {
  ## create some good and bad maps
  someN = 40
  somenames = paste0("S", 1:someN)
  badmap1 = matrix("a", nrow=someN, ncol=someN)
  rownames(badmap1) = colnames(badmap1) = somenames
  rownames(badmap1) = NULL
  goodmap = matrix(0, nrow=someN, ncol=someN)
  rownames(goodmap) = colnames(goodmap) = somenames
  ## test if getRepresentatives throws errors
  expect_error(MPgetRepresentatives(1:4))
  expect_error(MPgetRepresentatives(badmap1, k=4))
  expect_silent(MPgetRepresentatives(goodmap, k=4))
})


test_that("find representatives among a map", {
  mpmeta = MPgetAverageMetaDistance(mptest, verbose=FALSE)
  mapdist = dist(MPgetMap(mpmeta))
  ## with k<1, should select a fraction of available configurations 
  mpreps.c = MPgetRepresentatives(mapdist, method="complete", k=0.5)
  expect_equal(length(mpreps.c), length(mptest$configs)/2)
  ## with k>1, should select explicit number of configurations
  mpreps.c = MPgetRepresentatives(mapdist, method="complete", k=3)
  expect_equal(length(mpreps.c), 3)
  ## methods of selection should give different results (most often...)
  mpreps.e = MPgetRepresentatives(mapdist, method="extreme", k=3)
  expect_false(identical(mpreps.c, mpreps.e))
})


