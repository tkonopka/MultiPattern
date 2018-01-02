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
MPeasyConfig(mptest, 
             type=list(A="hclust", B=c("euclidean", "manhattan")),
             )





###############################################################################
## Tests computing distances


test_that("compute all distances in a MP analysis", {
  result = MPgetDistances(mptest, verbose=FALSE)
  expect_equal(sort(names(result)), sort(names(mptest$configs)))
})


test_that("compute some distances in an MP analysis", {
  myconf = names(mptest$configs)
  myconf = myconf[grep("clust.A", myconf)]
  result = MPgetDistances(mptest, configs=myconf, verbose=FALSE)
  expect_equal(sort(names(result)), sort(myconf))
})



###############################################################################
## Tests computing meta-distances


test_that("compute averge meta-dstances in an MP analysis (using subsampling)", {
  result1 = MPgetAverageMetaDistance(mptest, verbose=FALSE)
  expect_equal(colnames(result1), names(mptest$configs))
  expect_equal(rownames(result1), names(mptest$configs))
  ## should not give exactly same results on multiple repeats
  result2 = MPgetAverageMetaDistance(mptest, verbose=FALSE)
  expect_gt(sum(abs(result2-result1)), 0)
})


test_that("compute meta-distances using all configs/samples", {
  mpsims = MPgetDistances(mptest, verbose=FALSE)
  result.1 = MPgetMetaDistance(mpsims)
  result.2 = MPgetMetaDistance(mpsims)
  expect_equal(result.1, result.2)
  result.a2 = MPgetMetaDistance(mpsims, alpha=2)
  result.b1 = MPgetMetaDistance(mpsims, beta=1)
  expect_false(identical(result.1, result.a2))
  expect_false(identical(result.1, result.b1))
})




###############################################################################
## Tests for getting representatives


test_that("compute meta-distances using all configs/samples", {
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


