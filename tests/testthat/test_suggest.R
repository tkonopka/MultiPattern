## tests for creating series of MP conigurations with one command MPsuggestConfigs
## i.e. tests of automated choices of configurations

cat("\ntest_suggest.R ")


###############################################################################
## Prep data objects

num.samples = 10
snames = paste0("S", 1:num.samples)

## create dataset with real-valued fields
aa = MPdata4S[snames, c("D1", "D2")]
bb = MPdata6S[snames, c("D1", "D2")]
cc = MPdata9S[snames, c("D1", "D2")]
colnames(bb) = c("D3", "D4")
colnames(cc) = c("D5", "D6")
abc = cbind(aa, bb, cc)

## create dataset with binary fields
abc.bin = matrix(c(-1,9), ncol=4, nrow=num.samples)
colnames(abc.bin) = paste0("bin.", letters[1:4])
rownames(abc.bin) = snames

## create dataset withh binary fields with skew
abc.binskew = matrix(c(4, rep(8, num.samples-1)), ncol=4, nrow=num.samples)
colnames(abc.binskew) = paste0("skew.", letters[1:4])
rownames(abc.binskew) = snames

## create dataset with multi-valued fields
abc.multi = matrix(c(1,2,3,4), ncol=4, nrow=num.samples)
colnames(abc.multi) = paste0("multi.", letters[1:4])
rownames(abc.multi) = snames


## a large dataset (for triggering large object messages)
largeN = 2e4
largedata = cbind(A=1:largeN, B=1:largeN, C=1:largeN, D=1:largeN)
rownames(largedata) = paste0("S", 1:largeN)
mplarge = MPnew(rownames(largedata), data=list(large=largedata))



## extract configuration names
confNames = function(mp) {
  result = names(mp$configs)
  result = result[!grepl("^rnorm", result)]
  sort(result)
}




###############################################################################
## Tests for adding configs using automated suggestConfigs

## some labels for configurations that will appear in suggestConfigs
core.configs = c("euclidean", "canberra",
                 "clust.A2reg", "clust.A2alt",
                 "clust.C2reg", "clust.C2alt",
                 "clust.P2reg", "clust.P2alt",
                 "clust.S2reg", "clust.S2alt",
                 "clust.A3reg", "clust.A3alt",
                 "clust.C3reg", "clust.C3alt",
                 "clust.P3reg", "clust.P3alt",
                 "clust.S3reg", "clust.S3alt",
                 "subspaceR.1", "subspaceR.2",
                 "subspaceR.3", "subspaceR.4")
corepca.configs = c("PC1", "PC1.PC2")
coredbscan.configs = paste0("dbscan.", 1:4)
core.rnorm = c("rnorm.1", "rnorm.2")


test_that("suggest gives error when non-MP input", {
  expect_error(MPsuggestConfig(1:4, data="abc"))
})


test_that("suggest gives error when data is too small", {
  mpsmall = MPnew(snames[1:8], data=list(abc=abc[1:8,]))
  expect_error(MPsuggestConfig(mpsmall, data="abc"))
})


test_that("suggest displays message when in verbose mode", {
  mp = MPnew(snames, data=list(abc=abc))
  expect_silent(MPsuggestConfig(mp, data="abc", verbose=FALSE))
  mp = MPnew(snames, data=list(abc=abc))
  expect_message(MPsuggestConfig(mp, data="abc", verbose=TRUE))
  ## for large datasets expect time-warning
  expect_message(MPsuggestConfig(mplarge, data="large", verbose=TRUE), "wait")
})


test_that("suggest a series of configs (data with real-valued fields)", {
  mp = MPnew(snames, data=list(abc=abc))
  ## change settings (decreases default number of configuration based on randomness)
  mp$settings$num.random=2
  mp$settings$subspace.num.random=4
  mp$settings$num.PCs=2
  mp = MPsuggestConfig(mp, data="abc", verbose=FALSE)
  expected = c(paste0("AUTO:abc:real:", core.configs),
               paste0("AUTO:abc:real:", corepca.configs),
               paste0("AUTO:abc:real:", coredbscan.configs))
  expect_equal(confNames(mp), sort(expected))
})


test_that("suggest a series of configs (data with bin-valued fields)", {
  abc.binbin = cbind(abc.bin, abc.binskew)
  mp = MPnew(snames, data=list(abc.binbin=abc.binbin))
  ## change settings (decreases default number of configuration based on randomness)
  mp$settings$num.random=2
  mp$settings$subspace.num.random=4
  mp$settings$num.PCs=2
  mp = MPsuggestConfig(mp, data="abc.binbin", verbose=FALSE)
  expected = c(paste0("AUTO:abc.binbin:bin:", core.configs),
               paste0("AUTO:abc.binbin:binskew:", core.configs))
  expect_equal(confNames(mp), sort(expected))
})


test_that("suggest a series of configs (data with multi-valued fields)", {
  mp = MPnew(snames, data=list(abc.multi=abc.multi))
  ## change settings (decreases default number of configuration based on randomness)
  mp$settings$num.random=2
  mp$settings$subspace.num.random=4
  mp$settings$num.PCs=2
  mp = MPsuggestConfig(mp, data="abc.multi", verbose=FALSE)
  expected = c(paste0("AUTO:abc.multi:multi:", core.configs))
  expect_equal(confNames(mp), sort(expected))
})


test_that("suggest a series of configs (all data types at once)", {
  abc.all = cbind(abc, abc.bin, abc.binskew, abc.multi)
  mp = MPnew(snames, data=list(abc.all=abc.all))
  ## change settings (decreases default number of configuration based on randomness)
  mp$settings$num.random=2
  mp$settings$subspace.num.random=4
  mp$settings$num.PCs=2
  mp = MPsuggestConfig(mp, data="abc.all", verbose=FALSE)
  expected = c(paste0("AUTO:abc.all:real:", core.configs),
               paste0("AUTO:abc.all:real:", corepca.configs),
               paste0("AUTO:abc.all:real:", coredbscan.configs),
               paste0("AUTO:abc.all:bin:", core.configs),
               paste0("AUTO:abc.all:binskew:", core.configs),
               paste0("AUTO:abc.all:multi:", core.configs))
  expect_equal(confNames(mp), sort(expected))
})


test_that("suggest gives error when not MultiPattern object", {
  expect_error(MPsuggestConfig(1:4, data=c("abc"), verbose=FALSE))
})


test_that("suggest gives error with multiple datasets", {
  mp = MPnew(snames, data=list(abc.bin=abc.bin, abc.multi=abc.multi))
  expect_error(MPsuggestConfig(mp, data=c("abc.bin", "abc.multi"), verbose=FALSE))
})


test_that("suggest gives error when refers to inexistent dataset", {
  mp = MPnew(snames, data=list(abc.bin=abc.bin))
  expect_error(MPsuggestConfig(mp, data="abc.multi", verbose=FALSE))
})


test_that("suggest gives messages in verbose mode", {
  mp = MPnew(snames, data=list(abc.bin=abc.bin))
  expect_error(MPsuggestConfig(mp, data="abc.multi", verbose=FALSE))
})


test_that("suggest gives messages when data is single-values", {
  all0 = cbind(A=rep(0, 10), B=rep(0, 10))
  rownames(all0) = letters[1:10]
  mp = MPnew(rownames(all0), data=list(ZZ=all0))
  ## ask for automatic config suggestions
  expect_error(MPsuggestConfig(mp, data="ZZ", verbose=FALSE))
})


