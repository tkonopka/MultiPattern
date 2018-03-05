## tests for functions in MP_plot.R
## (These tests do not go very deep into the proper output of the plots.
## The tests are satisfied when the calls to plot do not yield errors 
## when a proper input is given)

cat("\ntest_plot.R ")


###############################################################################
## Data objects


## analyse a a dataset with several groups, but small number of samples
snames = rownames(MPdata4S)
mptest = MPnew(snames, data=list(A=MPdata4S[, 1:2]))
mptest$subsample.R = 4
mptest = MPeasyConfig(mptest, type=c("hclust", "euclidean"))
## to speed up examples, remove some configurations
mptest = MPremove(mptest,
                  config=c("A:clust.S2reg", "A:clust.S2alt", "A:clust.S3reg", "A:clust.S3alt"))
mpsims = MPgetDistances(mptest, verbose=FALSE)
mpmeta = MPgetAverageMetaDistance(mptest, verbose=FALSE)
mpmap = MPgetMap(mpmeta)

## for tests with point coloring
mapcolors = list()
mapcolors$reg = grep("reg", rownames(mpmap), value=T)
mapcolors$alt = grep("alt", rownames(mpmap), value=T)


## name of test file
testfile = "testplot.pdf"




###############################################################################
## tests for plot of cmdscale-like layout, e.g. meta-map

test_that("plot based on wrong input should give error", {
  bad.input = MPnew(letters)
  expect_error(MPplotmap(bad.input))
})

test_that("plot of one map should be silent", {
  good.input = cbind(A=1:4, B=1:4)
  pdf(file=testfile)
  expect_silent(MPplotmap(good.input, label=FALSE))
  expect_silent(MPplotmap(good.input, label=TRUE))
  dev.off()
  unlink(testfile)
})

test_that("plot of multiple maps in one file", {
  good.input = list(one=cbind(A=1:4, B=1:4), two=cbind(A=1:5, B=1:5))
  pdf(file=testfile)
  expect_silent(MPplotmap(good.input))
  expect_silent(MPplotmap(good.input, legend.separate=TRUE))
  dev.off()
  unlink(testfile)
})

test_that("plot of maps with legend on the side", {
  pdf(file=testfile)
  expect_silent(MPplotmap(mpmap, legend.separate=TRUE))
  dev.off()
  unlink(testfile)
})

test_that("plot of map with highlighted dots", {
  pdf(file=testfile)
  ## once with points
  expect_silent(MPplotmap(mpmap, highlight.points=rownames(mpmap)[c(1,3)]))
  ## once with labels
  expect_silent(MPplotmap(mpmap, highlight.points=rownames(mpmap)[c(1,3)],
                          label=TRUE))
  dev.off()
  unlink(testfile)
})

test_that("plot of map with different colors", {
  pdf(file=testfile)
  expect_silent(MPplotmap(mpmap, color=mapcolors))
  dev.off()
  unlink(testfile)
})

test_that("plot of maps with legend, colors, highlighted", {
  pdf(file=testfile)
  expect_silent(MPplotmap(mpmap, legend.separate=TRUE,
                          color=mapcolors,
                          highlight.points=rownames(mpmap)[1:4]))
  dev.off()
  unlink(testfile)
})



###############################################################################
## tests for plot with k colors


test_that("plot of scatter with k colors", {  
  mplayout1 = MPgetMap(mpsims[[1]])
  clusts = rep(c(1, 2), length(snames)/2)
  names(clusts) = snames
  
  pdf(file=testfile)
  expect_silent(MPplotScatterWithK(mplayout1, clusts))
  dev.off()
  unlink(testfile)
})




###############################################################################
## tests for plot with k colors


test_that("plot of scatter with links", {
  mpdist = mpsims[[1]]
  mplayout1 = MPgetMap(mpdist)
  clusts = rep(c(1, 2), length(snames)/2)
  names(clusts) = snames
  
  pdf(file=testfile)
  expect_silent(MPplotScatterWithLinks(mplayout1, mpdist))
  dev.off()
  unlink(testfile)
})



