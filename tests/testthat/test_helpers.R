## tests for functions in R/MP_helpers.R

cat("\ntest_helpers.R ")


###############################################################################
## Tests for subspaces


test_that("get subspace combinations (usual case)", {
  output = MPgetAll2Subspaces(letters[1:4])
  expected = list(a.b = c("a", "b"),
                  a.c = c("a", "c"),
                  a.d = c("a", "d"),
                  b.c = c("b", "c"),
                  b.d = c("b", "d"),
                  c.d = c("c", "d"))
  expect_equal(output, expected)
})

test_that("get subspace combinations (abnormal cases)", {
  expect_error(MPgetAll2Subspaces(letters[1]))
  expect_error(MPgetAll2Subspaces(NULL))
})

test_that("get subspace combinations (numeric)", {
  output = MPgetAll2Subspaces(1:3)
  expected = list("1.2" = c(1,2),
                  "1.3" = c(1,3),
                  "2.3" = c(2,3))
  expect_equal(output, expected)
})


test_that("get random subspaces", {
  output = MPgetRandomSubspaces(letters[1:8], 50)
  expect_equal(length(output), 50)
  ## ensure they are all different
  expect_mat = do.call(rbind, output)
  rownames(expect_mat) = NULL
  expect_equal(nrow(expect_mat), length(output))
  ## after so much sampling, should contain all the features
  out_items = sort(unique(unlist(output)))
  expect_equal(out_items, letters[1:8])
})


test_that("get random subspaces (special cases)", {  
  ## should not return duplicates (only four subspaces with 3/4)
  output4 = MPgetRandomSubspaces(letters[1:4], 10, d=3)
  expect_equal(length(output4), 4) 

  ## without oversampling, should sometimes get repeats and fail to return
  ## the appropriate number of subspaces
  numtries = 30
  out_wo = list()
  for (i in 1:numtries) {
    out_wo[[i]] = MPgetRandomSubspaces(letters[1:4], 4, d=2, oversample=1)
  }
  out_len = sum(sapply(out_wo, length))
  expect_lt(out_len, 4*numtries)

  ## with oversampling, should return appropriate number with near certainty
  out_w = list()
  for (i in 1:numtries) {
    out_wo[[i]] = MPgetRandomSubspaces(letters[1:4], 4, d=2, oversample=4)
  }
  out_len = sum(sapply(out_wo, length))
  expect_equal(out_len, 4*numtries)
})


test_that("get random subspaces can gives empty", {
  output = MPgetRandomSubspaces(letters[1:8], n=0)
  expected = list()
  expect_equal(output, expected)
})


test_that("get random subspaces gives error when requested d is too large", {
  expect_error(MPgetRandomSubspaces(letters[1:5], n=2, d=6))
})



###############################################################################
## Tests for top-feature selection

## create dataset with many constant featres, some varying features
mylen = 20
mydata = cbind(A=0,
               B=(1:mylen)/2,
               C=1:mylen,
               D=5,
               E=11,
               F=0.3,
               G=rep(c(0.01, 0.02, 0.03, 0.04), 5),
               H=rep(c(0.1, 0.2, 0.3, 0.4), each=5), 
               I=-100,
               J=999,
               K=(1:mylen)*(mylen:1)/1000
               )


test_that("identify feature that varies most", {
  ## create dataset with many constant features and one with some variation
  output = MPsuggestTopFeatures(mydata, ncomp=1)
  expected = c("B", "C")
  expect_equal(output, expected)
})

test_that("identify feature that varies most (fractional n)", {
  ## create dataset with many constant features and one with some variation
  output = MPsuggestTopFeatures(mydata, ncomp=0.05)
  expected = MPsuggestTopFeatures(mydata, ncomp=1)
  expect_equal(output, expected)
})



test_that("identify set of features", {
  ## create dataset with many constant features and one with some variation
  output = MPsuggestTopFeatures(mydata, ncomp=2)
  expected = c("B", "C", "H", "K")
  expect_equal(output, expected)
})

test_that("identify small set of features  ", {
  ## create dataset with many constant features and one with some variation
  output = MPsuggestTopFeatures(mydata, ncomp=function(x) { sqrt(x) } )
  expected = c("B", "C", "G", "H", "K")
  expect_equal(length(output), length(expected))
})




###############################################################################
## Tests for prep functions

test_that("make prep function using a function", {
  tempfun = function(x) {
    x+1
  }
  output = MPmakePrepFunction(tempfun)
  expect_equal(output, tempfun)
})

test_that("make prep function with wrong inputs", {
  expect_error(MPmakePrepFunction(c(1.2, 3.4)))
  myfactors = as.factor(letters[1:3])
  expect_error(MPmakePrepFunction(myfactors))
})

test_that("make prep function using character vector of features", {
  outfun = MPmakePrepFunction(c("D1", "D2"))
  outdata = outfun(MPdata4S)
  expected = MPdata4S[, c("D1", "D2")]
  expect_equal(outdata, expected)
})

test_that("make prep function using integer vector of features", {
  outfun = MPmakePrepFunction(1:2)
  outdata = outfun(MPdata4S)
  expected = MPdata4S[, c("D1", "D2")]
  expect_equal(outdata, expected)
})

test_that("make prep function using numeric vector with names", {
  myvec = c(D1=1.5, D2=0.5)
  outfun = MPmakePrepFunction(myvec)
  outdata = outfun(MPdata4S)
  expected = MPdata4S[, c("D1", "D2")]
  expected[, "D1"] = expected[, "D1"]*1.5
  expected[, "D2"] = expected[, "D2"]*0.5
  expect_equal(outdata, expected)
})




###############################################################################
## Tests for misc functions

test_that("creating square limits (simple case x)", {
  output= MPsquarelim(1:4, 1:5)
  expected = list(xlim=c(0.5, 4.5), ylim=c(1,5))
  expect_equal(output, expected)
})

test_that("creating square limits (simple case y)", {
  output= MPsquarelim(1:5, 1:4)
  expected = list(xlim=c(1, 5), ylim=c(0.5,4.5))
  expect_equal(output, expected)
})

test_that("creating square limits (some negative)", {
  output = MPsquarelim(1:4, -c(1,5))
  expected = list(xlim=c(0.5, 4.5), ylim=c(-5,-1))
  expect_equal(output, expected)
})


test_that("test cmd/tsne map", {
  d4 = dist(MPdata4S[, 1:2])
  output_cmd = MPgetMap(d4)
  expect_equal(dim(output_cmd), c(nrow(MPdata4S), 2))
  output_tsne = MPgetMap(d4, tsne=TRUE)
  expect_equal(dim(output_tsne), c(nrow(MPdata4S), 2))
})


test_that("normalization by column median/sd", {
  mm = as.matrix(cbind(A=c(1,5,9), B=c(12,15,16)))
  output = MPmatrixColNorm(mm)
  expected = mm
  expected[,1] = (mm[,1]-5)/sd(mm[,1])
  expected[,2] = (mm[,2]-15)/sd(mm[,2])
  expect_equal(output, expected)
})


test_that("randomize a matrix without replacement", {
  mm = matrix(1:100, ncol=10)
  output = MPrandomizeMatrix(mm)
  ## output should have as many unique element as input
  out.u = length(unique(as.vector(output)))
  expected.u = length(unique(as.vector(mm)))
  expect_equal(out.u, expected.u)
})

test_that("randomize a matrix with replacement", {
  mm = matrix(1:100, ncol=10)
  output = MPrandomizeMatrix(mm, perm.method="bootstrap")
  ## output is likely to pick an element twice,
  ## so number of unique elements should decrease
  out.u = length(unique(as.vector(output)))
  expected.u = length(unique(as.vector(mm)))
  expect_lt(out.u, expected.u)
})

test_that("randomize repeatedly should give different result", {
  mm = matrix(1:25, ncol=25)
  out1 = MPrandomizeMatrix(mm)
  out2 = MPrandomizeMatrix(mm)
  expect_gt(sum(abs(out1-out2)), 0)
})




## tests for MPcutree

