## tests for creating MP configurations using MPeasyConfig

cat("\ntest_easyconfig.R ")


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

## extract configuration names
confNames = function(mp) {
  result = names(mp$configs)
  result = result[!grepl("^rnorm", result)]
  sort(result)
}




###############################################################################
## Tests for adding configs using easyConfig


test_that("easy euclidean/manhattan/canberra/spearman", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc",
               type=c("euclidean", "manhattan", "canberra", "spearman"))
  expect_equal(length(mp$configs), 4)
})


test_that("easy subspace1", {
  ## should give a configuration corresponding to each feature
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="subspace1")
  expected = sort(paste0("abc:subspace1.", colnames(abc)))
  expect_equal(confNames(mp), expected)
})


test_that("easy subspace2", {
  ## should create one configuration for each pair of features 
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="subspace2")
  temp = expand.grid(colnames(abc), colnames(abc), stringsAsFactors=F)
  temp = temp[temp[,1] < temp[,2],]
  temp = apply(temp, 1, paste, collapse=".")
  expected = sort(paste0("abc:subspace2.", temp))
  expect_equal(confNames(mp), expected)
})




###############################################################################
## Tests for adding pca-themed configs using easyConfig

## expected configurations based on 3 and 4 PCA components
pca2 = sort(c("PC1", "PC1.PC2"))
pca3 = sort(c("PC1", "PC1.PC2", "PC1.PC3", "PC2.PC3", "PC1..PC3"))
pca4 = sort(c("PC1", "PC1.PC2", "PC1.PC3", "PC1.PC4",
              "PC2.PC3", "PC2.PC4", "PC3.PC4", "PC1..PC3", "PC1..PC4"))

test_that("easy PCA up to 4", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="pca")
  expected = paste0("abc:", pca4)
  expect_equal(confNames(mp), expected)
})


test_that("easy PCA up to 3", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="pca")
  expected = paste0("abc:", pca3)
  expect_equal(confNames(mp), expected)
})


test_that("easy PCA up to fraaction", {
  mp = MPnew(snames, data=list(abc=abc))
  ## ask for number of PCs to be half of the variance (here PC1 and PC2)
  mp = MPchangeSettings(mp, list(num.PCs=0.5))
  mp = MPeasyConfig(mp, data="abc", type="pca")
  expected = paste0("abc:", pca2)
  expect_equal(confNames(mp), expected)
})


test_that("easy robust PCA up to 3", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="rpca")
  expected = c(paste0("abc.rpcaS:", pca3), paste0("abc.rpcaL:", pca3))
  expect_equal(confNames(mp), sort(expected))
})




###############################################################################
## Tests for adding configurations to multiple datasets

test_that("easyconfig on two datasets", {
  mp = MPnew(snames, data=list(A=abc, B=abc))
  mp = MPeasyConfig(mp, type=list(A="euclidean", B="manhattan"))
  expected = c("A:euclidean", "B:manhattan")
  expect_equal(confNames(mp), expected)
})


test_that("easyconfig cannot understand data name and type in a list", {
  mp = MPnew(snames, data=list(A=abc, B=abc))
  expect_error(MPeasyConfig(mp, data="A", type=list(A="euclidean", B="manhattan")))
})


test_that("easyconfig type list must have names", {
  mp = MPnew(snames, data=list(A=abc, B=abc))
  expect_error(MPeasyConfig(mp, type=list("euclidean", "manhattan")))
})




###############################################################################
## Tests for adding clustering-themed configs using easyConfig

test_that("easy hclust", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="hclust")
  expected = sort(c("abc:clust.C2reg", "abc:clust.C2alt",
                    "abc:clust.A2reg", "abc:clust.A2alt",
                    "abc:clust.S2reg", "abc:clust.S2alt",
                    "abc:clust.C3reg", "abc:clust.C3alt",
                    "abc:clust.A3reg", "abc:clust.A3alt",
                    "abc:clust.S3reg", "abc:clust.S3alt"))  
  expect_equal(confNames(mp), expected)
})

test_that("easy pam", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="pam")
  expected = sort(c("abc:clust.P2reg", "abc:clust.P2alt",
                    "abc:clust.P3reg", "abc:clust.P3alt"))
  expect_equal(confNames(mp), expected)
})


