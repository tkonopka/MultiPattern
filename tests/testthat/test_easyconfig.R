## tests for creating MP configurations using MPeasyConfig

cat("\ntest_easyconfig.R\n")


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
## Adding configs based on subspaces


test_that("easy euclidean/manhattan/canberra/spearman", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc",
               type=c("euclidean", "manhattan", "canberra", "spearman", "pearson"))
  expect_equal(length(mp$configs), 7)
})


test_that("easy subspace1", {
  ## should give a configuration corresponding to each feature
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="subspace1")
  expected = sort(paste0("abc:subspace1.", colnames(abc)))
  expect_equal(confNames(mp), expected)
})


test_that("easy subspace1 on only some features", {
  ## should give a configuration corresponding to each feature
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="subspace1", preprocess=c("D2", "D4"))
  expected = sort(paste0("abc:subspace1.", c("D2", "D4")))
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


test_that("easy subspace2", {
  ## should create one configuration for each pair of features 
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="subspace2", preprocess=c("D2", "D4", "D5"))
  expected = sort(paste0("abc:subspace2.", c("D2.D4", "D2.D5", "D4.D5")))
  expect_equal(confNames(mp), expected)
})


test_that("easy subspaceR finite number", {
  mp = MPnew(snames, data=list(abc=abc))
  mp$settings$subspace.num.random = 3
  mp$settings$subspace.d.random = 2
  ## should create a finite number of subspace-like configurations
  mp = MPeasyConfig(mp, data="abc", type="subspacer")
  expect_equal(length(mp$configs), 3)
})


test_that("easy subspaceR skip", {
  mp = MPnew(snames, data=list(abc=abc))
  mp$settings$subspace.num.random = 0
  mp$settings$subspace.d.random = 2
  ## should run the subspacer plugin, but because subspace.num.random=0, nothing should change
  mp = MPeasyConfig(mp, data="abc", type="subspacer")
  expect_equal(length(mp$configs), 0)
})




###############################################################################
## Adding pca-themed configs using easyConfig

## expected configurations based on 3 and 4 PCA components
pca2 = sort(c("PCA.2", "PC1", "PC2", "PC1.PC2"))
pca3 = sort(c("PCA.3", "PC1", "PC2", "PC3",
              "PC1.PC2", "PC1.PC3", "PC2.PC3"))
pca4 = sort(c("PCA.4", "PC1", "PC2", "PC3", "PC4",
              "PC1.PC2", "PC1.PC3", "PC1.PC4",
              "PC2.PC3", "PC2.PC4", "PC3.PC4"))

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


test_that("easy PCA with too small dataset", {
  ## make data with one column only (too small for ICA)
  mp = MPnew(snames, data=list(abc=abc[,1,drop=FALSE]))
  mp = MPeasyConfig(mp, data="abc", type="pca")
  mp = MPeasyConfig(mp, data="abc", type="euclidean")
  expected = "abc:euclidean"
  expect_equal(confNames(mp), expected)
})


###############################################################################
## Adding ica-themed configs using easyConfig

test_that("easy RPCA up to 3", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="rpca")
  expected = c(paste0("abc.S:", pca3), paste0("abc.L:", pca3))
  expect_equal(confNames(mp), sort(expected))
})


test_that("easy RPCA aborts when NA", {
  abcna = abc
  abcna[1,1] = NA
  mp = MPnew(snames, data=list(abc=abcna))
  mp = MPchangeSettings(mp, list(num.PCs=3))
  mp = MPeasyConfig(mp, data="abc", type="rpca")
  mp = MPeasyConfig(mp, data="abc", type="euclidean")
  expected = c("abc:euclidean")
  expect_equal(confNames(mp), sort(expected))
})




###############################################################################
## Adding ica-themed configs using easyConfig


test_that("easy ICA up to fraaction", {
  mp = MPnew(snames, data=list(abc=abc))
  ## ask for number of PCs to be half of the variance (here PC1 and PC2)
  mp = MPchangeSettings(mp, list(num.ICs=0.2))
  mp = MPeasyConfig(mp, data="abc", type="ica")
  expected = gsub("PC", "IC", paste0("abc:", pca2))
  expect_equal(confNames(mp), expected)
})


test_that("easy ICA with too small dataset", {
  ## make data with one column only (too small for ICA)
  mp = MPnew(snames, data=list(abc=abc[,1,drop=FALSE]))
  mp = MPeasyConfig(mp, data="abc", type="ica")
  mp = MPeasyConfig(mp, data="abc", type="euclidean")
  expected = "abc:euclidean"
  expect_equal(confNames(mp), expected)
})




###############################################################################
## Adding dbscan-themed configs using easyConfig

test_that("easy hamming", {
  num = 10
  data.num = abc[1:num,1:2]
  data.mixed = data.frame(A=1:nrow(data.num), B= c("a", "b"))
  data.chars = data.frame(X=c("a", "b"), Z=letters[1:num], stringsAsFactors=F)
  data.matrix = as.matrix(data.chars)
  rownames(data.matrix) = rownames(data.mixed) = rownames(data.chars) = rownames(data.num)  
  mp = MPnew(snames, data=list(numeric=data.num, mixed=data.mixed,
                               chars=data.chars, matrix=data.matrix))
  ## with numeric data, hamming should skip and do nothing
  mp = MPeasyConfig(mp, data="numeric", type="hamming")
  ## with mixed data, hamming should pick out character/factor columns
  mp = MPeasyConfig(mp, data="mixed", type="hamming")
  ## with char data, hamming should use all and subspaces
  mp = MPeasyConfig(mp, data="chars", type="hamming")
  ## also with char data in matrix form
  mp = MPeasyConfig(mp, data="matrix", type="hamming")  
  expected = c("mixed:hamming",
               paste0("chars:hamming", c("", ".X", ".Z")),
               paste0("matrix:hamming", c("", ".X", ".Z")))
  expect_equal(confNames(mp), sort(expected))
})




###############################################################################
## Adding dbscan-themed configs using easyConfig

test_that("easy dbscan", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPeasyConfig(mp, data="abc", type="dbscan")
  expected = paste0("abc:dbscan.", c(1,2,3,4))
  expect_equal(confNames(mp), sort(expected))
})




###############################################################################
## Adding nmf-themed configs using easyConfig

test_that("easy nmf", {
  mp = MPnew(snames, data=list(abc=abc))
  mp = MPchangeSettings(mp, list(nmf.rank=4))
  mp = MPeasyConfig(mp, data="abc", type="nmf")
  expected = paste0("abc:nmf", c(2,4))
  expect_equal(confNames(mp), sort(expected))
})


test_that("easy nmf when data has NULL", {
  abcna = abc
  abcna[1,1] = NA
  mp = MPnew(snames, data=list(abc=abcna))
  mp = MPchangeSettings(mp, list(nmf.rank=4))
  mp = MPeasyConfig(mp, data="abc", type="nmf")
  mp = MPeasyConfig(mp, data="abc", type="euclidean")
  ## nmf should abort, euclidean will continue (will use NA imputation)
  expected = c("abc:euclidean")
  expect_equal(confNames(mp), sort(expected))
})


test_that("easy nmf when rank is too low", {
  small = abc[,1, drop=FALSE]
  mp = MPnew(snames, data=list(small=small))
  mp = MPchangeSettings(mp, list(nmf.rank=4))
  mp = MPeasyConfig(mp, data="small", type="nmf")
  mp = MPeasyConfig(mp, data="small", type="euclidean")
  ## nmf should abort, euclidean will continue (will use NA imputation)
  expected = c("small:euclidean")
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


