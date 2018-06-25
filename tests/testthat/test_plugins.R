## tests for functions in MP_plugins.R

cat("\ntest_plugins.R\n")


###############################################################################
## Tests


test_that("listing and creating plugins", {
  ## check in-built plugins
  output = MPlistPlugins()
  expected = c("canberra", "dbscan", "euclidean", "hamming", "manhattan",
               "spearman", "pearson",
               "subspace1", "subspace2", "subspacer", 
               "ica", "pca", "rpca", "nmf", "hclust", "pam")
  expect_equal(sort(output), sort(expected))
})


test_that("new custom plugin is detected", {
  ## create a new plugin (here dummy function)
  newplug.MultiPatternPlugin = function(mp, data.name) {
    
  }
  ## new listing should include a reference to newplug
  expect_true("newplug" %in% MPlistPlugins())
})

