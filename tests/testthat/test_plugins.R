## tests for functions in MP_plugins.R

cat("\ntest_plugins.R ")


###############################################################################
## Tests


test_that("listing and creating plugins", {
  ## check in-built plugins
  output = MPlistPlugins()
  expected = c("canberra", "dbscan", "euclidean", "manhattan", "spearman", 
               "subspace1", "subspace2", "subspacer", 
               "pca", "rpca",
               "hclust", "pam", "nmf")
  expect_equal(sort(output), sort(expected))

  ## create a new plugin (here dummy function)
  newplug.MultiPatternPlugin = function(mp, data.name) {
    
  }
  ## new listing should include a reference to newplug
  expect_true("newplug" %in% MPlistPlugins())
})


