## tests for creating MP conigurations

cat("\ntest_new.R ")


###############################################################################
## Data objects

## samples
num.samples = 10
snames = paste0("S", 1:num.samples)


## object used while testing adding configurations
prep.list = list(one="D1", both=c("D1", "D2"))
dist.list = list(euc=dist.euclidean, man=dist.manhattan)


## object with multiple datasets and configurations
mplarge = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
MPaddConfig(mplarge, "confA", data.name="A", dist.fun=dist.list)
MPaddConfig(mplarge, "confB", data.name="B", dist.fun=dist.list)




###############################################################################
## Tests for creating new MultiPattern objects

test_that("create empty new MP", {
  ## completely empty analysis is not possible, must specify items
  expect_error(MPnew())
  mp=MPnew(snames)
  expect_equal(mp$items, snames)
})


test_that("create new MP with empty itemnames", {
  ## perhaps create configuration with empty items vector? Allow it
  mp=MPnew(c())
  expect_equal(mp$items, c())
})


test_that("new creates object with correct structure", {
  ## test structure of creted object
  mp=MPnew(snames)
  expect_equal(class(mp), "MultiPattern")
  ## object should have the main components
  components = c("data", "items", "configs", "settings")
  expect_equal(length(intersect(components, names(mp))), length(components))
  ## data list should be empty
  expect_equal(length(mp$data), 0)
  ## settings should be annotated with a class for pretty-print
  expect_equal(class(mp$settings), "MultiPatternSettings")  
})


test_that("create object and attach datasets", {
  ## test structure of creted object
  mp=MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  ## data list should be empty
  expect_equal(length(mp$data), 2)
  expect_equal(names(mp$data), c("A", "B"))
})


###############################################################################
## Tests for modifying MultiPattern objects


test_that("add more datasets into an MP configuration (all at once)", {
  mp = MPnew(snames)
  MPaddData(mp, list(Four=MPdata4S, Six=MPdata6S))
  expect_equal(length(mp$data), 2)
})


test_that("add more datasets into an MP configuration (one at a time)", {
  mp = MPnew(snames)
  MPaddData(mp, list(Four=MPdata4S))
  MPaddData(mp, list(Six=MPdata6S))
  expect_equal(length(mp$data), 2)
})


test_that("attempt add data in non-list format", {
  ## attempt to add in non-list
  mp = MPnew(snames)
  expect_error(MPaddData(mp, MPdata4S))
})


test_that("attempt add dataset with a repeat name", {
  mp = MPnew(snames)
  MPaddData(mp, list(A=MPdata4S))
  expect_error(MPaddData(mp, list(A=MPdata6S)))
})


test_that("add individual analysis configuration", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  MPaddConfig(mp, "confA", data.name="A")
  expect_equal(length(mp$configs), 1)
  MPaddConfig(mp, "confB", data.name="B")
  expect_equal(length(mp$configs), 2)
  expect_equal(names(mp$configs), c("confA", "confB"))
})


test_that("add family of analysis configurations via preproces list", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  MPaddConfig(mp, "confA", data.name="A", preprocess=prep.list)
  expect_equal(length(mp$configs), 2)
  expect_equal(names(mp$configs), c("confA.one", "confA.both"))
})


test_that("add family of analysis configurations via dist list", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  MPaddConfig(mp, "confA", data.name="A", dist.fun=dist.list)
  expect_equal(length(mp$configs), 2)
  expect_equal(names(mp$configs), c("confA.euc", "confA.man"))
})


test_that("add family of analysis configurations (incorrect)", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  ## cannot specify both preprocess and dist.fun as complex objects
  expect_error(MPaddConfig(mp, "confA", data.name="A",
                           preprocess=prep.list, dist.fun=dist.list))
})




###############################################################################
## Tests for removing configurations/data from MP objects


test_that("check objects for tests are well-formed", {
  ## this is a test that ensures that subsequent tests start with correct object
  expect_equal(length(mplarge$configs), 4)
  expect_equal(length(mplarge$data), 2)
})


test_that("remove a dataset", {
  ## this tests starts with a pre-made objects
  mpnow = mplarge
  ## remove one dataset (should also remove associated configurations)
  MPremove(mpnow, data="B")
  expect_equal(names(mpnow$data), "A")
  expect_equal(names(mpnow$configs), c("confA.euc", "confA.man"))
  ## removal of second dataset should leave empty object
  MPremove(mpnow, data="A")
  expect_equal(length(mpnow$configs)+length(mpnow$data), 0)
})


test_that("remove a configuration", {
  ## this tests starts with a pre-made objects
  mpnow = mplarge
  ## remove one configuration 
  MPremove(mpnow, config="confB.euc")
  expect_equal(names(mpnow$data), c("A", "B"))
  expect_equal(names(mpnow$configs), c("confA.euc", "confA.man", "confB.man"))
  ## removal of all configurations for a dataset does not remove the dataset
  MPremove(mpnow, config="confB.man")
  expect_equal(names(mpnow$data), c("A", "B"))
  expect_equal(names(mpnow$configs), c("confA.euc", "confA.man"))
})


test_that("remove a configuration (warnings, errors)", {
  mpnow = mplarge
  ## attempt to remove non-existing configuration is silent
  MPremove(mpnow, config="does-not-exist")
  expect_equal(mplarge, mpnow)
  ## attempt to remove dataset that does not exist should give error
  expect_warning(MPremove(mpnow, data="ZZ"))
})




###############################################################################
## Tests for changing settings within a configuration


test_that("remove a configuration (warnings, errors)", {
  mpnow = mplarge
  mysettings = list(num.PCs=5, some.other=0)
  ## attempt to add a new setting (some.other) should give warning
  expect_warning(MPchangeSettings(mpnow, mysettings))
  ## can turn warnings off
  mpnow = mplarge
  MPchangeSettings(mpnow, mysettings, warn=FALSE)
  ## object should reflect changes values
  expect_equal(mpnow$settings$num.PCs, 5)
  expect_equal(mpnow$settings$some.other, 0)
})




###############################################################################
## Tests for shortcuts for adding configurations (easyConfig)


test_that("create sets of configurations with easy config", {
  mp = MPnew(snames, data=list(A=MPdata4S))
  MPeasyConfig(mp, data="A", type=c("euclidean", "spearman", "hclust"))
  expect_gt(length(mp$configs), 6)
})


test_that("create sets of configurations with easy config with prefix", {
  mp = MPnew(snames, data=list(A=MPdata4S))
  MPeasyConfig(mp, data="A", type=c("hclust"), config.prefix="zz")
  ## this should have all configurations start with zz
  expect_equal(sum(grepl("^zz", names(mp$configs))), length(mp$configs))
  ## this should not have any rnorm configurations
  ##expect_equal(sum(grepl("rnorm", names(mp$configs))), 0)
})


test_that("create sets of configurations for different datasets", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  MPeasyConfig(mp, 
               type=list(A="hclust", B=c("euclidean", "pam")))
  ## this should have all configurations start with zz
  expect_gt(sum(grepl("clust", names(mp$configs))), 2)
  expect_equal(sum(grepl("A:clust.P", names(mp$configs))), 0)
  ## this should not have any rnorm configurations
  expect_gt(sum(grepl("B.clust.P", names(mp$configs))), 2)
  expect_equal(sum(grepl("B:clust.C", names(mp$configs))), 0)
})


test_that("attempt easy with strange input", {
  mp = MPnew(snames, data=list(A=MPdata4S, B=MPdata6S))
  ## attempt for non-exisinting datasets
  expect_error(MPeasyConfig(mp, data="ZZ", type="euclidean"))
  expect_error(MPeasyConfig(mp, 
                            type=list(ZZ="hclust", B=c("euclidean", "pam")))
               )
})

