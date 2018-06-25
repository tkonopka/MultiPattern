## tests for functions in MP_print.R
## (These are rudimentary tests that just check if an outut was produced or not)

cat("\ntest_print.R\n")


###############################################################################
## Data objects

## samples
num.samples = 8
snames = paste0("S", 1:num.samples)

## object with a minimal number of samples and configurations
mptest = MPnew(snames, data=list(A=MPdata4S[, 1:2]))
MPeasyConfig(mptest, type=c("euclidean", "manhattan"))
mpsims = MPgetDistances(mptest, verbose=FALSE)



###############################################################################
## Tests that print functions produce output


test_that("print general info on MP analysis", {
  ## some kind of header message
  expect_message(print.MultiPattern(mptest), "configuration object")
  ## lines with particular information
  expect_message(print.MultiPattern(mptest), "observations")
  expect_message(print.MultiPattern(mptest), "features")
  expect_message(print.MultiPattern(mptest), "analysis")
})

test_that("print info on MP settings", {
  expect_message(print.MultiPatternSettings(mptest$settings, "settings"))
})

test_that("print info on MP similarities", {
  expect_message(print.MultiPatternSimilarities(mpsims), "similarities")
})

test_that("printing with wrong obects", {
  expect_error(print.MultiPattern(1:3))
  expect_error(print.MultiPatternSettings(1:3))
  expect_error(print.MultiPatternSimilarities(1:3))
})

