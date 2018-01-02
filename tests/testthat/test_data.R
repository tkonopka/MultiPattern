## tests for datasets
## the package should come with attached datasets

cat("\ntest_data.R")




test_that("datasets 2S, 2M, 2L exist", {
  ## test dimentions of datasets
  expect_equal(dim(MPdata2S), c(72, 3))
  expect_equal(dim(MPdata2M), c(360, 3))
  expect_equal(dim(MPdata2L), c(720, 3))
  ## make sure all datasets have rownames
  expect_false(is.null(rownames(MPdata2S)))
  expect_false(is.null(rownames(MPdata2M)))
  expect_false(is.null(rownames(MPdata2L)))
})


test_that("datasets 3S, 3M, 3L exist", {
  ## The 3-series have different item counts from the other series
  ## test dimentions of datasets
  expect_equal(dim(MPdata3S), c(84, 3))
  expect_equal(dim(MPdata3M), c(420, 3))
  expect_equal(dim(MPdata3L), c(840, 3))
  ## make sure all datasets have rownames
  expect_false(is.null(rownames(MPdata3S)))
  expect_false(is.null(rownames(MPdata3M)))
  expect_false(is.null(rownames(MPdata3L)))
})


test_that("datasets 4S, 4M, 4L exist", {
  ## test dimentions of datasets
  expect_equal(dim(MPdata4S), c(72, 3))
  expect_equal(dim(MPdata4M), c(360, 3))
  expect_equal(dim(MPdata4L), c(720, 3))
  ## make sure all datasets have rownames
  expect_false(is.null(rownames(MPdata4S)))
  expect_false(is.null(rownames(MPdata4M)))
  expect_false(is.null(rownames(MPdata4L)))
})


test_that("datasets 6S, 6M, 6L exist", {
  ## test dimenstions of datasets
  expect_equal(dim(MPdata6S), c(72, 3))
  expect_equal(dim(MPdata9M), c(360, 3))
  expect_equal(dim(MPdata9L), c(720, 3))
  ## make sure all datasets have rownames
  expect_false(is.null(rownames(MPdata6S)))
  expect_false(is.null(rownames(MPdata6M)))
  expect_false(is.null(rownames(MPdata6L)))
})


test_that("datasets 9S, 9M, 9L exist", {
  ## test dimentions of datasets
  expect_equal(dim(MPdata9S), c(72, 3))
  expect_equal(dim(MPdata9M), c(360, 3))
  expect_equal(dim(MPdata9L), c(720, 3))
  ## make sure all datasets have rownames
  expect_false(is.null(rownames(MPdata9S)))
  expect_false(is.null(rownames(MPdata9M)))
  expect_false(is.null(rownames(MPdata9L)))
})

