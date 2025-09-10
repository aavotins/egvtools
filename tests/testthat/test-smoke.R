test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
test_that("package loads", {
  expect_true(requireNamespace("egvtools", quietly = TRUE))
})
