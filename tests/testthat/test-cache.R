test_that("package cache can set, get, drop, and clear values", {
  egvtools:::.egv_cache_clear()
  on.exit(egvtools:::.egv_cache_clear(), add = TRUE)

  expect_null(egvtools:::.egv_cache_get("x"))
  expect_invisible(egvtools:::.egv_cache_set("x", 42))
  expect_equal(egvtools:::.egv_cache_get("x"), 42)

  egvtools:::.egv_cache_drop("x")
  expect_null(egvtools:::.egv_cache_get("x"))

  egvtools:::.egv_cache_set("a", 1)
  egvtools:::.egv_cache_set("b", 2)
  egvtools:::.egv_cache_clear()
  expect_null(egvtools:::.egv_cache_get("a"))
  expect_null(egvtools:::.egv_cache_get("b"))
})
