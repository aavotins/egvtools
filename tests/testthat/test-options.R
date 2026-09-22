test_that(".onLoad supplies package option defaults", {
  old_plan <- getOption("egvtools.future_plan")
  old_progress <- getOption("egvtools.progress")
  on.exit({
    options(egvtools.future_plan = old_plan)
    options(egvtools.progress = old_progress)
  }, add = TRUE)

  options(egvtools.future_plan = NULL, egvtools.progress = NULL)
  egvtools:::.onLoad(NULL, NULL)

  expect_equal(getOption("egvtools.future_plan"), "sequential")
  expect_true(getOption("egvtools.progress"))
})
