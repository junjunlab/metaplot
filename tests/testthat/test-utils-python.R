test_that("ensure_python_module reports and returns FALSE when python is unavailable", {
  testthat::local_mocked_bindings(
    py_config = function(...) invisible(NULL),
    py_available = function(...) FALSE,
    .package = "reticulate"
  )

  expect_message(
    result <- ensure_python_module("re"),
    "Please install python first"
  )
  expect_false(result)
})

test_that("ensure_python_module installs a missing module", {
  installed <- FALSE
  testthat::local_mocked_bindings(
    py_config = function(...) invisible(NULL),
    py_available = function(...) TRUE,
    py_module_available = function(...) FALSE,
    py_install = function(...) installed <<- TRUE,
    .package = "reticulate"
  )

  result <- ensure_python_module("re")

  expect_true(result)
  expect_true(installed)
})

test_that("ensure_python_module skips installation when the module is already present", {
  installed <- FALSE
  testthat::local_mocked_bindings(
    py_config = function(...) invisible(NULL),
    py_available = function(...) TRUE,
    py_module_available = function(...) TRUE,
    py_install = function(...) installed <<- TRUE,
    .package = "reticulate"
  )

  result <- ensure_python_module("re")

  expect_true(result)
  expect_false(installed)
})

test_that("ensure_python_module forwards a custom python path", {
  used_path <- NULL
  testthat::local_mocked_bindings(
    py_config = function(...) invisible(NULL),
    py_available = function(...) TRUE,
    use_python = function(python, ...) used_path <<- python,
    py_module_available = function(...) TRUE,
    .package = "reticulate"
  )

  ensure_python_module("re", pythonPath = "/usr/bin/python3")

  expect_equal(used_path, "/usr/bin/python3")
})
