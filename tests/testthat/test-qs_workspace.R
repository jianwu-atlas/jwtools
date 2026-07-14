test_that("qs_save_workspace + qs_load_workspace round-trip works", {
  skip_if_not_installed("qs2")

  # 在一个独立的小环境里放几个简单变量，避免污染全局环境
  e <- new.env()
  e$a <- 1:10
  e$b <- "hello"
  e$tmp_big <- matrix(0, 2, 2)

  tmp_file <- tempfile(fileext = ".qs2")
  on.exit(unlink(tmp_file), add = TRUE)

  saved <- qs_save_workspace(
    tmp_file,
    envir = e,
    nthreads = 1,
    exclude_pattern = "^tmp_",
    verbose = FALSE
  )
  expect_setequal(saved, c("a", "b"))

  e2 <- new.env()
  restored <- qs_load_workspace(tmp_file, nthreads = 1, envir = e2, verbose = FALSE)
  expect_setequal(restored, c("a", "b"))
  expect_identical(e2$a, 1:10)
  expect_identical(e2$b, "hello")
  expect_false(exists("tmp_big", envir = e2))
})

test_that("qs_save_workspace errors when nothing to save", {
  e <- new.env()
  tmp_file <- tempfile(fileext = ".qs2")
  expect_error(qs_save_workspace(tmp_file, envir = e, verbose = FALSE))
})

test_that("qs_load_workspace errors on missing file", {
  expect_error(qs_load_workspace(tempfile(fileext = ".qs2")))
})
