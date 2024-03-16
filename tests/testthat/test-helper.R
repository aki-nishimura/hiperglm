test_that("closeness check works", {
  # Test absolute tolerance check
  abs_err <- 1e-3
  a <- c(10, 10) # Any reasonable value > 1 should work
  b <- a
  b[1] <- b[1] + abs_err
  expect_true(are_all_close(a, b, abs_tol = 2 * abs_err, rel_tol = abs_err))
  expect_false(are_all_close(a, b, abs_tol = abs_err / 2, rel_tol = abs_err))

  # Test relative tolerance check
  rel_err <- 1e-3
  a <- c(0.1, 0.1) # Any reasonable value < 1 should work
  b <- a
  b[1] <- b[1] * (1 + rel_err)
  expect_true(are_all_close(a, b, abs_tol = rel_err, rel_tol = 2 * rel_err))
  expect_false(are_all_close(a, b, abs_tol = rel_err, rel_tol = rel_err / 2))
})
