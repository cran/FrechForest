context("Impurity")

test_that("Impurity works ", {
  expect_equal(length(impurity(list(type="scalar",Y=rnorm(100,0,1),id=c(1:100)))), 1)
})
