context("max fréchet tree length response")

test_that("correct length leaf fréchet tree", {
  expect_length(Tmax(X=matrix(runif(200),100,2), Y=runif(100),id=sort(rep(c(1:10),10)),time=rep(c(1:10),10), toPlot="none")$feuilles, 100)
})
