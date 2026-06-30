test_that("simu() works", {
  expect_error(
    dat <- simu(n=50,p=30,model="linear"),
    NA
  )
  expect_error(
    dat <- simu(n=50,p=30,model="binomial"),
    NA
  )
  expect_error(
    dat <- simu(n=50,p=30,model="cox"),
    NA
  )
  expect_error(
    dat <- simu(n=50,p=30,model="poisson"),
    NA
  )
  expect_error(
    dat <- simu(n=50,p=30,model="gee",geetype="poisson",m=3,corstr="exchangeable"),
    NA
  )
})
