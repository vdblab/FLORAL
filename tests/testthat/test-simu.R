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
})
