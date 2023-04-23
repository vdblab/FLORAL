test_that("FLORAL() works", {
  
  set.seed(23420)
  
  dat <- simu(n=50,p=100,model="linear")
  expect_error(
    fit <- FLORAL(dat$xcount,dat$y,family="gaussian",progress=FALSE,step2=TRUE),
    NA
  )
  
  dat <- simu(n=50,p=100,model="binomial")
  expect_error(
    fit <- FLORAL(dat$xcount,dat$y,family="binomial",progress=FALSE,step2=TRUE),
    NA
  )
  
  dat <- simu(n=50,p=100,model="cox")
  expect_error(
    fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d),family="cox",progress=FALSE,step2=TRUE),
    NA
  )
  
})
