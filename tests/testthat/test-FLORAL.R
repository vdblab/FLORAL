test_that("FLORAL() works", {
  
  set.seed(23420)
  
  dat <- simu(n=50,p=30,model="linear")
  expect_error(
    fit <- FLORAL(dat$xcount,dat$y,family="gaussian",progress=FALSE,step2=TRUE),
    NA
  )
  
  dat <- simu(n=50,p=30,model="binomial")
  expect_error(
    fit <- FLORAL(dat$xcount,dat$y,family="binomial",progress=FALSE,step2=TRUE),
    NA
  )
  
  dat <- simu(n=50,p=30,model="cox")
  expect_error(
    fit <- FLORAL(dat$xcount,survival::Surv(dat$t,dat$d),family="cox",progress=FALSE,step2=TRUE),
    NA
  )

  dat <- simu(n=50,p=30,model="poisson")
  expect_error(
    fit <- FLORAL(dat$xcount,dat$y,family="poisson",progress=FALSE,step2=TRUE),
    NA
  )

})

test_that("FLORAL() GEE Poisson works", {

  set.seed(23420)

  dat <- simu(n=50,p=30,model="gee",geetype="poisson",m=3,corstr="exchangeable")
  expect_error(
    fit <- FLORAL(x=cbind(dat$tvec, dat$xcount),y=dat$y,id=dat$id,family="poisson",
                  ncov=1,longitudinal=TRUE,intercept=TRUE,corstr="exchangeable",
                  lambda.min.ratio=1e-3,ncv=2,progress=FALSE,step2=FALSE),
    NA
  )

})
