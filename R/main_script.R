# # unconstrained data
# 
# dat <- simulation()
# 
# x <- dat$x
# y <- dat$y0
# xy <- as.vector(t(y) %*% x)
# xx <- t(x) %*% x
# beta1 <- gd_cov_al(xx, xy, nrow(x), 1, rep(0,1000), 1, 0, TRUE)
# beta2 <- coordinate_descent_cov_updating_augmented_lagrangian(xx, xy, nrow(x), 1, rep(0,1000), 1, 0, FALSE)
# 
# beta3 <- gd_cov(xx, xy, nrow(x), 1, rep(0,1000))
# beta4 <- gd_naive(x,y,1,rep(0,1000))
# 
# fit <- linear_lasso(dat$x,dat$y0,100)
# plot(log(fit$lambda),fit$mse)
# plot(log(fit$lambda),fit$loss)
# 
# fit1 <- linear_lasso(dat$x,dat$y0,length.lambda = 100, mu = 10^2)
# plot(log(fit$lambda),fit$lossout)
# plot(log(fit$lambda),fit$mseout)
# 
# fit <- cv_linear_lasso(dat$x,dat$y0,length.lambda = 10)
# plot(log(fit$lambda),fit$losscvout)
# plot(log(fit$lambda),fit$msecvout)
# plot(log(fit$lambda),fit$mseout)
# which.min(fit$msecvout)
# 
# library(glmnet)
# 
# glmnetfit <- cv.glmnet(dat$x,dat$y0)
# 
# constrainedfit <- linear_lasso_augmented_lagrangian(dat$x,dat$y0,length.lambda = 100, mu = 10^5)
# 
# library(logratiolasso)
# 
# model_fit <- glmnet.constr(dat$x, dat$y0, family = "gaussian")
# 
# ### constrained data
# 
# dat <- simulation_constrained(intercept = TRUE)
# # ptm <- proc.time()
# # fit <- linear_lasso_al(dat$x,dat$y0,len=100,mu=1,ub=100,intercept=FALSE)
# # proc.time() - ptm
# 
# ptm <- proc.time()
# res1 <- LogRatioLasso(dat$x,dat$y0,intercept=FALSE,progress=TRUE)
# proc.time() - ptm
# ptm <- proc.time()
# res2 <- LogRatioLasso(dat$x,dat$y,intercept=TRUE)
# proc.time() - ptm
# plot(log(res2$lambda),res2$cvmse.mean)
# 
# ### binary data
# 
# dat <- simulation_constrained_binary(n=100,p=200)
# ptm <- proc.time()
# res3 <- LogRatioLogisticLasso(dat$x,dat$y)
# proc.time() - ptm
# plot(log(res3$lambda),res3$cvmse.mean)
# res3$best.beta0
# res3$best.beta
# 
# 
# binaryfit <- logratiolasso_binary(dat$x,dat$y,mu=1,ub=50,tol=1e-7)
# plot(log(binaryfit$lambda),binaryfit$lossout)
# plot(log(binaryfit$lambda),binaryfit$mseout)
# plot(log(binaryfit$lambda),binaryfit$GIC)
# View(binaryfit$betaout)
# which.min(binaryfit$GIC)
# colSums(binaryfit$betaout != 0)[which.min(binaryfit$GIC)]
# 
# cvbinaryfit <- cv_logratiolasso_binary(dat$x,dat$y,mu=1,ub=50,tol=1e-7)
# plot(log(cvbinaryfit$lambda),cvbinaryfit$losscvout)
# plot(log(cvbinaryfit$lambda),cvbinaryfit$msecvout)
# plot(log(cvbinaryfit$lambda),cvbinaryfit$GIC)
# View(cvbinaryfit$betaout)
# which.min(cvbinaryfit$msecvout[1:80])
# colSums(cvbinaryfit$betaout != 0)[which.min(cvbinaryfit$msecvout[1:80])]
# 
# 
# model_fit <- glmnet.constr(dat$x, dat$y, family = "binomial")
# cvfit <- cv.glmnet.constr(model_fit,dat$x, dat$y)
# sum(cvfit$beta)
# sum(cvfit$beta!=0)
# 
# 
# dat <- simulation_constrained_cox(n=100,p=200)
# ptm <- proc.time()
# coxfit <- LogRatioCoxLasso(dat$x,survival::Surv(dat$t,dat$d))
# proc.time() - ptm
# plot(log(coxfit$lambda),coxfit$cvdev.mean)

# library(survival)
# coxfit <- coxph(Surv(dat$t,dat$d)~dat$x[,1:10])
# sum(coxfit$coefficients)
# 
# coxfit <- logratiolasso_cox(dat$x,dat$d,dat$t,mu=1,ub=50,tol=1e-7)
# plot(log(coxfit$lambda),coxfit$lossout)
# plot(log(coxfit$lambda),coxfit$mseout)
# plot(log(coxfit$lambda),coxfit$GIC)
# View(coxfit$betaout)
# which.min(coxfit$lossout)
# colSums(coxfit$betaout != 0)[which.min(coxfit$mseout)]
# 
# cvcoxfit <- cv_logratiolasso_cox(dat$x,dat$d,dat$t,mu=1,ub=50,tol=1e-7)
# plot(log(cvcoxfit$lambda),cvcoxfit$losscvout)
# plot(log(cvcoxfit$lambda),cvcoxfit$msecvout)
# plot(log(cvcoxfit$lambda),cvcoxfit$GIC)
# View(cvcoxfit$betaout)
# which.min(cvcoxfit$msecvout)
# colSums(cvcoxfit$betaout != 0)[which.min(cvcoxfit$msecvout)]
