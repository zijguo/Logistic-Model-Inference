######## Post Selection Method ###################

ps=function(X,y,loading)
{
  mod_1 <- cv.glmnet(X,y,family="binomial",alpha = 1)
  var <- unique(c(which(coef(mod_1,s="lambda.min")!=0)-1))
  mod <- glm(y~X[,var[-1]],family=binomial(link='logit'))
  loading_n <- c(1,loading[var[-1]])
  se <- sqrt(t(loading_n)%*%vcov(mod)%*%loading_n)
  ps_est=sum(loading_n*mod$coefficients)
  returnList<- list("ps.est"=as.vector(ps_est),"ps.se"=as.vector(se))
}

######### Plug-in Debiased using hdi ####################

deb.hdi=function(X,y,loading)
{
  fit.lasso <- hdi:::lasso.proj(x = X, y = y, standardize = F, family = "binomial", return.Z = T)
  pdata <- hdi:::prepare.data(x = X,y = y,standardize = F,family = "binomial")
  X <- pdata$x
  y <- pdata$y
  
  Z <- hdi:::score.rescale(Z = fit.lasso$Z, x = X)$Z
  Cov.est <- crossprod(Z)/(nrow(X)^2) 
  se <- sqrt(t(loading)%*%Cov.est%*%loading)
  
  hdi.est<-sum(loading*fit.lasso$bhat)
  returnList<- list("hdi.est"=as.vector(hdi.est),"hdi.se"=as.vector(se))
}

############# Plug-in Debiased using WLDP ####################

logistic.test <- function(x,y, nfolds=5, lambda = 0, 
                          tune.1 = 1.5, tune.2 = 1.01, 
                          intercept = F, fdr = 0.05){
  fi = function(x){
    exp(x)/(1+exp(x))
  }
  p = dim(x)[2]
  n = dim(x)[1]
  if(lambda == 0){
    logistic.cv = cv.glmnet(x = x, y = y, family = "binomial", alpha = 1,  
                            intercept=intercept, nfolds = nfolds, type.measure = "class")
    lambda = logistic.cv$lambda.min
  }
  my.logistic.fit = glmnet(x = x, y = y, family = "binomial", alpha = 1,  
                           intercept=intercept, lambda = lambda, standardize = F)
  b.hat = coef(my.logistic.fit)
  print("Estimation Succeeded!")
  
  
  W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try = matrix(ncol = p,nrow =5)
  tau.try = matrix(ncol = p,nrow = 5)
  
  V = matrix(ncol=n, nrow = p)
  tau = c()
  lambda.chosen=c()
  for(i in 1:p){
    nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for(lambda.i in 1:5){
      V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(sum((V[i,])^2*W.n1))))
      tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% X[,i])
    }
    zeta0 = sqrt(2*log(p))
    if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
    lambda.chosen[i] = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
    tau[i] = tau.try[lambda.chosen[i],i]
    lambda.chosen[i] = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
    tau[i] = tau.try[lambda.chosen[i],i]
    V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen[i]]
  }
  
  
  tau_cov=matrix(NA,p,p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      tau_cov[i,j]=sqrt(abs(sum((V[i,])*V[j,]*W.n1)/((V[i,]%*% X[,i])*(V[j,]%*% X[,j]))))
    }
  }
  
  V2 = t((t(V)*W.n1))
  #debaised estimator
  b.check = c()
  for(j in 1:p){
    b.check[j] = b.hat[j+1]+(V2[j,]%*%(y-fi(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
  }
  print("Bias Corrected!")
  
  M = b.check/tau
  
  return(list(b.check = b.check, M = M,cov=tau_cov))
}


deb.WLDP=function(X,y,loading)
{
  fit.deb.WLDP<-logistic.test(x=X,y=y)
  coef.est_WLDP<-fit.deb.WLDP$b.check
  deb.WLDP.est=sum(loading*coef.est.WLDP)
  Cov.est.WLDP <- fit.deb.WLDP$cov
  se.WLDP <- sqrt(t(loading)%*%Cov.est.WLDP%*%loading)
  returnList<- list("deb.WLDP.est"=as.vector(deb.WLDP.est),"deb.WLDP.se"=as.vector(se.WLDP))
}

############# Transformation Method #####################

Umethod = function(X,y,loading){
  H_star = loading%*%t(loading)/(sum(loading^2))
  D = diag(ncol(H_star)) - H_star
  
  eig = eigen(D)
  #eig$values
  
  P = eig$vectors
  P = P[,1:(ncol(H_star)-1)]
  
  Xt = matrix(NA, nrow=n, ncol = p)
  for(j in 1:n){
    Xt[j,1] = t(t(X[j,])%*%xnew/sum(xnew^2))
    Xt[j,-1]=t(t(X[j,])%*%P)
  }
  
  Est<-LF_logistic(Xt,y,loading=c(1,rep(0,(p-1))),weight = NULL)
}
