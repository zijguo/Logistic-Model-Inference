library(MASS)
library(glmnet)
library(CVXR)

A1gen <- function(rho,p){
A1=matrix(0,p,p)
for(i in 1:p){
   for(j in 1:p){
     A1[i,j]<-rho^(abs(i-j))
   }
  }
  A1
}
n = 400
p = 500
mu <- rep(0,p)
rho = 0.5
Cov <- (A1gen(rho,p))/2
Cov2<-matrix(NA,nrow=p,ncol=p)
for(i in 1:p){
 for(j in 1:p){
   Cov2[i,j]<-0.5^(1+abs(i-j))
 }
}
beta <- rep(0,p)
beta[1:10] <- 0.5*c(1:10)/10
a0 = 0
set.seed(2020)
loading <- MASS::mvrnorm(1,mu,Cov2)
loading[20:p]<-loading[20:p]/5
true <- expo(t(loading)%*%beta+a0)
X <- MASS::mvrnorm(n,mu,Cov)
exp_val <- X%*%beta+a0
prob <- exp(exp_val)/(1+exp(exp_val))
y <- rbinom(n,1,prob)
Est <- LF_logistic(X = X, y = y, loading = loading, intercept = TRUE, weight = NULL)

### outputs of LF_logistic
### Est$prop.est: the bias-corrected estimator of the case probability
### Est$CI: confidence intervals for the case probability 
### Est$proj: the p-dim projection direction for bias correction
### Est$plug.in: the plug-in estimator of the case probability
### Est$case: 1 denotes that this observation is labelled as case (case probability is above 0.5); 0 denotes control (case probability is below 0.5)
