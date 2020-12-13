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
n = 100
p = 400
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
beta[1:10] <- c(1:10)/5
X <- MASS::mvrnorm(n,mu,Cov)
exp_val <- X%*%beta
prob <- exp(exp_val)/(1+exp(exp_val))
y <- rbinom(n,1,prob)
loading <- MASS::mvrnorm(1,mu,Cov2)
Est <- LF_logistic(X = X, y = y, loading = loading, intercept = TRUE, weight = NULL)