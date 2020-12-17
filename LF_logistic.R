getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

expo <- function(z){
  g = exp(z)/(1+exp(z))
  return(g)
}

Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  #col.norm <- 1 / sqrt((1 / n) * diagXtX(X, MARGIN = 2))
  Xnor <- X %*% diag(col.norm)
  
  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)
  
  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta)
  return(returnList)
}

Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (is.null(lambda)) {
    outLas <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else {
    outLas <- glmnet(X, y, family = "binomial", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

Direction_fixedtuning_logistic<-function(X,loading,mu=NULL,weight,deriv.vec){
  pp<-ncol(X)
  n<-nrow(X)
  if(is.null(mu)){
    mu<-sqrt(2.01*log(pp)/n)
  }
  loading.norm<-sqrt(sum(loading^2))
  
  if (loading.norm==0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  
  v<-Variable(pp+1)
  
  obj<-1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob<-Problem(Minimize(obj))
  result<-solve(prob)
  print("fixed mu")
  print(mu)
  opt.sol<-result$getValue(v)
  cvxr_status<-result$status
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  returnList <- list("proj"=direction)
  return(returnList)
}

Direction_searchtuning_logistic<-function(X,loading,mu=NULL,weight,deriv.vec,resol=1.5, maxiter=10){     #included weight and deriv.vec
  pp<-ncol(X)
  n<-nrow(X)
  tryno = 1;
  opt.sol = rep(0,pp+1);
  lamstop = 0;
  cvxr_status = "optimal";
  
  mu = sqrt(2.01*log(pp)/n);
  while (lamstop == 0 && tryno < maxiter){
    ###### This iteration is to find a good tuning parameter
    lastv = opt.sol;
    lastresp = cvxr_status;
    loading.norm<-sqrt(sum(loading^2))
    
    if (loading.norm==0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    
    v<-Variable(pp+1)
    obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    obj<-1/4*sum(((X%*%H%*%v)^2)*weight*deriv.vec)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))    #######modified
    prob<-Problem(Minimize(obj))
    result<-solve(prob)
    cvxr_status<-result$status
    if(tryno==1){
      if(cvxr_status=="optimal"){
        incr = 0;
        mu=mu/resol;
        
        opt.sol<-result$getValue(v)
        
        temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd<-sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
        temp.sd<-initial.sd
      }else{
        incr = 1;
        mu=mu*resol;
      }
    }else{
      if(incr == 1){ ### if the tuning parameter is increased in the last step
        if(cvxr_status=="optimal"){
          lamstop = 1;
          opt.sol <- result$getValue(v)
          
        }else{
          mu=mu*resol;
        }
      }else{
        if(cvxr_status=="optimal"&&temp.sd<3*initial.sd){
          mu = mu/resol;
          
          opt.sol <- result$getValue(v)
          
          temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd<-sqrt(sum(((X%*% temp.vec)^2)*weight*deriv.vec)/(n)^2)*loading.norm
        }else{
          mu=mu*resol;
          opt.sol=lastv;
          lamstop=1;
          tryno=tryno-1
        }
      }
    }
    tryno = tryno + 1;
  }
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step<-tryno-1
  print(step)
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

LF_logistic<-function(X,y,loading,weight=NULL,intercept=TRUE,init.Lasso=NULL,lambda=NULL,mu=NULL,step=NULL,resol = 1.5,maxiter=10,alpha = 0.05){
  
  xnew <- loading
  X<-as.matrix(X)
  p <- ncol(X);
  n <- nrow(X);
  n_y <- length(y)
  
  if(n_y!=n)
  {
    print("Check dimensions of X and y")
  }
  else
  {
    data = na.omit(data.frame(y,X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X);
    n <- nrow(X);
    col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
    Xnor <- X %*% diag(col.norm);
    if(is.null(init.Lasso))
    {
      init.Lasso <- Initialization.step(X,y,lambda,intercept)
      htheta <- init.Lasso$lasso.est
    }
    else
    {
      htheta <- init.Lasso
    }
    
    
    if (intercept==TRUE){
      Xb <- cbind(rep(1,n),Xnor);
      Xc <- cbind(rep(1,n),X);
      col.norm <- c(1,col.norm);
      pp <- (p+1);
    } else {
      Xb <- Xnor;
      Xc <- X;
      pp <- p
    }
    
    sparsity <- sum(abs(htheta) > 0.001)
    sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / max(0.9*n, n - sparsity))
    
    if(intercept==TRUE){
      loading=rep(0,pp)
      loading[1]=1
      loading[-1]=xnew
    }else{
      loading=xnew
    }
    loading.norm<-sqrt(sum(loading^2))
    lasso.plugin<-sum(loading*htheta)
    deriv.vec <- exp(Xc%*%htheta)/(1+exp(Xc%*%htheta))^2
    
    if(is.null(weight)){
      weight <- 1/deriv.vec
    }
    
    count=0
    for(i in 1:ncol(X)){
      if(length(unique(X[,i]))==1){
        count=count+1
      }
    }
    if(count!=0 && intercept==TRUE)
    {
      print("Data is singular")
    }else{
      #####################################################################################################
      ################## Correction step
      
      
      if ((n>=6*p)){
        gamma.hat <- (1/n)*sum((Xc^2)*weight*deriv.vec)   ##### computed gamma.hat instead of sigma.hat WRONG
        tmp <- eigen(gamma.hat)
        tmp <- min(tmp$values)/max(tmp$values)
      }else{
        tmp <- 0
      }
      
      if ((n>=6*p)&&(tmp>=1e-4)){
        direction <- solve(gamma.hat)%*%loading
      }else{
        if(is.null(step)){
          step.vec<-rep(NA,3)
          for(t in 1:3){
            index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
            Direction.Est.temp<-Direction_searchtuning_logistic(Xc[index.sel,],loading,mu=NULL,weight = weight[index.sel],deriv.vec = deriv.vec[index.sel],resol,maxiter)
            step.vec[t]<-Direction.Est.temp$step
          }
          step<-getmode(step.vec)
        }
        print(paste("step is", step))
        Direction.Est<-Direction_fixedtuning_logistic(Xc,loading,mu=sqrt(2.01*log(pp)/n)*resol^{-(step-1)},weight = weight,deriv.vec = deriv.vec)
        
        while(is.na(Direction.Est)&&(step>0)){
          step<-step-1
          Direction.Est<-Direction_fixedtuning_logistic(Xc,loading,mu=sqrt(2.01*log(pp)/n)*resol^{-(step-1)},weight = weight,deriv.vec = deriv.vec)
        }
        direction<-Direction.Est$proj
      }
      exp_pred=Xc%*%(htheta)
      
      weighed.residual=(y - exp(exp_pred)/(1+ exp(exp_pred)))*weight   
      
      correction = sum((Xc%*%direction)*weighed.residual)/n;
      debias.est=lasso.plugin+correction*loading.norm
      se<-sqrt(mean((Xc%*%direction)^2*weight^2*deriv.vec))*loading.norm/sqrt(n) 
      CI <- c(debias.est - 2*qnorm(1-alpha/2)*se, debias.est + 2*qnorm(1-alpha/2)*se)
      if(debias.est - qnorm(1-alpha)*se > 0){
        case<-1
      }else{
        case<-0
      }
      
      returnList <- list("prop.est" = expo(debias.est),
                         "CI" = c(expo(CI[1]),expo(CI[2])),
                         "proj"=direction,
                         "plug.in"=expo(lasso.plugin),
                         "case"=case
      )
      return(returnList)
    }
    
  }
}
