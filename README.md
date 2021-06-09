# Logistic-Model-Inference
This repository demonstrates the implementation of LiVE method proposed in 

_Guo Z, Rakshit P, Herman DS, Chen J (2019). “Inference for Case Probability in High-dimensional
Logistic Regression.”_

"LF_logistic" contains all the source functions to implement LiVE method.

"Comparison" contains source functions to implement other methods of inference namely :

(1) ps : It implements post-selection inference where it first selects important predictors through penalized logistic regression and then fit a standard logistic regression with the selected predictors. The plug-in estimate of the linear combination of the regression coefficients is constructed using the coefficient estimates from the final step.

(2) deb.hdi : It computes the plug-in estimator of the linear form using the co-ordinate debiased Lasso estimator. It uses the R package hdi.

(3) deb.WLDP : It implements weighted LDP algorithm in _Rong Ma, T Tony Cai, Hongzhe Li (2018). "Global and simultaneous hypothesis testing forhigh-dimensional logistic regression models."_ to compute the plug-in estimator of the linear form.

(4) Umethod : It performs a generalization of the transformation method in _Yinchu Zhu and Jelena Bradic (2018). "Linear hypothesis testing in dense high-dimensional linear models."_ to compute a debiased estimator of the case probability, the conditional probability of the binary response variable taking value 1 given the predictors are assigned to \code{loading}.

"Example_LF_logistic" gives a simulated data example on inference for the case probability in high dimensional logistic regression with the source file "LF_logistic" and "Comparison".
