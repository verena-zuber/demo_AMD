#
# 11th April 2019
# Calculating the lasso estimate using the multivariable input data 


# cv = True:  Estimate the regularisation parameter bY cross-validation
# cv = False: Parameter can be specified manually

#Note: Intercept is always true





library(glmnet)


#
# class mv-mr input
#

setClass("mvMRInput",
         representation(betaX = "matrix",
                        betaY = "matrix",
                        betaXse = "matrix",
                        betaYse = "matrix",
                        exposure = "character",
                        outcome = "character",
                        snps = "character",
                        effect_allele = "character",
                        other_allele  = "character",
                        eaf           = "numeric",                        
			correlation = "matrix")
)


#
# output class
#



setClass("MRl1",
         representation(Exposure = "character",
                        Outcome = "character",
                        Estimate = "numeric",
			Lambda = "numeric")
)




setClass("MRenet",
         representation(Exposure = "character",
                        Outcome = "character",
                        Estimate = "numeric",
			Lambda1 = "numeric",
			Lambda2 = "numeric")
)








glmnet_l1_mr = function(object, cv=TRUE, lambda=0.1, cv.param="lambda.1se"){


	bX = object@betaX
	bY = object@betaY

	if(cv==TRUE){
		cv = cv.glmnet(bX, bY, family = "gaussian", nfold = 10, type.measure = "mse", intercept=FALSE, alpha = 1,  standardize = FALSE)
		if(cv.param=="lambda.1se"){bestlambda=cv$lambda.1se}
		if(cv.param=="lambda.min"){bestlambda=cv$lambda.min}
		g.out =  glmnet(bX, bY, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = 1,  standardize = FALSE)
		l1.coeff = coef(g.out)[2:(ncol(bX)+1)]
	}

	else{
		bestlambda = lambda
		g.out =  glmnet(bX, bY, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = 1,  standardize = FALSE)
		l1.coeff = coef(g.out)[2:(ncol(bX)+1)]
	}

	return(new("MRl1", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		Estimate =  l1.coeff,
		Lambda = bestlambda
	))

}








glmnet_enet_mr = function(object, cv=TRUE, lambda=0.1, alpha=0.1, cv.param="lambda.1se"){


	bX = object@betaX
	bY = object@betaY

	if(cv==TRUE){
		a = seq(0.1, 0.9, 0.1)
		if(cv.param=="lambda.1se"){
			search = foreach(i = a, .combine = rbind) %do% {
  				cv = cv.glmnet(bX, bY, family = "gaussian", nfold = 10, type.measure = "mse", intercept=FALSE, alpha = i,  standardize = FALSE)
  				data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
			}
			cv = search[search$cvm == min(search$cvm), ]
			bestlambda = cv$lambda.1se
			bestalpha = cv$alpha
		}
		if(cv.param=="lambda.min"){
			search = foreach(i = a, .combine = rbind) %do% {
  				cv = cv.glmnet(bX, bY, family = "gaussian", nfold = 10, type.measure = "mse", intercept=FALSE, alpha = i,  standardize = FALSE)
  				data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
			}
			cv = search[search$cvm == min(search$cvm), ]
			bestlambda = cv$lambda.min
			bestalpha = cv$alpha
		}

		enet.out = glmnet(bX, bY, family = "gaussian", intercept=FALSE,lambda = bestlambda, alpha = bestalpha,  standardize = FALSE)
		enet.coeff = coef(enet.out)[2:(ncol(bX)+1)]
	}

	else{
		bestlambda = lambda
		bestalpha = alpha
		enet.out =  glmnet(bX, bY, family = "gaussian", intercept=FALSE, lambda = bestlambda, alpha = bestalpha,  standardize = FALSE)
		enet.coeff = coef(enet.out)[2:(ncol(bX)+1)]
	}

	return(new("MRenet", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		Estimate =  enet.coeff,
		Lambda1 = bestlambda,
		Lambda2 = bestalpha
	))

}













