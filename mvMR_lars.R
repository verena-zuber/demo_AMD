#
# 11th April 2019
# Calculating the lasso estimate using the multivariable input data 
# based on the lars implementation 
# additional test for significance using the covTest package


# input
# object of class mvMRInput
# cv = True:  estimate the regularisation parameter bY cross-validation
# cv = False: parameter (fraction [0,1]) can be specified manually 


# intercept=T in lars has the effect of centering the x variables and y variable. 
# It doesn't include an explicit intercept term with a coefficient.



library(lars)




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



setClass("MRlars",
         representation(Exposure = "character",
                        Outcome = "character",
                        Estimate = "numeric",
                        Pvals = "numeric",
			Fraction = "numeric")
)




lars_mr = function(object, cv=TRUE, intercept=FALSE, fraction=0.1){


	bX = object@betaX
	bY = object@betaY

	if(intercept==TRUE){lars_out = lars(x=bX,y=bY, type = "lasso", intercept=TRUE)}
	else{lars_out = lars(x=bX,y=bY, type = "lasso", intercept=FALSE)}

	if(cv==TRUE){
		lars.cv = cv.lars(x=bX,y=bY, type = "lasso", mode = "fraction", plot = FALSE)
		best.fraction=lars.cv$index[which.min(lars.cv$cv.error)]
		lars.coef=coef(lars_out,s=best.fraction,mode="fraction")
	}
	else{
		best.fraction = fraction
		lars.coef=coef(lars_out,s=best.fraction,mode="fraction")
	}
	



	return(new("MRlars", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		Estimate =  lars.coef,
		Fraction = best.fraction
	))

}









