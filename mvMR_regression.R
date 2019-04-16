#
# 11th April 2019
# Calculating the inverse-variance weighted estimate using the multivariable input data based on a weighted linear model
#







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




setClass("mvREG",
         representation(
                        Estimate = "numeric",
                        StdError = "numeric",
                        Pval = "numeric")
	)






mvmr_wreg = function(object){

	bX = object@betaX
	bY = object@betaY
	#bXse = object@betaXse
	bYse = as.vector(object@betaYse)
	nsnps = nrow(bX)
	nrf = ncol(bX)
	if(length(bYse)==0){bYse=rep(1,nsnps)}
	lm_ivw = lm(bY ~ bX -1, weights= bYse^(-2)) 
	theta = summary(lm_ivw)$coefficients[,1]
	thetase = summary(lm_ivw)$coefficients[,2]	
	preg = summary(lm_ivw)$coefficients[,4]

	return(new("mvREG",
	Estimate = as.numeric(theta),
	StdError = as.numeric(thetase),
	Pval = as.numeric(preg)



	)
)


}






