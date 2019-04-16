#
# 11th April 2019
# Calculating the inverse-variance weighted estimate (and relevant statistics) using the multivariable input data 
# Based on the R package Mendelian Randomization
# https://cran.r-project.org/web/packages/MendelianRandomization/index.html
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


setClass("mvIVW",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",
                        Correlation = "matrix",
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
			Intercept = "numeric",
			Intercept_StdError = "numeric",
			Intercept_CILower = "numeric",
			Intercept_CIUpper = "numeric"
                        #SNPs = "character",
                        #Alpha = "numeric",
                        #Pvalue = "numeric",
                        #RSE = "numeric",
                        #Heter.Stat = "numeric")
	)
)








mvmr_ivw = function(object,
		model = "random",
		distribution = "normal",		
		correl = FALSE,
		alpha = 0.05,
		intercept = FALSE,
		...){

	bX = object@betaX
	bY = object@betaY
	#bXse = object@betaXse
	bYse = as.vector(object@betaYse)
	rho = object@correlation

	nsnps <- nrow(bX)
	nrf <- ncol(bX)

	#here we might check if there are enough SNPs to perform the multivariate analysis or if need regularised approaches

	if((correl == TRUE) &(is.na(sum(rho)))) {cat("Correlation matrix not given.")} 
	if(correl == FALSE){rho = diag(rep(1,nsnps))}
	else{rho=rho}


	if(length(bYse)==0){bYse=rep(1,nsnps)}
	omega <- bYse%o%bYse*rho

	if(intercept==FALSE) {
		bX = bX
		#bXse = bXse
		thetaIVW = solve(t(bX)%*%solve(omega)%*%bX)%*%t(bX)%*%(solve(omega)%*%bY)
		rse = bY - bX%*%thetaIVW
		}
	else {
		bX = cbind(rep(1,nsnps),bX)
		#bXse = c(0,bXse)
		thetaIVW = solve(t(bX)%*%solve(omega)%*%bX) %*%t(bX)%*%(solve(omega)%*%bY)
		rse = bY - bX%*%thetaIVW
		}



	if(model == "random") {thetaIVWse <- sqrt(diag(solve(t(bX)%*%solve(omega)%*%bX)))*max(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-nrf)),1)} 
	else if (model == "fixed"){thetaIVWse <- sqrt(diag(solve(t(bX)%*%solve(omega)%*%bX)))}

	if(distribution == "normal"){
		ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
		ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
		} 
	if (distribution == "t-dist"){
		ciLower <- ci_t("l", thetaIVW, thetaIVWse, nsnps - 1, alpha)
		ciUpper <- ci_t("u", thetaIVW, thetaIVWse, nsnps - 1, alpha)
		}

	#rse = sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-1))
	#heter.stat <- (nsnps - 1)*(rse^2)
	#pvalue.heter.stat <- pchisq(heter.stat, df = nsnps-1, lower.tail = F)
	#if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
	#if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=nsnps-1) }

	if(intercept == TRUE){
		thetaIVW = thetaIVW[2:(nrf+1)]
		thetaIVWse = thetaIVWse[2:(nrf+1)]
		ciLower = ciLower[2:(nrf+1)]
		ciUpper = ciUpper[2:(nrf+1)]
		theta_intercept = thetaIVW[1]
		interceptse = thetaIVWse[1]
		intercept_ciLower = ciLower[1]
		intercept_ciUpper = ciUpper[1]
		}
	else{
		theta_intercept=NaN
		interceptse=NaN
		intercept_ciLower = NaN
		intercept_ciUpper = NaN
		}

	return(new("mvIVW",
	Model = model,
	Exposure = object@exposure,
	Outcome = object@outcome,
	Correlation = object@correlation,
	Estimate = as.numeric(thetaIVW),
	StdError = as.numeric(thetaIVWse),
	CILower =  as.numeric(ciLower),
	CIUpper = as.numeric(ciUpper),
	Intercept = as.numeric(theta_intercept),
	Intercept_StdError = as.numeric(interceptse),
	Intercept_CILower = as.numeric(intercept_ciLower),
	Intercept_CIUpper = as.numeric(intercept_ciUpper)
	#SNPs = object@snps,
	#Pvalue = as.numeric(pvalue),
	#Alpha = alpha,
	#RSE = as.numeric(rse.corr),
	#Heter.Stat = c(heter.stat,pvalue.heter.stat)
	)
)


}







ci_normal <- function(type, mean, se, alpha){
  x <- 1 - alpha/2

  if(type == "l") return(mean - qnorm(x)*se)
  else if (type == "u") return(mean + qnorm(x)*se)
}


ci_t <- function(type, mean, se, df, alpha){
  x <- 1 - alpha/2

  if(type == "l") return(mean - qt(x, df = df)*se)
  else if (type == "u") return(mean + qt(x, df = df)*se)
}




