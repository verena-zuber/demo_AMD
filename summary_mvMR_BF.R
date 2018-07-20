#
# 12th April 2018
# risk factor selection for multivariable MR based on Bayes Factors
# computation is based on summary data which speeds up computation 
#
# functions
# summary_mvMR_BF: 	exhaustive evaluation of all combinations of risk factors
# logBF_summary: 	compute log10 Bayes factor for indicator on summary data
# beta_summary: 	compute causal estimate for indicator on summary data
# cooksD: 		Cooks distance
# beta_est, beta_gamma: compute causal estimate on data matrix
# logBF, logBF_gamma:	compute log10 Bayes factor for indicator on data matrix 

# now synced with manuscript and BF derivation
# input: object of class mvMRInput
#

library(combinat)


summary_mvMR_BF = function(object, sigma=0.5, prior_prob=0.5){


	bX = object@betaX
	bY = object@betaY
	n_SNPs = nrow(bX)
	n_RF = ncol(bX)
	XtX = t(bX) %*% bX		# dim  d=10  x d=10	 	
 	XtY = t(bX) %*% bY		# dim  d=10  x 1
	YtY = t(bY) %*% bY

	tupel_all = matrix(character(0),0,0)
	logBF_all = matrix(numeric(0), 0,0) 
	logprior_all = matrix(numeric(0), 0,0) 
	theta_all = matrix(numeric(0), (ncol(bX)),0)
	sigma_vec=rep(sigma, ncol(bX))

	for (i in 1:ncol(bX)){

		combi=matrix(combn(1:ncol(bX), i),nrow=i)
		tupel=apply(combi, 2, paste, collapse=",") 
		tupel_all = c(tupel_all, tupel)
		logBF = apply(combi, MARGIN=2, FUN = logBF_summary, XtY=XtY,XtX=XtX,YtY=YtY, sigma_vec=sigma_vec, n=n_SNPs)
		logBF_all = c(logBF_all,logBF)
		log_prior=log10(prior_prob^i*(1-prior_prob)^(n_RF-i))
		logprior_vec = rep(log_prior,length(logBF)) 
		logprior_all = c(logprior_all,logprior_vec)
		theta_est = apply(combi, MARGIN=2, FUN = beta_summary, XtY=XtY,XtX=XtX,sigma_vec=sigma_vec)
		theta_all = cbind(theta_all,theta_est) 
	}

	log_evidence = logBF_all + logprior_all

	#best model
	best_model = tupel_all[which.max(log_evidence)]
	best_model_est = theta_all[,which.max(log_evidence)]
	#posterior for each model
	max_evidence = max(log_evidence) 
	sum_calib = sum(10^(log_evidence-max_evidence)) * 10^max_evidence
	pp=10^log_evidence/sum_calib 
	pp_mat=matrix(pp, ncol=length(pp), nrow=ncol(bX), byrow=TRUE)
	#model averaging causal effect
	BMA_estimate=rowSums(pp_mat * theta_all)
	#marginal inclusion
	index_mat = matrix(0, ncol=length(pp), nrow=ncol(bX))
	index=lapply(unlist(lapply(tupel_all, FUN=strsplit, split=","), recursive=FALSE), FUN=helper)
	for(i in 1:length(pp)){index_mat[unlist(index[i]),i] = 1}
	pp_marginal =rowSums(pp_mat * index_mat)
	bma_coef = BMA_estimate
	best_model_coef = best_model_est


	return(new("MRBF", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		BMAve_Estimate = bma_coef,
		BestModel_Estimate = best_model_coef,
		BestModel = best_model,
		tupel = tupel_all,
		pp=pp,
		pp_marginal = pp_marginal

	))

}





setClass("MRBF",
         representation(Exposure = "character",
                        Outcome = "character",
                        BMAve_Estimate = "numeric",
			BestModel_Estimate = "numeric",
			BestModel = "character",
			tupel = "character",
			pp="numeric",
			pp_marginal="numeric" )
)

















# compute the logBF for a set of risk factors given gamma (indicator vector)

logBF_summary = function(XtY,XtX,YtY,sigma_vec, gamma, n){
	
	#print(gamma)
	XtY_g = XtY[gamma]
	XtX_g = XtX[gamma, gamma] 
	sigma_g = sigma_vec[gamma]	

	invnu=diag(sigma_g^-2, ncol=length(gamma))
	invOmega=(invnu + XtX_g)
	B=solve(invOmega)%*%XtY_g 
	logBF=(-0.5*log10(det(invOmega))  - log10(prod(sigma_g)) - (n/2) *(log10(YtY - t(B)%*%invOmega%*%B)-log10(YtY)))
	
	return(logBF)

}



# compute the causal estimates for a set of risk factors given gamma (indicator vector)

beta_summary = function(XtY=XtY,XtX=XtX,sigma_vec, gamma){

	XtY_g = XtY[gamma]
	XtX_g = XtX[gamma, gamma] 
	sigma_g = sigma_vec[gamma]

	invnu=diag(sigma_g^-2, ncol=length(gamma))
	theta_out = rep(0, nrow(XtY))
	invOmega=(invnu + XtX_g)
	B=solve(invOmega)%*%XtY_g 
	theta_out[gamma] = B

	return(theta_out)
}








# 
# Cooks distance
#

cooksD = function(y,x,sigma_vec){
	k=length(y)
	d=ncol(x)
	x_mat=x
	H_fm = x_mat %*% solve(t(x_mat) %*% x_mat + sigma_vec^{-2} ) %*% t(x_mat)
	h_i = diag(H_fm)
	e = y-(H_fm%*%y)
	s_sq = c(1/(k-d) * t(e) %*% e)
	cooksD = e^2/(s_sq * d) * (h_i/(1-h_i)^2 )

	return(cooksD)
}










#
# further functions to compute beta based on data matrix
#

# compute the logBF for a set of risk factors 

#y outcome
#x predictors (genetic associations with risk factors)
#xmat designmatrix
#sigma_vec vector or prior variances

logBF = function(y,x,sigma_vec){
	
	k=length(y)
	xmat=x
	invnu=diag(sigma_vec^-2, ncol=ncol(xmat))
	invOmega=(invnu + t(xmat)%*% xmat)
	#invOmega0=k
	B=solve(invOmega)%*%t(xmat)%*%y 
	logBF=(-0.5*log10(det(invOmega)) - log10(prod(sigma_vec)) - (k/2) *(log10(t(y)%*% y - t(B)%*%invOmega%*%B)-log10(t(y)%*%y)))
	return(logBF)
}



logBF_gamma = function(y,x,sigma_vec, gamma){
	
	k=length(y)
	sigma_g = sigma_vec[gamma]
	xmat=x[,gamma]
	invnu=diag(sigma_g^-2, ncol=length(gamma))
	invOmega=(invnu + t(xmat)%*% xmat)
	B=solve(invOmega)%*%t(xmat)%*%y 
	logBF=(-0.5*log10(det(invOmega))  - log10(prod(sigma_g)) - (k/2) *(log10(t(y)%*% y - t(B)%*%invOmega%*%B)-log10(t(y)%*%y)))
	
	return(logBF)

}




# compute the causal estimates for a set of risk factors given gamma (indicator vector)

beta_gamma = function(y,x,sigma_vec, gamma){

	k=length(y)
	sigma_g = sigma_vec[gamma]
	xmat=x[,gamma]
	invnu=diag(sigma_g^-2, ncol=length(gamma))
	theta_out = rep(0, ncol(x))
	invOmega=(invnu + t(xmat)%*% xmat)
	B=solve(invOmega)%*%t(xmat)%*%y 
	theta_out[gamma] = B	

	return(theta_out)
}


# compute the causal estimates for a set of risk factors

beta_est = function(y,x,sigma_vec, intercept){

	k=length(y)
	xmat=x
	invnu=diag(sigma_vec^-2, ncol=ncol(x))
	theta_out = rep(0, ncol(x))
	invOmega=(invnu + t(xmat)%*% xmat)
	B=solve(invOmega)%*%t(xmat)%*%y 
	theta_out = B

	return(theta_out)
}






