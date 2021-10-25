#
# 15th April 2019
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
#

# Class
# mvMRInput:		input format for summary_mvMR_BF
# MRBF:			output format of summary_mvMR_BF



library(combinat)





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
# class MRBF for output
#




setClass("MRBF",
         representation(Exposure = "character",
                        Outcome = "character",
                        BMAve_Estimate = "numeric",
			BestModel_Estimate = "numeric",
			BestModel = "character",
			tupel = "character",
			pp="numeric",
			pp_marginal="numeric",
			betaX="matrix",
			betaY="matrix")
)









#
# summary_mvMR_BF
#

summary_mvMR_BF = function(object, sigma=0.5, prior_prob=0.5){


	#bX = object@betaX
	#bY = object@betaY
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
	#for (i in 1:5){

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
	if(max_evidence>308){
		rescale_diff = max_evidence - 308
		log_evidence = log_evidence - rescale_diff
		max_evidence = 308
	}
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
		pp_marginal = pp_marginal,
		betaX = bX,
		betaY= bY

	))

}















#
# compute the logBF for a set of risk factors given gamma (indicator vector)
#

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


#
# compute the causal estimates for a set of risk factors given gamma (indicator vector)
#


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


# input: y = betaY, x = betaX, sigma_vec = prior variance
# output: cooksD = cooks Distance, cooksD_thresh = suggested threshold based on a F distribution with df1=d,df2=k-d


cooksD = function(y,x,sigma_vec){
	k=length(y)
	d=ncol(x)
	x_mat=x
	H_fm = x_mat %*% solve(t(x_mat) %*% x_mat + sigma_vec^{-2} ) %*% t(x_mat)
	h_i = diag(H_fm)
	e = y-(H_fm%*%y)
	s_sq = c(1/(k-d) * t(e) %*% e)
	cooksD = e^2/(s_sq * d) * (h_i/(1-h_i)^2 )
	thresh=qf(0.5,df1=d,df2=k-d)
	return(list(cooksD=cooksD,cooksD_thresh=thresh))
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












#
# make table: report the best individual models
#


#input
#BMA_output object of class MRBF
#prior_sigma needs to be specified as in main function
#top = 10 how many top models to be reported?
#write.out = FALSE write out table as csv file
#csv.file.name = "best_model_out"



report.best.model = function(BMA_output, prior_sigma=0.5, top = 10, digits = 3, write.out = TRUE, csv.file.name="best_model_out"){


	if(class(BMA_output)[1] !="MRBF"){
		message("Input needs to be of class MRBF")
		break
		}

	pp=BMA_output@pp
	models=BMA_output@tupel
	betaX=BMA_output@betaX
	betaY=BMA_output@betaY
	rf=BMA_output@Exposure	
	sort_pp_model_object=sort.int(pp, index.return=TRUE, decreasing=TRUE)
	grep_rf=models[sort_pp_model_object$ix][1:top]

	rf_top = list()
	tupel_top = list()
	for(i in 1:top){
		tupel_top[[i]] =as.numeric(unlist(strsplit(grep_rf[i], ",")))
		rf_top[i]=paste(rf[as.numeric(unlist(strsplit(grep_rf[i], ",")))],  collapse=",")
	}

	theta_top = list()
	Theta=lapply(tupel_top, FUN = beta_gamma, y=as.matrix(betaY),x=as.matrix(betaX), sigma_vec=rep(0.5, ncol(as.matrix(betaX))))

	for(i in 1:top){
		Theta_iter = Theta[[i]]
		Theta_iter[Theta_iter!=0]
		theta_top[[i]]= paste(round(Theta_iter[Theta_iter!=0],digits=digits),  collapse=",")
	}

	best_models_out=cbind(rf_top, round(pp[sort_pp_model_object$ix][1:top],digits=digits), theta_top)
	colnames(best_models_out) = c("rf combination", "posterior probability", "causal estimate")


	if(write.out == TRUE){write.csv(best_models_out, file=csv.file.name)}

	return(best_models_out)

}





#
# make table: report the model-averaged results (MR-BMA)
#


#input
#BMA_output object of class MRBF
#top = 10 how many top models to be reported?
#write.out = FALSE write out table as csv file
#csv.file.name = "mr_bma_out"



report.mr.bma = function(BMA_output,  top = 10, digits = 3, write.out = TRUE, csv.file.name="mr_bma_out"){

	if(class(BMA_output)[1] !="MRBF"){
		message("Input needs to be of class MRBF")
		break
		}

	pp_marginal = BMA_output@pp_marginal
	bma=BMA_output@BMAve_Estimate
	rf=BMA_output@Exposure	
	sort_pp_object=sort.int(pp_marginal, index.return=TRUE, decreasing=TRUE)
	marginal_out=cbind(rf[sort_pp_object$ix][1:top], round(pp_marginal[sort_pp_object$ix][1:top],digits=digits), round(bma[sort_pp_object$ix][1:top],digits=digits) )
	colnames(marginal_out)=c("rf", "marginal inclusion", "average effect")

	if(write.out == TRUE){write.csv(marginal_out, file=csv.file.name)}

	return(marginal_out)

}







#
# permutation procedure to compute empirical p-values
# 1. create permutations
#

#input
#BMA_output: object of class mvMR_SSS output of summarymvMR_SSS
#nrepeat: number of permutations (ideally 100k or larger)
#save.matrix = True saves the permuted prior probabilities for each repetition as Rdata file with file.names as specified

#output
#permute_bma: matrix of size nrepeat times number of risk facttors


create.permutations = function(BMA_output, nrepeat = 100000, save.matrix=TRUE, file.name = "permutation_mrBMA"){
	
	#retain the datasets
	betaX_ivw = BMA_output@betaX
	betaY_ivw = BMA_output@betaY

	#get the parameters of MR-BMA
	kmin = BMA_output@kmin
	kmax = BMA_output@kmax 
	max_iter = BMA_output@max_iter
	sigma = BMA_output@sigma
	prior_prob = BMA_output@prior_prob
	rf = BMA_output@Exposure
	rs = as.character(1:nrow(BMA_output@betaX))

	#set up the matrix where to store the permuted marginal inclusion probabilities
	permute_bma = matrix(0, nrow=nrepeat, ncol=ncol(betaX_ivw))

	#run the permutations
	for(i in 1:nrepeat){
  		permute_it = sample(1:nrow(betaX_ivw))
  		betaY_ivw_permute= betaY_ivw[permute_it]
  		permutation_input=new("mvMRInput", betaX = as.matrix(betaX_ivw), betaY = as.matrix(betaY_ivw_permute), snps=rs, exposure=rf, outcome = "permutation")
  		BMA_output=summarymvMR_SSS(permutation_input,kmin=kmin,kmax=kmax, max_iter, sigma=sigma, prior_prob=prior_prob)
  		permute_bma[i,] = BMA_output@pp_marginal
	}

	if(save.matrix==TRUE)
	{
		save(permute_bma, file=file.name)
	}

	return(permute_bma)

}




#
# permutation procedure to compute empirical p-values
# 2. calculatte p-value
#

#input
#BMA_output: object of class mvMR_SSS output of summarymvMR_SSS
#permute_bma: matrix of size nrepeat times number of risk facttors which contains the permuted marginal inclusion probabilities

#output
#res_out: table with row names equal to the risk factors, marginal inclusion probabilities, empirical p-values and Benjamini-Hochberg false discovery rates adjusting for multiple testing 


calculate.p = function(BMA_output, permute_bma){
	
	#save the observed marginal inclusion probabilities and names of risk factors
	mip_observed = BMA_output@pp_marginal
	rf = colnames(BMA_output@betaX)

	#calculate the p-value
	p_val = rep(1,length(mip_observed))
	for(i in 1:length(mip_observed)){
		p_val[i]=(sum(permute_bma[,i]>mip_observed[i])+1)/(length(permute_bma[,i])+1)
	}

	p_adjust=p.adjust(p_val, "BH")
	res_out = cbind(mip_observed,p_val,p_adjust)
	row.names(res_out) = rf
	colnames(res_out) = c("mip", "pval", "fdr")
	return(res_out)

}


















helper = function(x){
  as.numeric(x[, drop =  F])
}





