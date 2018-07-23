#
# 12th April 2018
# Risk factor selection for multivariable MR based on Bayes Factors
# Functions for simulating data for multivariable MR analysis

# eval.score:		ROC characteristics
# sim.mvnorm: 		fast multivariate normal
# generate_mv_input:	simulate synthetic input for simulation study
# generate_mv_input2: 	based on real risk factor associations simulate outcome associations







library(corpcor)




####
# class mv-mr input


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





# ROC characteristics and more for a classifier score vs the know truth

eval.score=function(score, truth){

	index=sort.int(score, decreasing=TRUE, index.return=TRUE)$ix
	truth[index]

	tp = cumsum(truth[index])
	fp = cumsum(1-truth[index])

	fn = sum(truth) - tp
	tn = sum(1-truth)- fp

	return(list(tp=tp, fp=fp, tn=tn, fn=fn, ppv = tp/(tp+fp), tpr=tp/(tp+fn), fpr=fp/(fp+tn)))
}








# simulate multivariate normal data
# mu:  mean
# sd:  standard deviation
# cR:  chol(R)

sim.mvnorm = function(n, mu, sd, cR)
{
   retval = matrix(rnorm(n * ncol(cR)), nrow = n) # independent unit normal
   retval = retval %*% cR                         # correlation
   retval = sweep(retval, 2, sd, "*")             # scale
   retval = sweep(retval, 2, mu, "+")             # location

   return(retval)
}





####
# simulate mv mr data outputs object of class mvMRInput
# parameters
###
# k: number of genetic variants
# theta: true causal association
# mu_beta: expectation of the genetic associations with the risk factors
# Sigma_beta: covariance between the genetic associations with the risk factors
# R2: proportion of variance of the outcome explained by the risk factors
# alpha=list(mean=0,sd=0): pleiotropy term



generate_mv_input = function(k,theta,mu_beta,Sigma_beta,R2,scale=TRUE,alpha=list(mean=0,sd=0)){
 
	d=length(mu_beta) # number of variables
	#check that dimension of theta, mu_beta and Sigma_beta matches
	dc=decompose.cov(Sigma_beta)
	cR=chol(dc$r)
	sd_beta = sqrt(dc$v)
	b_mat=sim.mvnorm(k, mu=mu_beta, sd=sd_beta, cR=cR)
	if(scale==TRUE) b_mat=scale(b_mat)
	else beta = b_mat
	cov_beta=cov(b_mat)
	#pleiotropy intercept
	if(alpha$mean==0 & alpha$sd==0) alpha_vec = rep(0,k)
	else alpha_vec = rnorm(k, mean=alpha$mean, sd = alpha$sd)
	#fix signal to noise
	var_explained = theta%*%cov_beta%*%theta
	var_epsilon = var_explained * (1-R2)/R2
	epsilon_sim=rnorm(k, mean=0, sd=sqrt(var_epsilon))
	#create the genetic variant to outcome association
	bY=alpha_vec+b_mat%*%theta+epsilon_sim
	#create names
	rs_nr = vector(mode="character", length=k)
	for(i in 1:k) rs_nr[i]=paste("rs",as.character(i),sep="")
	exposures = vector(mode="character", length=d)
	for(i in 1:d) exposures[i]=paste("E",as.character(i),sep="")
	sim_out=new("mvMRInput", betaX = b_mat, betaY = bY, snps=rs_nr, exposure=exposures, outcome = "0")
	return(sim_out)
 
}







####
# simulate mv mr data outputs object of class mvMRInput
# parameters
###
# beta_data: true risk factors (defining the number of genetic variants, correlation structure and beta distribution)
# theta: true causal association
# R2: proportion of variance of the outcome explained by the risk factors
# alpha=list(mean=0,sd=0): pleiotropy term



generate_mv_input2 = function(beta_data,theta,R2,scale=TRUE, alpha=list(mean=0,sd=0), corr_thresh=0.8){
 
	d=length(theta) # number of variables
	k=nrow(beta_data) # number of variants
	select = sample(1:ncol(beta_data),d, replace = FALSE)
	beta = scale(as.matrix(beta_data[,select]))
	#check that dimension of theta, mu_beta and Sigma_beta matches
	while(sum(sm2vec(cor(beta))>corr_thresh, na.rm=TRUE)>0){
		select = sample(1:ncol(beta_data),d, replace = FALSE)
		beta = scale(as.matrix(beta_data[,select]))		
	}
	if(scale==TRUE) beta=scale(beta)
	else beta = beta
	cov_beta = cov(beta)
	cor_beta = cor(beta)
	#pleiotropy intercept
	if(alpha$mean==0 & alpha$sd==0) alpha_vec = rep(0,k)
	else alpha_vec = rnorm(k, mean=alpha$mean, sd = alpha$sd)
	#fix signal to noise
	var_explained = theta%*%cov_beta%*%theta
	var_epsilon = var_explained * (1-R2)/R2
	epsilon_sim=rnorm(k, mean=0, sd=sqrt(var_epsilon))
	#create the genetic variant to outcome association
	bY=alpha_vec+beta%*%theta+epsilon_sim
	#create names
	rs_nr = vector(mode="character", length=k)
	for(i in 1:k) rs_nr[i]=paste("rs",as.character(i),sep="")
	exposures = vector(mode="character", length=d)
	for(i in 1:d) exposures[i]=paste("E",as.character(i),sep="")
	sim_out=new("mvMRInput", betaX = beta, betaY = bY, snps=rs_nr, exposure=exposures, outcome = "0")
	return(sim_out)
 
}





# small function to check proportion of variance explained



computeR2 = function(theta,Sigma_beta, var_epsilon){
 
	var_explained = theta %*% Sigma_beta %*%theta
	R2 = var_explained / (var_explained + var_epsilon)
	return(R2)
 
}




###
#variance decomposition
##
#R2 = var_explained / (var_explained + var_error)
#R2 var_explained +  R2 var_error = var_explained
#R2 var_error = 1-R2 var_explained
#var_error = var_explained (1-R2)/R2

# var explained = beta Sigma beta



