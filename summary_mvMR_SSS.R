#
# 15th April 2019
# risk factor selection for multivariable MR based on Bayes Factors
# computation is based on summary data which speeds up computation 
#
# functions
# summarymvMR_SSS: 		main function
# summary_stochastic_search: 	shotgun stochastic search
# create_environment: 		create new environment to search
# log_binom_gamma:		prior for model size
#

# Class
# mvMRInput:			input format for summarymvMR_SSS
# MR_SSS:			output format of summarymvMR_SSS


# now synced with manuscript and BF derivation
# input: object of class mvMRInput
# NOTE hash:: prevents clashes wit MendelianRandomization R package


library(combinat)
library(hash)







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


setClass("mvMR_SSS",
         representation(Exposure = "character",
                        Outcome = "character",
                        BMAve_Estimate = "numeric",
			BestModel_Estimate = "numeric",
			BestModel = "character",
			tupel = "character",
			pp="numeric",
			pp_marginal="numeric",
			betaX="matrix",
			betaY="matrix",
			kmin="numeric", kmax="numeric", max_iter="numeric", sigma="numeric", prior_prob="numeric"
			)
)





#
# call the stochastic search function and prepare output
#



summarymvMR_SSS = function(object, kmin=1, kmax=20, max_iter=1000, sigma=0.5, prior_prob=0.5, print=FALSE){


	bX = object@betaX
	bY = object@betaY


 
	sigma_vec=rep(sigma, ncol(bX))

	sss=summary_stochastic_search(y=bY, x=bX, sigma_vec=sigma_vec, prior_prob=prior_prob, kmin=kmin, kmax=kmax, max_iter=max_iter, print=print)

	tupel_all = keys(sss$hashlogBF)
	log10BF = hash::values(sss$hashlogBF)
	log10prior = hash::values(sss$hashlogprior)
	theta_all = hash::values(sss$hashTheta)

	if(is.list(theta_all)){
		log10BF = unlist(log10BF)
		theta_all=t(matrix(unlist(theta_all), ncol = ncol(bX), byrow = TRUE))
		}

	#log_evidence = log10BF + log10prior
	log_evidence = unlist(log10BF) + unlist(log10prior)
	#best model
	best_model = tupel_all[which.max(log_evidence)]
	best_model_est = theta_all[,which.max(log_evidence)]
	#posterior for each model
	max_evidence = max(log_evidence) 
	if(max_evidence>308){
		rescale_diff = max_evidence - 308
		log_evidence = log_evidence - (rescale_diff+1)
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

	return(new("mvMR_SSS", 
		Exposure = object@exposure,
		Outcome = object@outcome,
		BMAve_Estimate = bma_coef,
		BestModel_Estimate = best_model_coef,
		BestModel = best_model,
		tupel = tupel_all,
		pp=pp, 
		pp_marginal = pp_marginal,
		betaX = bX,
		betaY= bY,
		kmin=kmin, 
		kmax=kmax, 
		max_iter=max_iter, 
		sigma=sigma, 
		prior_prob=prior_prob
	))

}










#
# core stochastic search function
#



summary_stochastic_search = function(y,x,sigma_vec, prior_prob=0.5, kmin=1, kmax=20, max_iter=1000, print=FALSE){

	hashlogprior=hash()
	hashlogBF=hash()
	hashTheta=hash()
	n_SNPs = nrow(x)
	n_RF = ncol(x)
	XtX = t(x) %*% x		# dim  d=10  x d=10	 	
 	XtY = t(x) %*% y		# dim  d=10  x 1
	YtY = t(y) %*% y

	#deterministic part
	for (i in 1:kmin){
		#create all combinations of i variables
		comb = matrix(combn(1:ncol(x), i),nrow=i)	
		new_neighbourhood = lapply(seq_len(ncol(comb)), function(i) comb[,i])
		#compute log prior
		configlogprior = lapply(new_neighbourhood, FUN = log_binom_gamma, n_RF, prior_prob)
		#compute logBF / theta
		configlogBF=lapply(new_neighbourhood, FUN = logBF_summary, XtY=XtY,XtX=XtX,YtY=YtY, sigma_vec=sigma_vec, n=n_SNPs)
		configTheta=lapply(new_neighbourhood, FUN = beta_summary, XtY=XtY,XtX=XtX,sigma_vec=sigma_vec)
		#compute the log evidence
		#save as hash table with key=tupel - value (logBF or theta)
		if(i>1){tupel=apply(comb, MARGIN=2, FUN=paste, collapse=",") }
		else{tupel=new_neighbourhood}
		
		if(i==ncol(x)){
			hashlogprior[tupel] = configlogprior[[1]]
			hashlogBF[tupel] = configlogBF[[1]]
			hashTheta[tupel] = configTheta[[1]]
		}
		else{
			hashlogprior[tupel] = configlogprior
			hashlogBF[tupel] = configlogBF
			hashTheta[tupel] = configTheta
		}
	}



	#stochastic part
	iter = 1
	#only run stochastic search if kmax > kmin
	if(kmax > kmin){
		while (iter< max_iter){
			#print(paste("sss run: ", iter))
			#draw a random configuration with weights according to BF from the last set of tupel 
			log_evidence = unlist(configlogBF) + unlist(configlogprior)
			max_evidence = max(log_evidence) 
			if(max_evidence>308){
				rescale_diff = max_evidence - 308
				log_evidence = log_evidence - (rescale_diff+1)
				max_evidence = 308
			}
			sum_calib = sum(10^(log_evidence-max_evidence)) * 10^max_evidence
			evidence = 10^(log_evidence)
			random_nr=sample(1:length(evidence), size=1, prob = evidence/sum_calib)
			random_config=new_neighbourhood[[random_nr]]
			if(print==TRUE){print(random_config)}
			#create new neighbourhood
			new_neighbourhood = create_environment(random_config,ncol(x), kmin=kmin, kmax = kmax)
			tupel=lapply(new_neighbourhood, paste, collapse=",") 
			#check which ones we visited already and query their logBF
			visited_already=has.key(tupel, hashlogBF)
			tupel_visited2=tupel[visited_already] 
			#we need the correct keys - values order 
			tupel_visited=keys(hashlogBF[tupel_visited2])
			logBF_visited=hash::values(hashlogBF[tupel_visited2])
			logprior_visited=hash::values(hashlogprior[tupel_visited2])
			if(sum(visited_already)==length(visited_already)){
				configlogBF=logBF_visited
				configlogprior=logprior_visited
				test2=unlist(lapply(tupel_visited, FUN=strsplit, split=","), recursive=FALSE)
				new_neighbourhood = lapply(test2,FUN=helper)		
				}
			else{
				#compute log prior for the neighboorhood				
				configlogprior_new = lapply(new_neighbourhood[!visited_already], FUN = log_binom_gamma, n_RF, prior_prob)			
				#compute logBF / theta for the neighboorhood
				configlogBF_new=lapply(new_neighbourhood[!visited_already],  FUN = logBF_summary, XtY=XtY,XtX=XtX,YtY=YtY, sigma_vec=sigma_vec, n=n_SNPs)
				configTheta_new=lapply(new_neighbourhood[!visited_already],   FUN = beta_summary, XtY=XtY,XtX=XtX,sigma_vec=sigma_vec)
				#save as hash table with key=tupel - value (logBF or theta)
				tupel_insert=lapply(new_neighbourhood[!visited_already], paste, collapse=",")
				hashlogprior[tupel_insert] = configlogprior_new 
				hashlogBF[tupel_insert] = configlogBF_new
				hashTheta[tupel_insert] = configTheta_new
				#prepare new random draw
				#tupel=merge tupel_visited and tupel_insert 
				c_tupel = c(tupel_visited, tupel_insert)
				test2=unlist(lapply(c_tupel, FUN=strsplit, split=","), recursive=FALSE)
				new_neighbourhood = lapply(test2,FUN=helper)
				#configlogBF = merge logBF_visited and configlogBF_new
				configlogBF = c(logBF_visited, configlogBF_new) 
				configlogprior = c(logprior_visited, configlogprior_new)
				#and move on with the new neighbourhood and logBF as starting point for the random draw
				}
			iter = iter +1 
			}
		}

	return(list(hashlogBF=hashlogBF, hashTheta=hashTheta, hashlogprior=hashlogprior))

}













#
# function that creates a new envirnoment given a starting configuraiton current config, 
# the size of the design matrix m, 
# and the minimum (kmin=1) and maximum size (kmax=20) of models considered
#



create_environment=function(current_config, m, kmin=1, kmax=20){

	actual_kmax=min(kmax,m)
	current_size=length(current_config)
	all_members = 1:m
	possible_new_members = setdiff(all_members,current_config)

	new_neighbourhood = list()

	if(current_size == kmin){

		#swap move new_size = current_size
		for(i in 1:current_size){
			j=length(new_neighbourhood)
			loop_config=current_config[-i]
			for(k in 1:length(possible_new_members)){
				new_neighbourhood[[j+k]] = sort(c(loop_config, possible_new_members[k]))
			}
		}
		
		#add move new_size = current_size + 1
		j=length(new_neighbourhood)
		for(i in 1:length(possible_new_members)){
			new_neighbourhood[[j+i]] = sort(c(current_config, possible_new_members[i]))
		}

		j=length(new_neighbourhood)
		new_neighbourhood[[j+1]] = current_config	

	}

	else if(current_size == kmax){

		#delete move new_size = current_size - 1
		for(i in 1:current_size){
			new_neighbourhood[[i]] = current_config[-i]
		}

		#swap move new_size = current_size
		for(i in 1:current_size){
			j=length(new_neighbourhood)
			loop_config=current_config[-i]
			for(k in 1:length(possible_new_members)){
				new_neighbourhood[[j+k]] = sort(c(loop_config, possible_new_members[k]))
			}
		}

		j=length(new_neighbourhood)
		new_neighbourhood[[j+1]] = current_config	

	}

	else{
		#delete move new_size = current_size - 1
		for(i in 1:current_size){
			new_neighbourhood[[i]] = current_config[-i]
		}

		#swap move new_size = current_size
		for(i in 1:current_size){
			j=length(new_neighbourhood)
			loop_config=current_config[-i]
			for(k in 1:length(possible_new_members)){
				new_neighbourhood[[j+k]] = sort(c(loop_config, possible_new_members[k]))
			}
		}

		#add move new_size = current_size + 1
		j=length(new_neighbourhood)
		for(i in 1:length(possible_new_members)){
			new_neighbourhood[[j+i]] = sort(c(current_config, possible_new_members[i]))
		}
		
		j=length(new_neighbourhood)
		new_neighbourhood[[j+1]] = current_config		


	}

	return(new_neighbourhood)

}




#
# function to compute a binomial prior with prior prob
#


log_binom_gamma = function(gamma, n_x, prior_prob){

	n_model_size = length(gamma)	
	log_prior = log10(prior_prob^n_model_size * (1-prior_prob)^(n_x-n_model_size))

return(log_prior)

}







#
# make table: report the best individual models
#


#input
#BMA_output object of class mvMR_SSS
#prior_sigma needs to be specified as in main function
#top = 10 how many top models to be reported?
#write.out = FALSE write out table as csv file
#csv.file.name = "best_model_out"



sss.report.best.model = function(BMA_output, prior_sigma=0.5, top = 10, digits = 3, write.out = TRUE, csv.file.name="best_model_out"){


	if(class(BMA_output)[1] !="mvMR_SSS"){
		message("Input needs to be of class mvMR_SSS")
		break
		}

	pp=BMA_output@pp
	models=BMA_output@tupel
	betaX=BMA_output@betaX
	betaY=BMA_output@betaY
	rf=BMA_output@Exposure
	prior_sigma=BMA_output@sigma
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
#BMA_output object of class mvMR_SSS
#top = 10 how many top models to be reported?
#write.out = FALSE write out table as csv file
#csv.file.name = "mr_bma_out"



sss.report.mr.bma = function(BMA_output,  top = 10, digits = 3, write.out = TRUE, csv.file.name="mr_bma_out"){

	if(class(BMA_output)[1] !="mvMR_SSS"){
		message("Input needs to be of class mvMR_SSS")
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


