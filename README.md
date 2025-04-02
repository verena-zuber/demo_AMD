# Selecting causal risk factors from high-throughput experiments using multivariable Mendelian randomization

(20th of July 2018, updated 14th of August 2021)


MR-BMA:
-------------
MR-BMA is a Bayesian algorithm to perform risk factor selection in multivariable MR. 
The main function to run MR-BMA [Zuber et al (2020)] is summarymvMR_SSS() which needs an mvMRInput object:

amd_nmr_input=new("mvMRInput", betaX = as.matrix(betaX_ivw), betaY = as.matrix(amd_beta_ivw), snps=rs, exposure=rf, outcome = "amd")

BMA_output=summarymvMR_SSS(amd_nmr_input,kmin=1,kmax=12, prior_prob=0.1, max_iter=100)

Please take care: It is necessary to adjust the following parameters according to your data and how you want to run the search algorithm (exhaustive search of models (all combinations of risk factors) or stochastic search):

- kmin: The minimum model size considered.
- kmax: The maximum model size considered. Please note computing all possible combinations of more than 12 risk factors may become computationally infeasible. If kmin = kmax an exhaustive search is performed, if kmin<kmax a stochastic search is performed. Stochastic search is recommended for more than 12-15 risk factors.

! stochastic search may be very time consuming where realisitic and feasible we recommend the exhaustive search by settin kmin and kmax to the same value. Importantly even if there are more risk factors in the model the maximum model size may be much smaller.

- prior_prob: Please adjust according to prior knowledge. the smaller the prior_prob, the sparser the models.
- max_iter: Number of stochastic search runs. If stochastic search is used, it is recommended to set initially to a small number to make sure the algorithm is running (eg 100 as specified here). NOTE: for the final output we recommend at least 10k, ideally 100k runs.

The output of the MR-BMA function is complex and includes all combination of risk factors visited. Additionally, this packages includes two functions to summarize the output and create output tables which can be saved as text files, one for the best model (combination of risk factors) and one for the model-averaged results :

best.model.out = sss.report.best.model(BMA_output, top = 10, write.out = TRUE, csv.file.name="amd_best_10models_n145")

best.model.out

mr.bma.out = sss.report.mr.bma(BMA_output, top = 10, write.out = TRUE, csv.file.name="amd_mr_bma_n145")

mr.bma.out


The outout of the sss.report.mr.bma can be further analyzed for diagnostic tests to identify influential and outlying IVs and calculate empirical p-values.

First, the diagnostic function detects influential genetic variants and outliers for all models with the largest posterior probabilities (larger than diag_ppthresh = 0.02). 

diagnostics_output = diagnostics(BMA_output, diag_ppthresh = 0.02)

diagnostics_output$rmCD

diagnostics_output$rmO

Finally, we have proposed a permutation procedure [Levin et al (2021)] to calculate empirical p-values which is a computationally very intensive (may be hours of runtime). Please make sure to run first with few runs (nrepeat = 10 or 100) and evaluate the run time. Consider running this command on a remote server and saving the permutation p-values as provided by the function create.permutations(), if save.matrix=TRUE the permuted posterior probabilities are saved for all repetitions. The function calculate.p calculates the actual p-values, here they are saved in the object empirical.p. 
NOTE: For stable results we recommend nrepeat = 100000.

permute_bma = create.permutations(BMA_output, nrepeat = 100, save.matrix=TRUE, file.name = "permutation_mrBMA")

empirical.p = calculate.p(BMA_output, permute_bma)



Application example: NMR metabolites as risk factors for age-related macular degeneration
--------------------------
This repository includes the documentation of the multivariable MR analysis to prioritise NMR metabolites as risk factors for age-related macular degeneration (AMD) using MR-BMA as detailed in the manuscript Selecting causal risk factors from high-throughput experiments using multivariable Mendelian randomization by Verena Zuber, Johanna Maria Colijn, Caroline Klaver, and Stephen Burgess [Zuber et al 2020].
This package provides the summary-level data on genetic association on NMR metabolites and AMD as example data in the Rdata file "amd_example".  Alternatively, the summary-level data on genetic association on NMR metabolites and AMD used for our study is now part of the data challenge of the MR conference 2019 in Bristol
https://www.mendelianrandomization.org.uk/the-mr-data-challenge-2019/
and available in the R-package MRChallenge2019
https://github.com/WSpiller/MRChallenge2019


Compilation and replication of the application example:
------------------------------------------------------

1. all genetic instrumental variables included (n=148, runtime: 14 minutes):

knit("run-BMA-AMD-n148.Rmd")

markdownToHTML('run-BMA-AMD-n148.md', 'run-BMA-AMD-n148.html', options=c("use_xhml"))

2. after excluding outliers and influential points (n=145, runtime: 11 minutes):

knit("run-BMA-AMD-n145.Rmd")

markdownToHTML('run-BMA-AMD-n145.md', 'run-BMA-AMD-n145.html', options=c("use_xhml"))




References:
-----------
Zuber, V., Colijn, J. M., Klaver, C. & Burgess, S. Selecting likely causal risk factors from high-throughput experiments using multivari- able mendelian randomization. Nature Communications 11, 29 (2020). URL https://doi.org/10.1038/s41467-019-13870-3.

Levin, M. G. et al. Prioritizing the role of major lipoproteins and subfractions as risk factors for peripheral artery disease. Circulation (2021). URL https://doi.org/10.1161/CIRCULATIONAHA.121. 053797.



Data availability:
------------------
We thank the authors and contributors of the following studies for publishing genome-wide summary data:


AMD
---
http://amdgenetics.org/

Fritsche, L. G., Igl, W., Bailey, J. N. C., Grassmann, F., Sengupta, S., Bragg-Gresham, J. L., Burdon, K. P., Hebbring, S. J., Wen, C., Gorski, M., et al. (2015). A large genome-wide association study of age-related macular degeneration highlights contributions of rare and common variants. Nature Genetics 48, 134 

NMR metabolites
---------------
http://computationalmedicine.fi/data#NMR_GWAS

Kettunen, J., Demirkan, A., Wuertz, P., Draisma, H. H. M., Haller, T., Rawal, R., Vaarhorst, A., Kangas, A. J., Lyytikaeinen, L.-P., Pirinen, M., et al. (2016). Genome-wide study for circulating metabolites identifies 62 loci and reveals novel systemic effects of lpa. Nature Communications 7.



R version:
----------
Developed under: R version 3.4.2 (2017-09-28) -- "Short Summer"
Tested under: R version 3.6.1 (2019-07-05) -- "Action of the Toes"



Installation:
-------------

source("summary_mvMR_SSS.R")

source("summary_mvMR_BF.R")

Installation time: < 1 sec




Package dependencies:
---------------------

require(knitr)

require(markdown)

require(ggplot2)

require(combinat)

require(hash)

require(corpcor)



Licence:
---------

The MIT Licence
https://opensource.org/licenses/MIT
