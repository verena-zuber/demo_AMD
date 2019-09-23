# Selecting causal risk factors from high-throughput experiments using multivariable Mendelian randomization

(20th of July 2018, updated 23rd of September 2019)

Documentation of the multivariable MR analysis to prioritise NMR metabolites as risk factors for age-related macular degeneration (AMD) using MR-BMA as detailed in Selecting causal risk factors from high-throughput experiments using multivariable Mendelian randomization by Verena Zuber, Johanna Maria Colijn, Caroline Klaver, and Stephen Burgess.


The summary-level data on genetic association on NMR metabolites and AMD used for our study is now part of the data challenge of the MR conference 2019 in Bristol

https://www.mendelianrandomization.org.uk/the-mr-data-challenge-2019/

and available in the R-package MRChallenge2019

https://github.com/WSpiller/MRChallenge2019


Alternatively, use:

load("amd_example")


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





Compilation and replication of the application example:
------------------------------------------------------

1. all genetic instrumental variables included (n=148, runtime: 14 minutes):

knit("run-BMA-AMD-n148.Rmd")

markdownToHTML('run-BMA-AMD-n148.md', 'run-BMA-AMD-n148.html', options=c("use_xhml"))

2. after excluding outliers and influential points (n=145, runtime: 11 minutes):

knit("run-BMA-AMD-n145.Rmd")

markdownToHTML('run-BMA-AMD-n145.md', 'run-BMA-AMD-n145.html', options=c("use_xhml"))




Package dependencies:
---------------------

require(MRChallenge2019)

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