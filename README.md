# demo_AMD
# 20th of July 2018

Documentation of the multivariable MR analysis to prioritise NMR metabolites as risk factors for AMD using MR-BMA as detailed in Selecting causal risk factors from high-throughput experiments using multivariable Mendelian randomization by Verena Zuber, Johanna Maria Colijn, Caroline Klaver, and Stephen Burgess.

Example data is based on publicely available summary data on: 

AMD
---
http://amdgenetics.org/

Fritsche, L. G., Igl, W., Bailey, J. N. C., Grassmann, F., Sengupta, S., Bragg-Gresham, J. L., Burdon, K. P., Hebbring, S. J., Wen, C., Gorski, M., et al. (2015). A large genome-wide association study of age-related macular degeneration highlights contributions of rare and common variants. Nature Genetics 48, 134 

NMR metabolites
---------------
http://computationalmedicine.fi/data#NMR_GWAS

Kettunen, J., Demirkan, A., Wuertz, P., Draisma, H. H. M., Haller, T., Rawal, R., Vaarhorst, A., Kangas, A. J., Lyytikaeinen, L.-P., Piri- nen, M., et al. (2016). Genome-wide study for circulating metabolites identifies 62 loci and reveals novel systemic effects of lpa. Nature Communications 7.



Compilation:

1. all genetic instrumental variables included:
knit("run-BMA-AMD-n148.Rmd")
markdownToHTML('run-BMA-AMD-n148.md', 'run-BMA-AMD-n148.html', options=c("use_xhml"))

2. after excluding outliers and influential points:
knit("run-BMA-AMD-n145.Rmd")
markdownToHTML('run-BMA-AMD-n145.md', 'run-BMA-AMD-n145.html', options=c("use_xhml"))





Package dependencies:

require(knitr)
require(markdown)
require(ggplot2)
require(combinat)
require(hash)
require(corpcor)
