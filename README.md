# DNM_ANM_PGS
This repository contains analyses we did to assess the association between phased de novo SNV mutations and a polygenic score for age at natural menopause (ANM). 

Phased de novo calls were from Joanna Kaplani who produce them for her 2022 paper: https://doi.org/10.1038/s41586-022-04712-2. Read-based approach was used to phase a de novo SNV, using heterozygous variants within 500 bp of the DNM. Phasing of mutations was performed with a custom Python (3) script available at GitHub (https://github.com/queenjobo/PhaseMyDeNovo).

We performed possion regression, correcting for the same QC covariates following Kaplani et al 2022 paper:
maternally phased de novo SNVs ~ maternal PGS + maternal age + 
                 proband mean coverage + mother mean coverage + 
                 proband prop aligned reads + mother prop aligned reads + 
                 proband snvs + mother snvs + median VAF + median bayes

We used unrelated individuals with European ancestry defined by GEL. For PGS, we regressed out 20 PCs and scaled the residuals. Higher PGS indicates later ANM.

Sanity check: the paternal model, regressing paternal counts against maternal PGS, and regressing maternal counts against paternal PGS.

For rare variant burden score analysis, we calculated a weighted burden score considering all ANM-associated genes using effect sizes estimated in the exome-wide association analysis. We weighted missense variants and PTVs by their corresponding masks. Higher weighted burden score indicates later ANM. In this analysis, we replaced the PGS with the burden score in the above models and ran poisson regression. We also fitted each ANM gene separately where carriers are coded as 1 and non-carriers are coded as 0. 
