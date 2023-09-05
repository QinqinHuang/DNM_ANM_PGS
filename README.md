# DNM_ANM_PGS
This repository contains analyses we did to assess the association between phased de novo SNV mutations and a polygenic score for age at natural menopause. 

Phased de novo calls were from Joanna Kaplani who produce them for her 2022 paper: https://doi.org/10.1038/s41586-022-04712-2. Read-based approach was used to phase a de novo SNV, using heterozygous variants within 500 bp of the DNM. Phasing of mutations was performed with a custom Python (3) script available at GitHub (https://github.com/queenjobo/PhaseMyDeNovo).

We performed possion regression, correcting for the same QC covariates following Kaplani et al 2022 paper:
maternally phased de novo SNVs ~ maternal PRS + maternal age + 
                 proband mean coverage + mother mean coverage + 
                 proband prop aligned reads + mother prop aligned reads + 
                 proband snvs + mother snvs + median VAF + median bayes

Sanity check: the paternal model, regressing paternal counts against maternal PGS, and regressing maternal counts against paternal PGS.

