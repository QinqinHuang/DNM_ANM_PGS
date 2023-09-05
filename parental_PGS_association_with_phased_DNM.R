#------------------------------------------------------------
# 2022-03-04
#
# Run the following regression models:
# (1) Paternally phased DNMs ~ father_PRS + cov
#     Maternally phased DNMs ~ mother_PRS + cov
# (2) Sanity check: Paternally phased DNMs ~ mother_PRS + cov
#                   Maternally phased DNMs ~ father_PRS + cov
#
# cov: using all that Joanna used
#
# PRS - using PC adjusted PRS
# 
#----
# Regression model:
# phased DMNs: glm.nb breaks when trying the full model,
#   so using possion regression
#------------------------------------------------------------
library(data.table)
library(foreach)
library(MASS)

setwd("/re_gecip/paediatrics/QQHuang/JohnPerry/Results_2022Mar/PRS_association_unrelated/")


## read phenotype data
DNM = fread("../pheno_table/DNM_unrelateds_PRS_2022_03_04.tab")

# keep unrelated probands
DNM = DNM[related_proband == "unrelated"]

# v14 sample list
v14full = fread("/re_gecip/paediatrics/QQHuang/Sample_stats/v14/v14_sampleinfo.txt")

# PCs calculated within EUR (N=61946)
PCs_EUR = fread("/re_gecip/paediatrics/QQHuang/Sample_stats/v14/PCs_20_withinEUR.txt")


## PRS
PRS_mother = fread("../pheno_table/PRS_raw_mother_all.txt")
PRS_father = fread("../pheno_table/PRS_raw_father_all.txt")




#----- define models ------
# model: age + all QC metrics
model_fa = formula("fphased ~ PRS_fa + father_age + 
                 proband_mean_coverage + father_mean_coverage +
                 proband_aligned_reads + father_aligned_reads +
                 proband_snvs + father_snvs + medVAF + medbayes")
model_mo = formula("mphased ~ PRS_mo + mother_age + 
                 proband_mean_coverage + mother_mean_coverage + 
                 proband_aligned_reads + mother_aligned_reads + 
                 proband_snvs + mother_snvs + medVAF + medbayes")

# Sanity check
model_fDNM_mPRS = formula("fphased ~ PRS_mo + father_age + 
                 proband_mean_coverage + father_mean_coverage +
                 proband_aligned_reads + father_aligned_reads +
                 proband_snvs + father_snvs + medVAF + medbayes")
model_mDNM_fPRS = formula("mphased ~ PRS_fa + mother_age + 
                 proband_mean_coverage + mother_mean_coverage + 
                 proband_aligned_reads + mother_aligned_reads + 
                 proband_snvs + mother_snvs + medVAF + medbayes")

# model to adjust PCs for PRS
model_PRS_PC = formula("rawPRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")



##--- My function to get summary statistics for four models ---
getstatistics = function(ddmo, ddfa) {
  
  # father PRS ~ fDNM
  fit_fa = glm(model_fa, data = ddfa, family = poisson())
  # mother PRS ~ mDNM
  fit_mo = glm(model_mo, data = ddmo, family = poisson())
  # Sanity check: father PRS ~ mDNM
  fit_mDNM_fPRS = glm(model_mDNM_fPRS, data = ddfa, family = poisson())
  # Sanity check: mother PRS ~ fDNM
  fit_fDNM_mPRS = glm(model_fDNM_mPRS, data = ddmo, family = poisson())
  
  # summary stats
  gettable_fa = function(myfit, label) {
    d = as.data.frame(summary(myfit)$coefficients)["PRS_fa", ]
    names(d)[3] = gsub(" ", "", names(d)[3])
    names(d)[-3] = c("beta", "se", "Pval")
    d$model = label
    d$variable = c("PRS_fa")
    d$N = length(myfit$residuals)
    return(d)
  }
  gettable_mo = function(myfit, label) {
    d = as.data.frame(summary(myfit)$coefficients)["PRS_mo", ]
    names(d)[3] = gsub(" ", "", names(d)[3])
    names(d)[-3] = c("beta", "se", "Pval")
    d$model = label
    d$variable = c("PRS_mo")
    d$N = length(myfit$residuals)
    return(d)
  }
  
  returndd = rbind(gettable_fa(fit_fa, "father_model"), 
                   gettable_fa(fit_mDNM_fPRS, "mDNM_fPRS"), 
                   gettable_mo(fit_mo, "mother_model"), 
                   gettable_mo(fit_fDNM_mPRS, "fDNM_mPRS") )
  
  return(returndd)
}


PRSlist = c("ANM_EUR", "ANM_excl7_EUR")
all(PRSlist %in% names(PRS_mother))


## run the analysis
# normalise PRS column name
ddPRSmo = PRS_mother[, .(mother_id, ID, superpop)]
ddPRSmo$rawPRS = PRS_mother[, get(ii)]
ddPRSmo = ddPRSmo[!is.na(rawPRS)]

ddPRSfa = PRS_father[, .(father_id, ID, superpop)]
ddPRSfa$rawPRS = PRS_father[, get(ii)]
ddPRSfa = ddPRSfa[!is.na(rawPRS)]
  
# using within EUR PCs to adjust PRS
ddPRSmo_PCs = merge(ddPRSmo, PCs_EUR, by = "ID")
lm_PRS = lm(model_PRS_PC, data = ddPRSmo_PCs)
ddPRSmo_PCs$PRS_corr = as.numeric(scale(lm_PRS$residuals))

ddPRSfa_PCs = merge(ddPRSfa, PCs_EUR, by = "ID")
lm_PRS = lm(model_PRS_PC, data = ddPRSfa_PCs)
ddPRSfa_PCs$PRS_corr = as.numeric(scale(lm_PRS$residuals))

# add adjusted PRS into the phenotype table
# mothers: keep only unrelated
ddmo = merge(DNM[related_mother == "unrelated"], 
             ddPRSmo_PCs[,.(mother_id, PRS_mo = PRS_corr)], by = "mother_id")
# fathers: keep only unrelated
ddfa = merge(DNM[related_father == "unrelated"], 
             ddPRSfa_PCs[,.(father_id, PRS_fa = PRS_corr)], by = "father_id")

# results
PRSstats = getstatistics(ddmo, ddfa)


fwrite(PRSstats, "PhasedDNM_faPRS_or_moPRS_20PCadj_poisson_stats_unrelated.txt", sep = "\t")







