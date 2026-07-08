# applying the DK PACNI algorithm to ADNI ADSP data
# ethan whitman
# 1/10/24

library(plm)
library(lmtest)
library(dplyr)
library(reshape2)
library(grid)
library(ggplot2)
library(haven)
library(labelled)
library(kableExtra)
library(sjPlot)
library(Hmisc)
library(broom)
library(performance)
library(arm)
library(faraway)
library(pscl)
library(nnet)
library(brant)

# load in ADNI data
load('/Users/ew198/Documents/brainpace/data/adni/inhouse_organized/adni_adsp.Rdata')

# load in model weights for DK+GWR version
load('/Users/ew198/Documents/brainpace/results/full_model_train/dk_gwr/dk_gwr_model_Wed_Sep_20_13_44_09_2023_533.Rdata')
dk_gwr_model <- enet
coefs_dk_gwr <- coef(dk_gwr_model$finalModel, dk_gwr_model$bestTune$lambda)
rownames(coefs_dk_gwr) <- gsub("^(Right|Left)\\.(.*)$", "\\2_\\1", rownames(coefs_dk_gwr))

coefs_dk_gwr <- coefs_dk_gwr[order(rownames(coefs_dk_gwr)),]

coefs_dk_gwr <- coefs_dk_gwr[!grepl("WM.hypointensities_", names(coefs_dk_gwr))]

data.frame(coefs=c(names(coefs_dk_gwr)[2:length(coefs_dk_gwr)], NA, NA, NA, NA, NA, NA),
           adni=colnames(adni_adsp)[1:(ncol(adni_adsp))])

adni_adsp$pacni <- rep(NA, nrow(adni_adsp))
x<-1
for (s in 1:nrow(adni_adsp)){
  adni_adsp$pacni[x] <- sum(as.numeric(adni_adsp[s,])[1:c(ncol(adni_adsp)-7)]*as.numeric(coefs_dk_gwr)[2:length(coefs_dk_gwr)], coefs_dk_gwr[1], na.rm=T)
  x<-x+1
}


adni_adsp_merge_dx$age_at_scan <- adni_adsp_merge_dx$age_at_scan / 12

adni_adsp <- merge(adni_adsp, 
                             adni_adsp_merge_dx[,c('RID', 'date_scanned', 'age_at_scan', 'baseline_date', 'age_at_baseline', 'DX', 'dx_recode', 'EXAMDATE', 'time_diff', 'APOE4')], 
                             by = c('RID', 'age_at_scan'))

adni_adsp_pacni <- adni_adsp

#save(adni_adsp_pacni, file='/Users/ew198/Documents/brainpace/data/adni/inhouse_organized/adni_adsp_pacni.Rdata')

