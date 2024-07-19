# applying the DK PACNI algorithm to HCP test retest data
# ethan whitman
# 9/25/2023

library(dplyr)
library(psych)

# load in HCP data
load('/Users/ew198/Documents/brainpace/data/hcp/inhouse_organized/hcp.Rdata')

# load in model weights for DK+GWR version
load('/Users/ew198/Documents/brainpace/results/full_model_train/dk_gwr/dk_gwr_model_Wed_Sep_20_13_44_09_2023_533.Rdata')
dk_gwr_model <- enet
coefs_dk_gwr <- coef(dk_gwr_model$finalModel, dk_gwr_model$bestTune$lambda)
rownames(coefs_dk_gwr) <- gsub("^(Right|Left)\\.(.*)$", "\\2_\\1", rownames(coefs_dk_gwr))

coefs_dk_gwr <- coefs_dk_gwr[order(rownames(coefs_dk_gwr)),]

coefs_dk_gwr <- coefs_dk_gwr[!grepl("WM.hypointensities_", names(coefs_dk_gwr))]

hcp$pacni <- rep(NA, nrow(hcp))
x<-1
for (s in 1:nrow(hcp)){
  hcp$pacni[x] <- sum(as.numeric(hcp[s,])[1:c(ncol(hcp)-5)]*as.numeric(coefs_dk_gwr)[2:length(coefs_dk_gwr)], coefs_dk_gwr[1], na.rm=T)
  x<-x+1
}


# load in HCP data
load('/Users/ew198/Documents/brainpace/data/hcp/inhouse_organized/hcp_dk.Rdata')

# load in model weights for DK version
load('/Users/ew198/Documents/brainpace/results/full_model_train/dk/model_Tue_Sep_26_10_11_19_2023_885.Rdata')
dk_model <- enet
coefs_dk <- coef(dk_model$finalModel, dk_model$bestTune$lambda)
rownames(coefs_dk) <- gsub("^(Right|Left)\\.(.*)$", "\\2_\\1", rownames(coefs_dk))

coefs_dk <- coefs_dk[order(rownames(coefs_dk)),]

coefs_dk <- coefs_dk[!grepl("WM.hypointensities_", names(coefs_dk))]

hcp_dk$pacni <- rep(NA, nrow(hcp_dk))
x<-1
for (s in 1:nrow(hcp_dk)){
  hcp_dk$pacni[x] <- sum(as.numeric(hcp_dk[s,])[1:c(ncol(hcp_dk)-5)]*as.numeric(coefs_dk)[2:length(coefs_dk)], coefs_dk[1], na.rm=T)
  x<-x+1
}


####################################################
############## test retest reliability #############
####################################################


# GWR version
ICC(data.frame(time1 = hcp$pacni[1:45], time2 = hcp$pacni[46:90]))
cor.test(hcp$pacni[1:45], hcp$pacni[46:90])

# no GWR version
ICC(data.frame(time1 = hcp_dk$pacni[1:45], time2 = hcp_dk$pacni[46:90]))
cor.test(hcp_dk$pacni[1:45], hcp_dk$pacni[46:90])


