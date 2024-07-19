# running grid search using DK atlas instead of DKT
# 9/8/23
# ethan whitman

library(readr)
library(caret)
library("matrixStats")
library(tidyr)
library(plyr)
library(foreign)
library(dplyr)
library(glmnet)
library(glmnetUtils) # for more memory-efficient cv ridge regression. must be loaded after glmnet!
suppressMessages(require(optparse))

################### OPTIONS #############################

root <- "/cifs/hariri-long/"
behavvar <- "PaceOfAgingP45" 
workdir <- paste0(root, "Projects/ethan/pacni/performance_output/cluster/tuning_parameters_dk_gwr/")
# ICCthr <- 0.75
# each brainvarlist is a list of data frame names. They should have corresponding variables for ROIs within the DF, named ROIs_<DF>
brainvarlists <- list(list("asegALL", "CT", "SA", "GMV", "GWR")) 
# brainvarlists <- list(list("asegALL", "gCT", "gSA", "gGMV")) 
# brainvarlists <- list(list("asegALL", "gCT", "gSA", "gGMV", "FA","MD","AD","RD","WMH","GFC")) 

################### END OPTIONS #############################

# Read Data

## trying the model using DKT atlas for completeness and also for potentially needing this in ADNI
CT <- read.csv(paste0(root, "/Database/DBIS/Imaging/FreeSurfer/v6.0/FreeSurfer_aparc_ThickAvg.csv"), na.strings=".")
SA <- read.csv(paste0(root, "/Database/DBIS/Imaging/FreeSurfer/v6.0/FreeSurfer_aparc_SurfArea.csv"), na.strings=".")
GMV <- read.csv(paste0(root, "/Database/DBIS/Imaging/FreeSurfer/v6.0/FreeSurfer_aparc_GrayVol.csv"), na.strings=".")
aseg <- read.csv(paste0(root, "/Database/DBIS/Imaging/HCPMPP/aseg_FreeSurfer6.0.csv"), na.strings=".")
GWR_temp <- read.csv(paste0(root, "/Database/DBIS/Imaging/FreeSurfer/v6.0/FreeSurfer_w-g.pct_Mean.csv"), na.strings=".")
GWR <- GWR_temp[,1:(ncol(GWR_temp)-1)]
CT$ID <- as.numeric(sub("DMHDS","",CT$ID)); names(CT)[1] <- "snum"
SA$ID <- as.numeric(sub("DMHDS","",SA$ID)); names(SA)[1] <- "snum"
GMV$ID <- as.numeric(sub("DMHDS","",GMV$ID)); names(GMV)[1] <- "snum"
GWR$ID <- as.numeric(sub("DMHDS","",GWR$ID)); names(GWR)[1] <- "snum"


## add Glasser and WM measures for "delux" version

# gCT <- read.csv(paste0(root, "/Database/DBIS/Imaging/HCPMPP/corrCT_HCPMPP.csv"), header = TRUE, na.strings=".")
# gSA <- read.csv(paste0(root, "/Database/DBIS/Imaging/HCPMPP/SA_HCPMPP.csv"), header = TRUE, na.strings=".")
# gGMV <- read.csv(paste0(root, "/Database/DBIS/Imaging/HCPMPP/volume_HCPMPP.csv"), header = FALSE, na.strings=".")
# WMH <- read.csv(paste0(root, "/Database/DBIS/Imaging/WMHvolume_UBO_052419.csv"), na.strings=".")
# FA <- read.csv(paste0(root, "/Database/DBIS/Imaging/DTI_ENIGMA_ROIs_averageFA.csv"), na.strings=".")
# MD <- read.csv(paste0(root, "/Database/DBIS/Imaging/DTI_ENIGMA_ROIs_averageMD.csv"), na.strings=".")
# AD <- read.csv(paste0(root, "/Database/DBIS/Imaging/DTI_ENIGMA_ROIs_averageAD.csv"), na.strings=".")
# RD <- read.csv(paste0(root, "/Database/DBIS/Imaging/DTI_ENIGMA_ROIs_averageRD.csv"), na.strings=".")
### retests
# gCT_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/HCPMPP/corrCT_HCPMPP_retest.csv"), header = TRUE, na.strings=".")
# gSA_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/HCPMPP/SA_HCPMPP_retest.csv"), header = TRUE, na.strings=".")
# gGMV_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/HCPMPP/volume_HCPMPP_retest.csv"), header = FALSE, na.strings=".")
# gCT <- merge(gCT, gCT_retest, all=TRUE)
# gSA <- merge(gSA, gSA_retest, all=TRUE)
# gGMV <- merge(gGMV, gGMV_retest, all=TRUE)
################## need WMH
# FA_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/DTI_ENIGMA_ROIs_averageFA_retests.csv"), na.strings=".")
# MD_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/DTI_ENIGMA_ROIs_averageMD_retests.csv"), na.strings=".")
# AD_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/DTI_ENIGMA_ROIs_averageAD_retests.csv"), na.strings=".")
# RD_retest <- read.csv(paste0(root, "/Database/DBIS/Imaging/retests/DTI_ENIGMA_ROIs_averageRD_retests.csv"), na.strings=".")
# FA_retest$ID <- as.numeric(sub("DMHDS","",FA_retest$ID)); names(FA_retest)[1] <- "snum"
# MD_retest$bidsid <- as.numeric(sub("sub-","",MD_retest$bidsid)); names(MD_retest)[1] <- "snum"
# AD_retest$bidsid <- as.numeric(sub("sub-","",AD_retest$bidsid)); names(AD_retest)[1] <- "snum"
# RD_retest$bidsid <- as.numeric(sub("sub-","",RD_retest$bidsid)); names(RD_retest)[1] <- "snum"
# FA <- merge(FA, FA_retest, all=TRUE)
# MD <- merge(MD, MD_retest, all=TRUE)
# AD <- merge(AD, AD_retest, all=TRUE)
# RD <- merge(RD, RD_retest, all=TRUE)

## ROIs
# ## reformat DTI names
# names(FA)[9:ncol(FA)] <- paste0(names(FA)[9:ncol(FA)], "_FA")
# names(MD)[9:ncol(MD)] <- paste0(names(MD)[9:ncol(MD)], "_MD")
# names(RD)[9:ncol(RD)] <- paste0(names(RD)[9:ncol(RD)], "_RD")
# names(AD)[9:ncol(AD)] <- paste0(names(AD)[9:ncol(AD)], "_AD")
# ROIs_FA <- c(names(FA)[grepl("_L|_R|CorpCal|MidCer|Pont",names(FA))],"Fornix_FA")
# ROIs_MD <- c(names(MD)[grepl("_L|_R|CorpCal|MidCer|Pont",names(MD))],"Fornix_MD")
# ROIs_AD <- c(names(AD)[grepl("_L|_R|CorpCal|MidCer|Pont",names(AD))],"Fornix_AD")
# ROIs_RD <- c(names(RD)[grepl("_L|_R|CorpCal|MidCer|Pont",names(RD))],"Fornix_RD")
# ROIs_WMH <- c("wmhVol")
## store names
ROIs_CT <- paste0("CT_",names(CT)[!grepl("snum",names(CT))]); names(CT)[!grepl("snum",names(CT))] <- ROIs_CT # need to add "CT" so unique from others
ROIs_SA <- paste0("SA_",names(SA)[!grepl("snum",names(SA))]); names(SA)[!grepl("snum",names(SA))] <- ROIs_SA # need to add "SA" so unique from others
ROIs_GMV <- paste0("GMV_",names(GMV)[!grepl("snum",names(GMV))]); names(GMV)[!grepl("snum",names(GMV))] <- ROIs_GMV # need to add "GMV" so unique from others
# ROIs_gCT <- names(gCT)[grepl("_ROI",names(gCT))] #region_labels[1:360,"V1"]
# ROIs_gSA <- names(gSA)[grepl("_ROI",names(gSA))] #region_labels[1:360,"V1"]
ROIs_aseg <- names(aseg)[grepl("Cerebellum.Cortex|Caudate|Putamen|Pallidum|Brain.Stem|Hippocampus|Amygdala|Accumbens|Thalamus|VentralDC",names(aseg))]
asegALL <- aseg
ROIs_asegALL <- names(asegALL)[8:ncol(asegALL)] # ALL aseg values
ROIs_GWR <- paste0("GWR_",names(GWR)[!grepl("snum",names(GWR))]); names(GWR)[!grepl("snum",names(GWR))] <- ROIs_GWR
############### gGMV, NEEDS UPDATED!
# names(gGMV) <- c("id","avgVol",sub("corrCT","vol",ROIs_gCT))
# ROIs_gGMV <- names(gGMV)[grepl("_ROI",names(gGMV))] 
# gGMV$snum <- as.numeric(sub("sub-","",gGMV$id))



## behavioral data
behavdata <- read.spss(paste0(root,'/Projects/Annchen/DBIS/Gradients/GradientReliabilityPaper/data/Cortex_Cogn_Checking2023.sav'),to.data.frame = TRUE,use.value.labels=FALSE) # for POA
### join datasets to pare down to SMs who have all brain data (there is no one with FC that doesn't have struct)
behav_merged <- plyr::join_all(list(
  behavdata,
  aseg[,c("snum","HCPMPP_inclusive")]),
  by="snum",type="full") 
behav_merged <- behav_merged[which( behav_merged$HCPMPP_inclusive %in% c("1","1f") & 
                                      !is.na(behav_merged[,paste(behavvar)])), ]




# PREDICT

for (iter in 1:100){
  
  for ( brainvarlist in brainvarlists ) {
    
    # make ROI list
    braindat <- data.frame(snum=aseg$snum) # might not need to initialize this way, but just in case
    ROIs_full <- c()      
    for (i in 1:length(brainvarlist)) {
      brainvarcur <- brainvarlist[[i]]
      ROIs_cur <- get(paste0("ROIs_", brainvarcur))
      # ICCs_cur <- get(paste0("ICCs_", brainvarcur))
      # threshold ROIs by ICC
      # ROIs_full <- c(ROIs_full, ROIs_cur[!is.na(ICCs_cur) & ICCs_cur > ICCthr])
      ROIs_full <- c(ROIs_full, ROIs_cur)
      braindat <- merge(braindat, get(brainvarcur), by="snum")
      if (i>1) { brainvar <- paste(brainvar, brainvarcur, sep=".") } else { brainvar <- brainvarcur }
    } # end loop thru brainvar
    braindat <- braindat[which(braindat$HCPMPP_inclusive %in% c("1","1f")),]
    
    data <- merge(braindat[,c("snum",ROIs_full)], behav_merged[, c("snum","sex", behavvar)], by="snum")
    data <- data[complete.cases(data),]
    
    training.samples <- data[, paste(behavvar)] %>% createDataPartition(p = 0.9, list = FALSE)
    train.data <- data[training.samples, ]
    test.data  <- data[-training.samples, ]
    
    # regress sex and motion from training set
    lm <- lm(train.data[,paste(behavvar)] ~ train.data$sex)
    train.data$behav_resids <- scale(lm$residuals) #### sure we want to scale this way?
    # adjust using same parameters in test set
    coefs <- lm$coefficients # 1=intercept, 2=sex
    test_fitted_covars <- coefs[1] + coefs[2]*(as.numeric(test.data$sex)-1)
    test.data$behav_adj <- scale( test.data[,paste(behavvar)] - test_fitted_covars )
    
    # keep only IVs and DV of interest
    ROIs <- ROIs_full
    train.data <- train.data[, c("behav_resids", ROIs)]
    test.data <- test.data[, c("behav_adj", ROIs)]
    
    # Setup a grid range of lambda values:
    lambdas <- 10^seq(-3, 3, length = 25)
    alphas <- seq(0,1,length = 15)
    
    ## predict with elastic net regression; !!! or should i use tuneGrid =expand.grid(alpha=seq(0,1,length=10) ...
    train.data <- train.data %>% select(-Left.WM.hypointensities, -Right.WM.hypointensities, -Left.non.WM.hypointensities, -Right.non.WM.hypointensities)
    enet <- train(
      as.formula(paste("behav_resids", "~ .")), data = train.data, method = "glmnet",
      trControl = trainControl("cv", number = 10),
      tuneGrid = expand.grid(alpha = alphas, lambda = lambdas)
    )
    
    ## save out everything
    outname <- paste0(gsub(" ", "_", gsub(":","_",date())), "_", round(runif(1,100,999),0))
    save(enet, file=paste0(workdir,"/grids_",outname,".Rdata"))
    
    print(paste0('ENET CV DONE ', iter))
    
  } # end loop through brainvar
  
}




