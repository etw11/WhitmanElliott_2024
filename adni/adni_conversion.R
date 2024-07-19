# finding how many ADNI ADSP participants converted to MCI/dementia
# ethan whitman
# 1/10/24

library(lubridate)
library(survival)
library(ggfortify)
library(survival)
library(survminer)
library(cowplot)
library(intsurv)
library(patchwork)

# load in ADN PACNI scores
load('/Users/ew198/Documents/brainpace/data/adni/inhouse_organized/adni_adsp_pacni.Rdata')

# load brainageR scores from ADNI
adni_brainager <- read.csv('/Users/ew198/Documents/brainpace/data/adni/brainageR.csv')
adni_brainager$RID <- adni_brainager$RID_DATE
adni_brainager$date_scanned <- adni_brainager$RID_DATE

adni_brainager$RID <- gsub("_.*", "", adni_brainager$RID)
adni_brainager$date_scanned <- gsub(".*_", "", adni_brainager$date_scanned)
adni_brainager$date_scanned <- as.Date(adni_brainager$date_scanned, format = "%Y%m%d")

# merge
adni_adsp_pacni_brainager <- merge(adni_adsp_pacni, adni_brainager[c('RID', 'brain.predicted_age', 'date_scanned')], by=c('RID', 'date_scanned'))

# bag calculation
adni_adsp_pacni_brainager$bag <- as.numeric(adni_adsp_pacni_brainager$brain.predicted_age) - adni_adsp_pacni_brainager$age_at_scan

# pacni age residualization
adni_adsp_pacni_brainager$pacni_resid <- scale(resid(lm(pacni~age_at_scan, data = adni_adsp_pacni_brainager)))

# scale
adni_adsp_pacni_brainager$pacni_scaled <- scale(adni_adsp_pacni_brainager$pacni)
adni_adsp_pacni_brainager$bag_scaled <- scale(adni_adsp_pacni_brainager$bag)

# baseline PACNI values
adni_adsp_pacni_baseline_temp <- adni_adsp_pacni_brainager[order(adni_adsp_pacni_brainager$RID, adni_adsp_pacni_brainager$age_at_scan),]
adni_adsp_pacni_baseline <- adni_adsp_pacni_baseline_temp[!duplicated(adni_adsp_pacni_baseline_temp$RID),]

# only pulling participants in ADSP 
adnimerge <- read.csv("/Users/ew198/Documents/brainpace/data/adni/ADNIMERGE_29May2023.csv")
adnimerge_adsp_inc <- adnimerge[adnimerge$RID %in% adni_adsp_precombat$RID,]

# dividing into CN stable, CN -> MCI, CN/MCI -> AD

long_adsp_dx <- data.frame(RID = unique(adni_adsp_pacni_baseline$RID), 
                           group = c(rep(NA, nrow(adni_adsp_pacni_baseline))),
                           dx_date = c(rep(NA, nrow(adni_adsp_pacni_baseline))))
adnimerge_adsp_inc <- adnimerge[adnimerge$RID %in% adni_adsp_precombat$RID,]

x <- 1
for (i in unique(adni_adsp_pacni_baseline$RID)){
  temp_df <- as.data.frame(adnimerge_adsp_inc[adnimerge_adsp_inc$RID == i,])
  temp_df$DX[temp_df$DX == ""] <- NA
  temp_df <- temp_df[!is.na(temp_df$DX),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'CN' && temp_df$DX[timepoints] == 'CN'){
      long_adsp_dx$group[x] <- 'sCN'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'CN' && temp_df$DX[timepoints] == 'MCI'){
      long_adsp_dx$group[x] <- 'CN_to_MCI'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'MCI' && temp_df$DX[timepoints] == 'CN'){
      long_adsp_dx$group[x] <- 'MCI_to_CN'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'CN' && temp_df$DX[timepoints] == 'Dementia'){
      long_adsp_dx$group[x] <- 'CN_to_Dem'
      log_index <- temp_df$DX == 'Dementia'
      long_adsp_dx$dx_date[x] <- temp_df$EXAMDATE[which(log_index)[1]]
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'MCI' && temp_df$DX[timepoints] == 'MCI'){
      long_adsp_dx$group[x] <- 'sMCI'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'MCI' && temp_df$DX[timepoints] == 'Dementia'){
      long_adsp_dx$group[x] <- 'MCI_to_Dem'
      log_index <- temp_df$DX == 'Dementia'
      long_adsp_dx$dx_date[x] <- temp_df$EXAMDATE[which(log_index)[1]]
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'Dementia' && temp_df$DX[timepoints] == 'Dementia'){
      long_adsp_dx$group[x] <- 'sDem'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'Dementia' && temp_df$DX[timepoints] == 'MCI'){
      long_adsp_dx$group[x] <- 'Dem_to_MCI'
    } else if (adni_adsp_pacni_baseline$DX[adni_adsp_pacni_baseline$RID == i] == 'Dementia' && temp_df$DX[timepoints] == 'CN'){
      long_adsp_dx$group[x] <- 'Dem_to_CN'
    }
    if (long_adsp_dx$group[x] != 'MCI_to_Dem' | long_adsp_dx$group[x] != 'CN_to_Dem'){
      long_adsp_dx$dx_date[x] <-  temp_df$EXAMDATE[timepoints]
    }
  }
  x <- x+1
}


# order by baseline
adni_adsp_or_index <- with(adni_adsp_pacni_brainager, order(RID, date_scanned))
adni_adsp_precombat_or <- adni_adsp_pacni_brainager[adni_adsp_or_index,]                          
adni_adsp_precombat_or_baseline <- adni_adsp_precombat_or[!duplicated(adni_adsp_precombat_or$RID),]                          

adni_adsp_precombat_or_baseline_conv <- merge(adni_adsp_precombat_or_baseline, long_adsp_dx, by = 'RID')

# calculating age and time to event
adni_adsp_precombat_or_baseline_conv$age_at_baseline_yr <- adni_adsp_precombat_or_baseline_conv$age_at_baseline/12

adni_adsp_precombat_or_baseline_conv$age_at_finaldx <- (adni_adsp_precombat_or_baseline_conv$age_at_baseline + 
                                                          interval(ymd(adni_adsp_precombat_or_baseline_conv$baseline_date), ymd(adni_adsp_precombat_or_baseline_conv$dx_date)) %/% months(1))

adni_adsp_precombat_or_baseline_conv$age_at_finaldx_yr <- adni_adsp_precombat_or_baseline_conv$age_at_finaldx/12

adni_adsp_precombat_or_baseline_conv$time_to_event <- adni_adsp_precombat_or_baseline_conv$age_at_finaldx_yr - adni_adsp_precombat_or_baseline_conv$age_at_scan

# dividing into dx groups
adni_adsp_precombat_or_baseline_conv$dem_conversion <- as.numeric(dplyr::recode(adni_adsp_precombat_or_baseline_conv$group, 
                                                                                'sMCI' = '0', 'sCN' = '0',  'MCI_to_Dem' = '1', 'CN_to_Dem'='1'))
adni_adsp_precombat_or_baseline_conv$mci_conversion <- as.numeric(dplyr::recode(adni_adsp_precombat_or_baseline_conv$group, 
                                                                                'sCN' = '0', 'CN_to_Dem'='1', 'CN_to_MCI'='1'))
adni_adsp_precombat_or_baseline_conv$any_conversion <- as.numeric(dplyr::recode(adni_adsp_precombat_or_baseline_conv$group, 
                                                                                'sMCI' = '0', 'sCN' = '0', 'MCI_to_Dem' = '1', 'CN_to_Dem'='1', 'CN_to_MCI'='1'))
adni_adsp_precombat_or_baseline_conv$mcidem_conversion <- as.numeric(dplyr::recode(adni_adsp_precombat_or_baseline_conv$group, 
                                                                                   'sMCI' = '0', 'MCI_to_Dem' = '1'))


# removing people without follow up
adni_adsp_paired_pacni_dx_baseline_conv_pt <- adni_adsp_precombat_or_baseline_conv[adni_adsp_precombat_or_baseline_conv$time_to_event >0,]


# pacni
# CN
pacni_cox_cn <- coxph(Surv(time_to_event, mci_conversion) ~ pacni_scaled + sex + age_at_scan,
                      data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mci_conversion),])

# MCI
pacni_cox_mci <- coxph(Surv(time_to_event, mcidem_conversion) ~ pacni_scaled + sex + age_at_scan,
                       data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mcidem_conversion),])

# BAG
# CN
bag_cox_cn <- coxph(Surv(time_to_event, mci_conversion) ~ bag_scaled + sex + age_at_scan,
                    data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mci_conversion),])

# MCI
bag_cox_mci <- coxph(Surv(time_to_event, mcidem_conversion) ~ bag_scaled + sex + age_at_scan,
                     data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mcidem_conversion),])


# COMB
# CN
comb_cox_cn <- coxph(Surv(time_to_event, mci_conversion) ~ pacni_scaled + bag_scaled + sex + age_at_scan,
                     data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mci_conversion),])

# MCI
comb_cox_mci <- coxph(Surv(time_to_event, mcidem_conversion) ~ pacni_scaled + bag_scaled + sex + age_at_scan,
                      data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mcidem_conversion),])

# controlling for APOE4
# pacni
# CN
pacni_cox_cn_apoe <- coxph(Surv(time_to_event, mci_conversion) ~ pacni_scaled + sex + age_at_scan + factor(APOE4),
                           data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mci_conversion),])
# MCI
pacni_cox_mci_apoe <- coxph(Surv(time_to_event, mcidem_conversion) ~ pacni_scaled + sex + age_at_scan + factor(APOE4),
                            data=adni_adsp_paired_pacni_dx_baseline_conv_pt[!is.na(adni_adsp_paired_pacni_dx_baseline_conv_pt$mcidem_conversion),])



