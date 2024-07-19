# ADNI longitudinal analysis
# ethan whitman
# 2/9/24

library(lme4)
library(nlme)
library(viridis)
library(longCombat)
library(dplyr)

# load PACNI scores from ADNI
load('/Users/ew198/Documents/brainpace/data/adni/inhouse_organized/adni_adsp_pacni.Rdata')

# load brainageR scores from ADNI
adni_brainager <- read.csv('/Users/ew198/Documents/brainpace/data/adni/brainageR.csv')
adni_brainager$RID <- adni_brainager$RID_DATE
adni_brainager$date_scanned <- adni_brainager$RID_DATE

adni_brainager$RID <- gsub("_.*", "", adni_brainager$RID)
adni_brainager$date_scanned <- gsub(".*_", "", adni_brainager$date_scanned)
adni_brainager$date_scanned <- as.Date(adni_brainager$date_scanned, format = "%Y%m%d")

adni_brainager$brain.predicted_age <- as.numeric(adni_brainager$brain.predicted_age)

# merge
adni_adsp_pacni_brainager <- merge(adni_adsp_pacni, adni_brainager[c('RID', 'brain.predicted_age', 'date_scanned')], by=c('RID', 'date_scanned'))

# bag calculation
adni_adsp_pacni_brainager$bag <- adni_adsp_pacni_brainager$brain.predicted_age - adni_adsp_pacni_brainager$age_at_scan

# scale
adni_adsp_pacni_brainager$pacni_scaled <- scale(adni_adsp_pacni_brainager$pacni)
adni_adsp_pacni_brainager$bag_scaled <- scale(adni_adsp_pacni_brainager$bag)

# restrict to 2 or more timepoints
RID_repeat_index <- names(table(adni_adsp_pacni_brainager$RID)[table(adni_adsp_pacni_brainager$RID) > 1])
adni_adsp_pacni_brainager_long <- adni_adsp_pacni_brainager[adni_adsp_pacni_brainager$RID %in% RID_repeat_index,]

# check that people don't change sites or scanner strength

for (i in c(unique(adni_adsp_pacni_brainager_long$RID))){
  temp_df <- adni_adsp_pacni_brainager_long[adni_adsp_pacni_brainager_long$RID == i ,c('RID', 'site', 'field_strength')]
  print(i)
  if (length(unique(temp_df$site)) > 1){
    print('SITE CHANGE')
  } else if (length(unique(temp_df$field_strength)) > 1){
    print('SCANNER CHANGE')
  }
}

# need to run longComBat to harmonize before next step
adni_adsp_pacni_brainager_long$dx_recode <- dplyr::recode("CN" = "0", "MCI" = "1", "Dementia" = "2", adni_adsp_pacni_brainager_long$DX)
fs_vars <- colnames(adni_adsp_pacni_brainager_long %>% dplyr::select(-RID, -date_scanned, -age_at_scan, -site, -sex, -dataset, -field_strength, -pacni, -baseline_date, -age_at_baseline, -DX, -dx_recode,
                                                              -EXAMDATE, -time_diff, -APOE4, -brain.predicted_age, -bag, -pacni_scaled, -bag_scaled))

adni_adsp_pacni_brainager_long$batch <- adni_adsp_pacni_brainager_long$site
adni_adsp_pacni_brainager_long$batch[adni_adsp_pacni_brainager_long$field_strength == '3'] <- adni_adsp_pacni_brainager_long$site[adni_adsp_pacni_brainager_long$field_strength == '3'] + 2000

adni_adsp_longcombat <- longCombat(idvar='RID', 
                          timevar='age_at_scan',
                          batchvar='batch', 
                          features=fs_vars, 
                          formula='sex',
                          ranef='(1+age_at_scan|RID)',
                          data=adni_adsp_pacni_brainager_long %>% dplyr::select(-date_scanned, -site, -dataset, -field_strength, -pacni, -baseline_date, -age_at_baseline, -DX, -dx_recode, 
                                                                         -EXAMDATE, -time_diff, -APOE4, -brain.predicted_age, -bag, -pacni_scaled, -bag_scaled))

adni_adsp_pacni_brainager_long_combat <- left_join(adni_adsp_pacni_brainager_long, 
                                            adni_adsp_longcombat$data_combat %>% dplyr::select(-batch),
                                            by = c('RID', 'age_at_scan'))

# order by scandate
adni_adsp_pacni_brainager_long_combat <- adni_adsp_pacni_brainager_long_combat[order(adni_adsp_pacni_brainager_long_combat$RID, adni_adsp_pacni_brainager_long_combat$age_at_scan),]

# generate growth curves
adni_adsp_pacni_brainager_long_combat$hc_bilat <- rowMeans(data.frame(left=adni_adsp_pacni_brainager_long_combat$Hippocampus_Left.combat, 
                                                                      right=adni_adsp_pacni_brainager_long_combat$Hippocampus_Right.combat))
hc_bilat_curve <-  coef(lmer(scale(hc_bilat)~1 + (scale(age_at_scan)|RID), data = adni_adsp_pacni_brainager_long_combat))$RID[,1]

# get observation length
mri_obs_length <- rep(NA, length(unique(adni_adsp_pacni_brainager_long_combat$RID)))
baseline_age <- rep(NA, length(unique(adni_adsp_pacni_brainager_long_combat$RID)))
x <- 1
for (i in unique(adni_adsp_pacni_brainager_long_combat$RID)){
  temp_df <- adni_adsp_pacni_brainager_long_combat[adni_adsp_pacni_brainager_long_combat$RID == i,]
  mri_obs_length[x] <- max(temp_df$age_at_scan) - min(temp_df$age_at_scan)
  baseline_age[x] <- min(temp_df$age_at_scan)
  x <- x + 1
}

adni_adsp_randomslopes <- data.frame(RID = unique(adni_adsp_pacni_brainager_long_combat$RID), 
                                     sex = adni_adsp_pacni_brainager_long_combat$sex[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     age_at_scan = adni_adsp_pacni_brainager_long_combat$age_at_scan[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     dx_recode =  adni_adsp_pacni_brainager_long_combat$dx_recode[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     APOE4 =  adni_adsp_pacni_brainager_long_combat$APOE4[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     DX =  adni_adsp_pacni_brainager_long_combat$DX[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     obs_length = mri_obs_length,
                                     hc_bilat_curve = hc_bilat_curve, 
                                     baseline_pacni_scaled = adni_adsp_pacni_brainager_long_combat$pacni_scaled[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     baseline_pacni_resid = adni_adsp_pacni_brainager_long_combat$pacni_resid[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     baseline_bag_scaled = adni_adsp_pacni_brainager_long_combat$bag_scaled[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)],
                                     baseline_bag_resid = adni_adsp_pacni_brainager_long_combat$bag_resid[!duplicated(adni_adsp_pacni_brainager_long_combat$RID)])

# PACNI predicting decline

# HC
summary(lm(scale(hc_bilat_curve)~baseline_pacni_scaled+factor(sex)+age_at_scan+obs_length, data = adni_adsp_randomslopes))
summary(lm(scale(hc_bilat_curve)~baseline_bag_scaled+factor(sex)+age_at_scan+obs_length, data = adni_adsp_randomslopes))

# comb model
summary(lm(scale(hc_bilat_curve)~baseline_pacni_scaled+baseline_bag_scaled+factor(sex)+age_at_scan+obs_length, data = adni_adsp_randomslopes))

# controlling for APOE
summary(lm(scale(hc_bilat_curve)~baseline_pacni_scaled+factor(sex)+age_at_scan+obs_length+factor(APOE4), data = adni_adsp_randomslopes))
summary(lm(scale(hc_bilat_curve)~baseline_bag_scaled+factor(sex)+age_at_scan+obs_length+factor(APOE4), data = adni_adsp_randomslopes))
summary(lm(scale(hc_bilat_curve)~baseline_bag_scaled+baseline_pacni_scaled+factor(sex)+age_at_scan+obs_length+factor(APOE4), data = adni_adsp_randomslopes))

