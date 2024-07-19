# getting ADNI education for grace
# 5/13/24


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

# residualize
adni_adsp_pacni_brainager$pacni_resid <- scale(resid(lm(pacni~age_at_baseline+DX+factor(sex), data = adni_adsp_pacni_brainager) ))
adni_adsp_pacni_brainager$bag_resid <- rep(NA, nrow(adni_adsp_pacni_brainager))
adni_adsp_pacni_brainager$bag_resid[!is.na(adni_adsp_pacni_brainager$bag)] <- scale(resid(lm(bag~age_at_baseline+DX+factor(sex), data = adni_adsp_pacni_brainager) ))

# merge in education from ADNIMERGE

adnimerge <- read.csv("/Users/ew198/Documents/brainpace/data/adni/ADNIMERGE_29May2023.csv")
adni_edu <- adnimerge[!duplicated(adnimerge$RID), c('RID', 'PTEDUCAT')]
adni_adsp_pacni_brainager <- merge(adni_adsp_pacni_brainager, adni_edu, by = 'RID')

adni_adsp_pacni_brainager$edu_scaled <- scale(adni_adsp_pacni_brainager$PTEDUCAT)

adni_adsp_pacni_brainager_firstobs <- dplyr::arrange(adni_adsp_pacni_brainager, RID, age_at_scan)
adni_adsp_pacni_brainager_firstobs <- adni_adsp_pacni_brainager_firstobs[!duplicated(adni_adsp_pacni_brainager_firstobs$RID),]
adni_adsp_pacni_brainager_firstobs_bag <- adni_adsp_pacni_brainager_firstobs[!is.na(adni_adsp_pacni_brainager_firstobs$bag),]

# analysis
pacni_edu <- summary(lm(scale(PTEDUCAT)~pacni_scaled+factor(sex)+age_at_scan, adni_adsp_pacni_brainager_firstobs))
bag_edu <- summary(lm(scale(PTEDUCAT)~bag_scaled+factor(sex)+age_at_scan, adni_adsp_pacni_brainager_firstobs))
comb_edu <- summary(lm(scale(PTEDUCAT)~pacni_scaled+bag_scaled+factor(sex)+age_at_scan, adni_adsp_pacni_brainager_firstobs))

# BAG sample
pacni_edu_bagsample <- summary(lm(scale(PTEDUCAT)~pacni_scaled+factor(sex)+age_at_scan, adni_adsp_pacni_brainager_firstobs_bag))




