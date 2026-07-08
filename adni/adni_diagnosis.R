# ADNI brainageR cross sectional association comparison
# ethan whitman
# 1/10/23

library(plm)
library(lmtest)
library(reshape2)
library(ggridges)

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
adni_adsp_pacni_brainager$bag <- as.numeric(adni_adsp_pacni_brainager$brain.predicted_age) - adni_adsp_pacni_brainager$age_at_scan

# scale
adni_adsp_pacni_brainager$pacni_scaled <- scale(adni_adsp_pacni_brainager$pacni)
adni_adsp_pacni_brainager$bag_scaled <- scale(adni_adsp_pacni_brainager$bag)

# residualize
adni_adsp_pacni_brainager$pacni_resid <- scale(resid(lm(pacni~age_at_baseline+DX+factor(sex), data = adni_adsp_pacni_brainager) ))
adni_adsp_pacni_brainager$bag_resid <- rep(NA, nrow(adni_adsp_pacni_brainager))
adni_adsp_pacni_brainager$bag_resid[!is.na(adni_adsp_pacni_brainager$bag)] <- scale(resid(lm(bag~age_at_baseline+DX+factor(sex), data = adni_adsp_pacni_brainager) ))

adni_adsp_pacni$dx_recode <- dplyr::recode("CN" = "0", "MCI" = "1", "Dementia" = "2", adni_adsp_pacni$DX)

# BAG and pacni correlation
cor.test(adni_adsp_pacni_brainager$pacni, adni_adsp_pacni_brainager$bag)
cor.test(adni_adsp_pacni_brainager$age_at_scan, adni_adsp_pacni_brainager$pacni)

adni_adsp_pacni_brainager_bag <- adni_adsp_pacni_brainager[!is.na(adni_adsp_pacni_brainager$bag),]

# DX group differences
summary(lm(scale(pacni)~as.factor(dx_recode)+as.factor(sex)+age_at_scan, data = adni_adsp_pacni))

# pacni robust errors
summary(lm(scale(pacni)~as.factor(dx_recode)+as.factor(sex)+age_at_scan, data = adni_adsp_pacni))
#adni_adsp_pacni_brainager$dx_recode_mci_ref <- factor(adni_adsp_pacni_brainager$dx_recode, levels = c("1", "0", "2"))
group_pacni_betas_adni_adsp <- as.data.frame(summary(lm(pacni_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan, data = adni_adsp_pacni_brainager))$coefficients)
group_pacni_betas_adni_adsp$dx <- c('int', 'MCI', 'Dementia', 'sex', 'age_at_scan')
colnames(group_pacni_betas_adni_adsp)[1:2] <- c("beta", "std_error")

pm1 <- plm(pacni_scaled~as.factor(dx_recode)+sex+as.numeric(age_at_scan), model = "pooling", data = pdata.frame(adni_adsp_pacni_brainager), na.action = na.omit)
G <- length(unique(adni_adsp_pacni_brainager$RID))
N <- length(adni_adsp_pacni_brainager$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
adsp.pacni.dx.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

group_pacni_betas_adni_adsp$robust_error <- adsp.pacni.dx.plm.res[,2]
group_pacni_betas_adni_adsp$ci_l <- group_pacni_betas_adni_adsp$beta - (1.96*group_pacni_betas_adni_adsp$robust_error)
group_pacni_betas_adni_adsp$ci_h <- group_pacni_betas_adni_adsp$beta + (1.96*group_pacni_betas_adni_adsp$robust_error)

# BAG robust errors
summary(lm(bag_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan, data = adni_adsp_pacni_brainager_bag))
group_brainage_betas_adni_adsp <- as.data.frame(summary(lm(bag_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan, data = adni_adsp_pacni_brainager_bag))$coefficients)
group_brainage_betas_adni_adsp$dx <- c('int', 'MCI', 'Dementia', 'sex', 'age_at_scan')
colnames(group_brainage_betas_adni_adsp)[1:2] <- c("beta", "std_error")

pm1 <- plm(bag_scaled~as.factor(dx_recode)+sex+as.numeric(age_at_scan), model = "pooling", data = pdata.frame(adni_adsp_pacni_brainager_bag), na.action = na.omit)
G <- length(unique(adni_adsp_pacni_brainager_bag$RID))
N <- length(adni_adsp_pacni_brainager_bag$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
adsp.bag.dx.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

group_brainage_betas_adni_adsp$robust_error <- adsp.bag.dx.plm.res[,2]
group_brainage_betas_adni_adsp$ci_l <- group_brainage_betas_adni_adsp$beta - (1.96*group_brainage_betas_adni_adsp$robust_error)
group_brainage_betas_adni_adsp$ci_h <- group_brainage_betas_adni_adsp$beta + (1.96*group_brainage_betas_adni_adsp$robust_error)

# PACNI combined robust errors
summary(lm(pacni_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan+bag_scaled, data = adni_adsp_pacni_brainager_bag))
group_pacni_bagres_betas_adni_adsp <- as.data.frame(summary(lm(pacni_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan+bag_scaled, data = adni_adsp_pacni_brainager_bag))$coefficients)
group_pacni_bagres_betas_adni_adsp$dx <- c('int', 'MCI', 'Dementia', 'sex', 'age_at_scan', 'bag_scaled')
colnames(group_pacni_bagres_betas_adni_adsp)[1:2] <- c("beta", "std_error")

pm1 <- plm(pacni_scaled~as.factor(dx_recode)+sex+as.numeric(age_at_scan)+bag_scaled, model = "pooling", data = pdata.frame(adni_adsp_pacni_brainager_bag), na.action = na.omit)
G <- length(unique(adni_adsp_pacni_brainager_bag$RID))
N <- length(adni_adsp_pacni_brainager_bag$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
adsp.pacni.bres.dx.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

group_pacni_bagres_betas_adni_adsp$robust_error <- adsp.pacni.bres.dx.plm.res[,2]
group_pacni_bagres_betas_adni_adsp$ci_l <- group_pacni_bagres_betas_adni_adsp$beta - (1.96*group_pacni_bagres_betas_adni_adsp$robust_error)
group_pacni_bagres_betas_adni_adsp$ci_h <- group_pacni_bagres_betas_adni_adsp$beta + (1.96*group_pacni_bagres_betas_adni_adsp$robust_error)

# BAG combined robust errors
summary(lm(bag_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan+pacni_scaled, data = adni_adsp_pacni_brainager_bag))
group_brainage_pacnires_betas_adni_adsp <- as.data.frame(summary(lm(bag_scaled~as.factor(dx_recode)+as.factor(sex)+age_at_scan+pacni_scaled, data = adni_adsp_pacni_brainager_bag))$coefficients)
group_brainage_pacnires_betas_adni_adsp$dx <- c('int', 'MCI', 'Dementia', 'sex', 'age_at_scan', 'pacni_scaled')
colnames(group_brainage_pacnires_betas_adni_adsp)[1:2] <- c("beta", "std_error")

pm1 <- plm(bag_scaled~as.factor(dx_recode)+sex+as.numeric(age_at_scan)+pacni_scaled, model = "pooling", data = pdata.frame(adni_adsp_pacni_brainager_bag), na.action = na.omit)
G <- length(unique(adni_adsp_pacni_brainager_bag$RID))
N <- length(adni_adsp_pacni_brainager_bag$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
adsp.bag.pres.dx.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

group_brainage_pacnires_betas_adni_adsp$robust_error <- adsp.bag.pres.dx.plm.res[,2]
group_brainage_pacnires_betas_adni_adsp$ci_l <- group_brainage_pacnires_betas_adni_adsp$beta - (1.96*group_brainage_pacnires_betas_adni_adsp$robust_error)
group_brainage_pacnires_betas_adni_adsp$ci_h <- group_brainage_pacnires_betas_adni_adsp$beta + (1.96*group_brainage_pacnires_betas_adni_adsp$robust_error)



