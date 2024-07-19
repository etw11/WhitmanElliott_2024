# testing pacni functions with UKB data

library(tidyverse)
##load functions
source("/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/scripts/PACNI_Rpackage_240518/LoadFreeSurferStats_fast.R")
source("/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/scripts/PACNI_Rpackage_240518/DunedinHarmonization.R")
source("/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/scripts/PACNI_Rpackage_240518/ExportDunedinPACNI.R")

###Load UKB phenos
imagingPhenos <- read_csv("/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/ukb_pacni_pheno_imagingSubs.csv")
##Load list of all FS dirs
subList <- read_csv("/ncf/sba08/Triarchy_OA/Max/UKBanalyses/superBrainAging/lists/UKB_allFSdirs_T1.csv",col_names = F)

#Remove subs missing stats files in FS directories
subList <- subList %>% filter(!X1 %in% c("1353005_20263_2_0","2202905_20263_2_0", "2899188_20263_2_0", "4190529_20263_2_0", "5023579_20263_2_0"))
#Remove subs missing some ROIs in stats files
subList <- subList %>% filter(!X1 %in% c("5746349_20263_2_0"))



##Rename for cov file
imagingPhenos$site <- imagingPhenos$f.54.2.0
imagingPhenos$ID <- paste0(imagingPhenos$f.eid,"_20263_2_0")
#Setup birthday by combining month and year
imagingPhenos$birthMonthNumber <- match(imagingPhenos$f.52.0.0, month.name)
imagingPhenos <- imagingPhenos %>% mutate(birthday = ymd(paste0(f.34.0.0,"-",as.numeric(birthMonthNumber),"-15"))) ##Assume day 15 to be less biased
imagingPhenos<-imagingPhenos %>% mutate(ageAtScan = as.numeric((f.53.2.0 - birthday)/365.25))
imagingPhenos$age <- imagingPhenos$ageAtScan
imagingPhenos$sex <- imagingPhenos$f.31.0.0

covFile <- imagingPhenos %>% select(ID, age, sex, site) %>% filter(ID %in% subList$X1)
imagingPhenos <- imagingPhenos %>% filter(ID %in% subList$X1)
write_csv(covFile, "/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/covariates.csv")
write_rds(imagingPhenos, "/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/ukb_allPhenos_imagingSubs_240605.rds")

# running functions on UKB baseline FS data downloaded from AMS
date()
df <- LoadFreeSurferStats_fast(fsdir = '/ncf/cnl07/Users/UKB/freesurfer/data/',
                          covdir = '/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/',
                          missing_gwr = TRUE)
date()
#Took ~50 minutes to load and 
newName.df <- df
df.save <- df

newName.df$sex<-tolower(newName.df$sex) #Male and female have to be lower case 
newName.df$sex[newName.df$sex == "female"] <- "1"
newName.df$sex[newName.df$sex == "male"] <- "2"

df <- newName.df

data_pacni <- ExportDunedinPACNI(data = df,
                                 modeldir = '/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/scripts/PACNI_Rpackage_240518/',
                                 outdir = '/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/')

data_pacni_noGWR <- ExportDunedinPACNI(data = df,
                                 modeldir = '/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/scripts/PACNI_Rpackage_240518/',
                                 outdir = '/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/',
                                 harmonization = FALSE,
                                 missing_gwr = T)

write_csv(data_pacni,"/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/pacni_ukb_fs6_240605.csv")
write_csv(data_pacni_noGWR,"/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/generatePACNI_ukb_240605/pacni_ukb_fs6_noHarmonization_noGWR_240605.csv")

# comparing with previously run PACNI versions
data_pacni_FS6_new<-data_pacni
data_pacni_FS6_new$ID <- as.numeric(gsub("_20263_2_0", "", data_pacni_FS6_new$ID))

#orig FS 6 PACNI
load('/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/ukb_pacniMon_Oct__9_17_02_35_2023_355.Rdata')
data_pacni_FS6_old<-ukb_postcombat

# homemade FS 7 PACNI
load('/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/ukb_pacniTue_Apr__9_08_59_44_2024_340.Rdata')
data_pacni_FS7_old<-ukb_postcombat

##Combine
data_pacni_FS6_old$pacni_FS6_old<-data_pacni_FS6_old$pacni
data_pacni_FS6_old <- data_pacni_FS6_old %>% select(eid,pacni_FS6_old)
data_pacni_FS7_old$pacni_FS7<-data_pacni_FS7_old$pacni
data_pacni_FS7_old <- data_pacni_FS7_old %>% select(eid,pacni_FS7)


data_pacni_all<-full_join(data_pacni_FS6_new,data_pacni_FS7_old, join_by(ID == eid))
data_pacni_all<-full_join(data_pacni_all,data_pacni_FS6_old, join_by(ID == eid))
#write_csv(data_pacni,"/ncf/sba08/Triarchy_OA/Max/UKBanalyses/PACNI/data/ukb_PACNI_oldFS6andFS7_240417.csv")

###Exact same as my version...
cor.test(data_pacni_all$pacni, data_pacni_all$pacni_FS7)
cor.test(data_pacni_all$pacni, data_pacni_all$pacni_FS6_old)
cor.test(data_pacni_all$pacni_FS7, data_pacni_all$pacni_FS6_old)

ggplot(data_pacni, 
       aes(x=pacni, y=pacni_FS7)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)+
  labs(x = "FS6",
       y = "FS7") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))


# testing with ADNI stats files


df <- LoadFreeSurferStats(fsdir = '/Users/ew198/Documents/brainpace/data/adni/freesurfer_stats/',
                          covdir = '/Users/ew198/Documents/brainpace/data/adni/pacni_package/')

data_harmonized <- DunedinHarmonization(data = df,
                                        cov_list = c('sex', 'age', 'dx'),
                                        synthdir = "/Users/ew198/Documents/brainpace/results/synth/",
                                        outdir = '/Users/ew198/Documents/brainpace/data/adni/pacni_package/')

data_pacni <- ExportDunedinPACNI(data = data_harmonized,
                                 modeldir = '/Users/ew198/Documents/brainpace/scripts/pacni_package/',
                                 outdir = "/Users/ew198/Documents/brainpace/data/adni/pacni_package/")


