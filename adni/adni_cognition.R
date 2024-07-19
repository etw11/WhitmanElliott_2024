# relationship between PACNI and cognitive scores in ADNI
# 11/7/23
# ethan whitman

library(ggplot2)
library(gridExtra)
library(beepr)
library(colorspace)
library(ggsci)

# load ADNI PACNI scores
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

# adding columns 
adni_adsp_pacni_brainager$SCANDATE <- adni_adsp_pacni_brainager$date_scanned
adni_adsp_pacni_brainager$RID_DATE <- adni_adsp_pacni_brainager$date_scanned

#################################################
############## PACNI AND COGNITION ##############
#################################################

### get cognitive scores from ADNIMERGE
adnimerge <- read.csv("/Users/ew198/Documents/brainpace/data/adni/ADNIMERGE_29May2023.csv")
adnimerge_cognitive <- adnimerge[,c('RID', 'EXAMDATE', 'ADAS11', 'ADAS13', 'ADASQ4', 'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_forgetting', 'RAVLT_perc_forgetting', 'LDELTOTAL', 'DIGITSCOR', 'TRABSCOR', 'FAQ', 'MOCA')]
adnimerge_cognitive <- adnimerge_cognitive[adnimerge_cognitive$RID %in% adni_adsp_pacni_brainager$RID,]

# cog testing QC filtering

# trails
TRABSCOR_outlier_index <- (scale(adnimerge_cognitive$TRABSCOR) > 4) + (scale(adnimerge_cognitive$TRABSCOR) < -4)
adnimerge_cognitive$TRABSCOR[as.logical(TRABSCOR_outlier_index)] <- NA
adnimerge_cognitive$TRABSCOR[adnimerge_cognitive$TRABSCOR == 0] <- NA

# RAVLT_forgetting
RAVLT_forgetting_outlier_index <- (scale(adnimerge_cognitive$RAVLT_forgetting) > 4) + (scale(adnimerge_cognitive$RAVLT_forgetting) < -4)
adnimerge_cognitive$RAVLT_forgetting[as.logical(RAVLT_forgetting_outlier_index)] <- NA

# RAVLT_perc_forgetting
RAVLT_perc_forgetting_outlier_index <- (scale(adnimerge_cognitive$RAVLT_perc_forgetting) > 4) + (scale(adnimerge_cognitive$RAVLT_perc_forgetting) < -4)
adnimerge_cognitive$RAVLT_perc_forgetting[as.logical(RAVLT_perc_forgetting_outlier_index)] <- NA


## as with DX, look for most proximate test, either before or after scan
## unlike DX, test must be within 6 months, and can only be assigned to one scan
## just starting a new N=6204 dataframe without all the parcel data, will add columns for each paired test

data_merged <- adni_adsp_pacni_brainager[,c("RID","RID_DATE","sex","age_at_scan","date_scanned","SCANDATE","pacni","bag","pacni_scaled","bag_scaled","dx_recode")] 
tests <- c('ADAS13', 'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting', 'LDELTOTAL', 'DIGITSCOR', 'TRABSCOR', 'FAQ', 'MOCA')
for(test in tests){
  
  print(test)
  
  data_merged[,paste(test)] <- NA
  data_merged[,paste0(test,"_interval")] <- NA
  data_merged[,paste0(test,"_date")] <- NA
  
  for(s in unique(data_merged$RID)){ # loop thru individuals
    
    img_cur <- data_merged %>% filter(RID==s & is.na(get(test))) %>% dplyr::select(c("SCANDATE", paste(test))) %>% arrange(SCANDATE) # scans for current subj
    cog_cur <- adnimerge %>% filter(RID==s & !is.na(get(test))) %>% dplyr::select(c(paste(test),"EXAMDATE")) %>% arrange(EXAMDATE) # tests for current subj
    
    if( nrow(cog_cur) > 0 & nrow(img_cur) > 0 ){
      
      cog_cur$pairedScanDate <- NA
      img_cur_unpaired <- img_cur
      
      for(thr in c(30,60,183)) { # first require a closer match, then increase intervals
        
        # loop thru cog and match scans to each one
        for(date in cog_cur[is.na(cog_cur$pairedScanDate),"EXAMDATE"]){ 
          
          if(nrow(img_cur_unpaired)==0){ break } # if there aren't any unpaired scans left, we are done
          
          closestImg <- which.min( abs(as.Date(img_cur_unpaired$SCANDATE) - as.Date(date)) )
          
          if( abs(as.Date(img_cur_unpaired[closestImg,"SCANDATE"]) - as.Date(date)) < thr ){ # make sure less than <thr> apart
            
            row <- which(data_merged$RID==s & data_merged$SCANDATE==img_cur_unpaired[closestImg,"SCANDATE"]) # position in final data frame
            data_merged[row, paste(test)] <- cog_cur[cog_cur$EXAMDATE==date, paste(test)]
            data_merged[row, paste0(test,"_interval")] <- as.Date(img_cur_unpaired[closestImg,"SCANDATE"]) - as.Date(date) 
            data_merged[row, paste0(test,"_date")] <- date
            cog_cur[which(cog_cur$EXAMDATE==date), "pairedScanDate"] <- as.Date(img_cur_unpaired[closestImg,"SCANDATE"]) # just for checks
            img_cur_unpaired <- img_cur_unpaired[-closestImg, ]
            
          }
          
        } # end loop thru cog 
        # # check for ties ... maybe come back to this, prob needs to be adapted
        # if( sum( abs(as.Date(cog_cur$EXAMDATE) - as.Date(date)) == abs(data_merged[i, paste0(test,"_interval")]) ) > 1 ) { print (paste("Tie for",i, RID_cur, date, cog_cur )) }
      } # end loop thru thr
      
      # some checks
      img_cur <- data_merged %>% filter(RID==s) %>% dplyr::select("SCANDATE", paste0(test,"_date"), paste(test), paste0(test,"_interval") )
      ## no clue why date formatting doesn't work one row at a time in above loop, do it here so i can easily look at the dates
      cog_cur$pairedScanDate <- as.Date(cog_cur$pairedScanDate) 
      img_cur[,paste0(test,"_date")] <- as.Date(img_cur[,paste0(test,"_date")])
      
      ## check if a cog test is closer to a scan that it's not paired with
      cog_cur$closestScanDate <- NA
      for(i in 1:nrow(cog_cur)) {
        cog_cur[i,"closestScanDate"] <- as.Date( img_cur[which.min( abs(as.Date(cog_cur[i,"EXAMDATE"]) - as.Date(img_cur$SCANDATE)) ), "SCANDATE"] )
      }
      cog_cur$closestScanDate <- as.Date(cog_cur$closestScanDate) # still no clue why this doesn't work one row at a time
      cog_cur$closestInterval <- cog_cur$closestScanDate - as.Date(cog_cur$EXAMDATE)
      cog_cur$pairedInterval <- cog_cur$pairedScanDate - as.Date(cog_cur$EXAMDATE)
      
      if( sum( (abs(cog_cur$closestInterval) < 183) & (abs(cog_cur$closestInterval) < abs(cog_cur$pairedInterval)), na.rm=TRUE ) > 0 ){
        print(paste(s, "has",test,"closer to a scan it's not paired with"))
      }      
      
      ## if there are unpaired tests with scans within 6 months, check whether it's already paired with another scan
      if( sum( abs(cog_cur[is.na(cog_cur$pairedInterval), "closestInterval"]) < 183 ) > 0 ) {
        for(i in 1:nrow(img_cur)){
          if(is.na(img_cur[i,paste0(test,"_date")])){
            closestCog <- which.min ( abs(as.Date(img_cur[i,"SCANDATE"]) - as.Date(cog_cur$EXAMDATE)) )
            if((closestCog>0) & abs(as.Date(img_cur[i,"SCANDATE"]) - as.Date(cog_cur[closestCog, "EXAMDATE"])) < 183){
              if( cog_cur[closestCog, "EXAMDATE"] %in% img_cur[,paste0(test,"_date")] ) {
                print(paste(s, "does not have",test,"paired with it but there is a paired test within 6 mos"))
              } else {
                print(paste(s, "does not have",test,"paired with it but there is an UNPAIRED test within 6 mos"))
              }
            }
          }
        }
      }
    } # end if this sub has data for this test
  } # end loop thru subs
} # end loop thru tests


### run associations
adni_adsp_cog_test_prox <- data_merged
adni_adsp_cog_test_prox <- adni_adsp_cog_test_prox[!is.na(adni_adsp_cog_test_prox$bag),]

cog_assoc_results_test <- data.frame(test=character(), model=character(), comb=numeric(), n=numeric(), time_diff=numeric(), beta=numeric(), error=numeric(), rqsuare=numeric())

tests <- c('ADAS13', 'MMSE', 'RAVLT_immediate', 'RAVLT_learning', 'RAVLT_perc_forgetting', 'LDELTOTAL', 'DIGITSCOR', 'TRABSCOR', 'FAQ', 'MOCA')
for(test in tests){
  
    # subsetting data
    data_temp <- adni_adsp_cog_test_prox
    data_temp$var <- scale(data_temp[,c(test)])
    data_temp <- data_temp[!is.na(data_temp$var), ]
    data_temp$time_diff <- data_temp[,c(paste0(test, "_interval"))]
    
    # PACNI
    pacni_model_orig <- summary(lm(var~pacni_scaled+factor(sex)+age_at_scan, data=data_temp ) )
    pm1 <- plm(var~pacni_scaled+factor(sex)+age_at_scan, 
               model = "pooling", 
               data = pdata.frame(data_temp), 
               na.action = na.omit)
    G <- length(unique(data_temp$RID))
    N <- length(data_temp$RID)
    dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
    firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
    pacni_model.plm <- coeftest(pm1, vcov = firm_c_vcov)
    
    
    cog_assoc_results_test <- rbind(cog_assoc_results_test, 
                               data.frame(test=test,
                                          model='pacni',
                                          comb = 0,
                                          n_obs=nrow(data_temp),
                                          n_ind=length(unique(data_temp$RID)),
                                          time_diff=median(abs(data_temp$time_diff)),
                                          beta=pacni_model.plm[2,1],
                                          error=pacni_model.plm[2,2],
                                          p=pacni_model.plm[2,4],
                                          rsquare=pacni_model_orig$adj.r.squared))
    
    # BAG
    bag_model_orig <- summary(lm(var~bag_scaled+factor(sex)+age_at_scan, data=data_temp ) )
    pm1 <- plm(var~bag_scaled+factor(sex)+age_at_scan, 
               model = "pooling", 
               data = pdata.frame(data_temp), 
               na.action = na.omit)
    G <- length(unique(data_temp$RID))
    N <- length(data_temp$RID)
    dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
    firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
    bag_model.plm <- coeftest(pm1, vcov = firm_c_vcov)
    
    cog_assoc_results_test <- rbind(cog_assoc_results_test, 
                                    data.frame(test=test,
                                               model='bag',
                                               comb = 0,
                                               n_obs=nrow(data_temp),
                                               n_ind=length(unique(data_temp$RID)),
                                               time_diff=median(abs(data_temp$time_diff)),
                                               beta=bag_model.plm[2,1],
                                               error=bag_model.plm[2,2],
                                               p=bag_model.plm[2,4],
                                               rsquare=bag_model_orig$adj.r.squared))

    
    ## combined results
    comb_model_orig <- summary(lm(var~pacni_scaled+bag_scaled+factor(sex)+age_at_scan, data=data_temp ) )
    pm1 <- plm(var~pacni_scaled+bag_scaled+factor(sex)+age_at_scan, 
               model = "pooling", 
               data = pdata.frame(data_temp), 
               na.action = na.omit)
    G <- length(unique(data_temp$RID))
    N <- length(data_temp$RID)
    dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
    firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
    comb_model.plm <- coeftest(pm1, vcov = firm_c_vcov)
    
    cog_assoc_results_test <- rbind(cog_assoc_results_test, 
                                    data.frame(test=test,
                                               model='pacni',
                                               comb = 1,
                                               n_obs=nrow(data_temp),
                                               n_ind=length(unique(data_temp$RID)),
                                               time_diff=median(abs(data_temp$time_diff)),
                                               beta=comb_model.plm[2,1],
                                               error=comb_model.plm[2,2],
                                               p=comb_model.plm[2,4],
                                               rsquare=comb_model_orig$adj.r.squared))
    
    
    cog_assoc_results_test <- rbind(cog_assoc_results_test, 
                                    data.frame(test=test,
                                               model='bag',
                                               comb = 1,
                                               n_obs=nrow(data_temp),
                                               n_ind=length(unique(data_temp$RID)),
                                               time_diff=median(abs(data_temp$time_diff)),
                                               beta=comb_model.plm[3,1],
                                               error=comb_model.plm[3,2],
                                               p=comb_model.plm[3,4],
                                               rsquare=comb_model_orig$adj.r.squared))
    
}


cog_assoc_results_test$ci_l <- cog_assoc_results_test$beta - (1.96*cog_assoc_results_test$error)
cog_assoc_results_test$ci_h <- cog_assoc_results_test$beta + (1.96*cog_assoc_results_test$error)

