# ethan whitman
# 9/6/24
# reading brainlat data

library(DunedinPACNI)
library(stringr)
library(car)

brainlat_df <- LoadFreeSurferStats(fsdir = '/Users/ew198/Documents/brainpace/data/brainlat/freesurfer',
                                   sublist = '/Users/ew198/Documents/brainpace/data/brainlat/gid.csv')

brainlat_pacni <- ExportDunedinPACNI(data = brainlat_df,
                                     outdir = '/Users/ew198/Documents/brainpace/data/brainlat')


# check demog data

brainlat_demog <- read.csv('/Users/ew198/Documents/brainpace/data/brainlat/demog/BrainLat_Demographic_MRI.csv')

brainlat_pacni_demog_all <- left_join(brainlat_pacni, brainlat_demog, by = 'ID') # N = 486
brainlat_pacni_demog <- brainlat_pacni_demog_all[!is.na(brainlat_pacni_demog_all$Age),] # 52 people missing age
brainlat_pacni_demog$pacni_resid <- scale(resid(lm(pacni~Age+sex, brainlat_pacni_demog)))
brainlat_pacni_demog$country <- gsub("sub\\-", "", brainlat_pacni_demog$ID)
brainlat_pacni_demog$country <- str_replace_all(brainlat_pacni_demog$country, "[:digit:]", "")

# read brainageR

brainlat_bag <- read.csv('/Users/ew198/Documents/brainpace/data/brainlat/brainageR/BrainLAT_brainageR.csv')
brainlat_pacni_bag <- left_join(brainlat_pacni_demog, brainlat_bag, by = 'ID')
brainlat_pacni_bag$bag <- c(brainlat_pacni_bag$brain.predicted_age - brainlat_pacni_bag$Age)

# filter by diagnosis
brainlat_pacni_bag$diagnosis <- factor(brainlat_pacni_bag$diagnosis, levels = c('CN', 'AD', 'FTD', 'MS', 'PD'))
brainlat_pacni_bag_dem <- brainlat_pacni_bag[brainlat_pacni_bag$diagnosis %in% c('CN', 'AD', 'FTD'),]

# descriptive plots
cor.test(brainlat_pacni_bag_dem$pacni, brainlat_pacni_bag_dem$bag)

brainlat_pacni_bag_scatter <- ggplot(data=brainlat_pacni_bag_dem, 
                                     aes(x=pacni, y=bag)) +
  geom_point() +
  geom_smooth(method='lm', color = 'red') +
  labs(x = "DunedinPACNI",
       y = "Brain age gap") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15))

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/figs_bl_pacni_bag_scatter.png", 
       brainlat_pacni_bag_scatter,
       width = 4, height = 4, units = 'in', dpi = 500)

cor.test(brainlat_pacni_bag_dem$pacni, brainlat_pacni_bag_dem$Age)

brainlat_pacni_age <- ggplot(data=brainlat_pacni_bag_dem, 
                             aes(x=Age, y=pacni)) +
  geom_point() +
  geom_smooth(method='lm', color = 'red') +
  labs(x = "Age",
       y = "DunedinPACNI") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  scale_x_continuous(breaks = c(50, 70, 90),
                       labels = c(50, 70, 90))

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/figs_bl_pacni_age_scatter.png", 
       brainlat_pacni_age,
       width = 4, height = 4, units = 'in', dpi = 500)



# associations with diangosis
brainlat_pacni_dx <- summary(lm(scale(pacni)~diagnosis+Age+sex, brainlat_pacni_bag_dem))
brainlat_bag_dx <- summary(lm(scale(bag)~diagnosis+Age+sex, brainlat_pacni_bag_dem))

brainlat_pacni_dx_comb <- summary(lm(scale(pacni)~diagnosis+Age+sex+bag, brainlat_pacni_bag_dem))
brainlat_bag_dx_comb <- summary(lm(scale(bag)~diagnosis+Age+sex+pacni, brainlat_pacni_bag_dem))

brainlat_pacni_bag_dem$pacni_resid <- scale(resid(lm(pacni~Age+sex, brainlat_pacni_bag_dem)))
brainlat_pacni_bag_dem$bag_resid <- scale(resid(lm(bag~Age+sex, brainlat_pacni_bag_dem)))


brainlat_pacni_dx_plot <- ggplot(data = brainlat_pacni_bag_dem,
       aes(y = pacni_resid, x = factor(diagnosis, levels = c('CN', 'AD', 'FTD')),
           fill = factor(diagnosis))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7)+
  labs(x = 'Diagnosis', y = 'DunedinPACNI', fill = 'Diagnosis') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 0),
        title = element_text(size = 0),
        legend.text = element_text(size = 15),
        legend.position = 'none')+ 
  scale_fill_manual(values = c('#709AE1FF', '#C80813FF', '#C80813FF'),
                    labels = c('CN', 'AD', 'FTD'))


# read cognition data

brainlat_cog <- read.csv('/Users/ew198/Documents/brainpace/data/brainlat/demog/BrainLat_Cognition_MRI.csv')

brainlat_pacni_cog <- left_join(brainlat_pacni_bag_dem, brainlat_cog, by = c('ID' = 'MRI_ID'))

brainlat_pacni_moca <- summary(lm(scale(moca_total)~scale(pacni)+Age+sex, brainlat_pacni_cog))
brainlat_bag_moca <- summary(lm(scale(moca_total)~scale(bag)+Age+sex, brainlat_pacni_cog))

brainlat_pacni_moca_comb <- summary(lm(scale(moca_total)~scale(pacni)+Age+sex+scale(bag), brainlat_pacni_cog))
brainlat_bag_moca_comb <- summary(lm(scale(moca_total)~scale(bag)+Age+sex+scale(pacni), brainlat_pacni_cog))

brainlat_moca_table <- data.frame(test=rep(c('MoCA'),4),
                                model=c(rep('PACNI', 2), rep('BAG', 2)),
                                comb=c(0,1,0,1),
                                beta=c(brainlat_pacni_moca$coefficients[2,1], brainlat_pacni_moca_comb$coefficients[2,1],
                                       brainlat_bag_moca$coefficients[2,1], brainlat_bag_moca_comb$coefficients[2,1]),
                                
                                se=c(brainlat_pacni_moca$coefficients[2,2], brainlat_pacni_moca_comb$coefficients[2,2],
                                     brainlat_bag_moca$coefficients[2,2], brainlat_bag_moca_comb$coefficients[2,2]))

brainlat_moca_table$ci_l <- c(brainlat_moca_table$beta - (1.96*brainlat_moca_table$se))
brainlat_moca_table$ci_h <- c(brainlat_moca_table$beta + (1.96*brainlat_moca_table$se))

#save(brainlat_moca_table, file = '/Users/ew198/Documents/brainpace/results/brainlat/brainlat_moca_table.Rdata')
#write.csv(brainlat_moca_table, file = '/Users/ew198/Documents/brainpace/results/brainlat/brainlat_moca_table.csv')

brainlat_pacni_cog_mocasample <- brainlat_pacni_cog[!is.na(brainlat_pacni_cog$moca_total),]
brainlat_pacni_cog_mocasample$pacni_scaled <- scale(brainlat_pacni_cog_mocasample$pacni)
brainlat_pacni_cog_mocasample$bag_scaled <- scale(brainlat_pacni_cog_mocasample$bag)

# get combined beta coefficients
pacni_bag_moca_model <- lm(scale(moca_total)~pacni_scaled+Age+sex+bag_scaled, brainlat_pacni_cog_mocasample)

# additive moca effects
brainlat_moca_add_table <- data.frame(test=rep(c('MoCA'),3),
                                  model=c('PACNI', 'BAG', 'PACNI+BAG'),
                                  beta=c(brainlat_pacni_moca$coefficients[2,1],
                                         brainlat_bag_moca$coefficients[2,1],
                                         as.numeric(deltaMethod(pacni_bag_moca_model, "pacni_scaled+bag_scaled")[1])),
                                  
                                  se=c(brainlat_pacni_moca$coefficients[2,2],
                                       brainlat_bag_moca$coefficients[2,2],
                                       as.numeric(deltaMethod(pacni_bag_moca_model, "pacni_scaled+bag_scaled")[2])))

brainlat_moca_add_table$ci_l <- c(brainlat_moca_add_table$beta - (1.96*brainlat_moca_add_table$se))
brainlat_moca_add_table$ci_h <- c(brainlat_moca_add_table$beta + (1.96*brainlat_moca_add_table$se))

#save(brainlat_moca_add_table, file = '/Users/ew198/Documents/brainpace/results/brainlat/brainlat_moca_add_table.Rdata')

# ADNI and BrainLat overlay

load('/Users/ew198/Documents/brainpace/data/adni/adni_adsp_moca_prox.Rdata')
adni_adsp_moca_prox$pacni_resid <- scale(resid(lm(pacni~age_at_scan+sex, adni_adsp_moca_prox)))

moca_overlay <- ggplot() +
  geom_point(data = adni_adsp_moca_prox, aes(x = MOCA, y = pacni_resid), color = '#075149FF')+
  geom_point(data = brainlat_pacni_cog, aes(x = moca_total, y = pacni_resid), color = '#FD7446FF')+
  geom_smooth(data = adni_adsp_moca_prox, aes(x = MOCA, y = pacni_resid), method = 'lm', color = '#075149FF')+
  geom_smooth(data = brainlat_pacni_cog, aes(x = moca_total, y = pacni_resid), method = 'lm', color = '#FD7446FF')+
  labs(x = 'MoCA', y = 'DunedinPACNI') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))

# make forest plot

brainlat_dx_table <- data.frame(test=rep(c('AD', 'FTD', 'MS'),2),
                               model=c(rep('PACNI', 3), rep('BAG', 3)),
                               beta=c(brainlat_pacni_dx$coefficients[2,1], brainlat_pacni_dx$coefficients[3,1], brainlat_pacni_dx$coefficients[4,1],
                                      brainlat_bag_dx$coefficients[2,1], brainlat_bag_dx$coefficients[3,1], brainlat_bag_dx$coefficients[4,1]),
                               
                               se=c(brainlat_pacni_dx$coefficients[2,2], brainlat_pacni_dx$coefficients[3,2], brainlat_pacni_dx$coefficients[4,2],
                                    brainlat_bag_dx$coefficients[2,2], brainlat_bag_dx$coefficients[3,2], brainlat_bag_dx$coefficients[4,2]))

brainlat_dx_table$ci_l <- c(brainlat_dx_table$beta - (1.96*brainlat_dx_table$se))
brainlat_dx_table$ci_h <- c(brainlat_dx_table$beta + (1.96*brainlat_dx_table$se))

# comb forest plot

brainlat_dx_table <- data.frame(test=rep(c('AD', 'FTD'),4),
                                model=c(rep('PACNI', 4), rep('BAG', 4)),
                                comb=rep(c(0,0,1,1),2),
                                beta=c(brainlat_pacni_dx$coefficients[2,1], brainlat_pacni_dx$coefficients[3,1],
                                       brainlat_pacni_dx_comb$coefficients[2,1], brainlat_pacni_dx_comb$coefficients[3,1],
                                       brainlat_bag_dx$coefficients[2,1], brainlat_bag_dx$coefficients[3,1],
                                       brainlat_bag_dx_comb$coefficients[2,1], brainlat_bag_dx_comb$coefficients[3,1]),
                                
                                se=c(brainlat_pacni_dx$coefficients[2,2], brainlat_pacni_dx$coefficients[3,2],
                                     brainlat_pacni_dx_comb$coefficients[2,2], brainlat_pacni_dx_comb$coefficients[3,2],
                                     brainlat_bag_dx$coefficients[2,2], brainlat_bag_dx$coefficients[3,2],
                                     brainlat_bag_dx_comb$coefficients[2,2], brainlat_bag_dx_comb$coefficients[3,2]))

brainlat_dx_table$ci_l <- c(brainlat_dx_table$beta - (1.96*brainlat_dx_table$se))
brainlat_dx_table$ci_h <- c(brainlat_dx_table$beta + (1.96*brainlat_dx_table$se))

#save(brainlat_dx_table, file = '/Users/ew198/Documents/brainpace/results/brainlat/brainlat_dx_table.Rdata')
#write.csv(brainlat_dx_table, file = '/Users/ew198/Documents/brainpace/results/brainlat/brainlat_dx_table.csv')


# make forest plot comparing ADNI and BrainLat effects

load('/Users/ew198/Documents/brainpace/results/adni/cognition/group_pacni_betas_adni_adsp.Rdata')

brainlat_comp_table <- data.frame(test=c('Dementia (ADNI)', 'AD (BrainLat)', 'FTD (BrainLat)'),
                                  dataset=c('ADNI', 'BrainLat', 'BrainLat'),
                                beta=c(group_pacni_betas_adni_adsp[3,1], brainlat_pacni_dx$coefficients[2,1], brainlat_pacni_dx$coefficients[3,1]),
                                
                                se=c(group_pacni_betas_adni_adsp[3,6], brainlat_pacni_dx$coefficients[2,2], brainlat_pacni_dx$coefficients[3,2]))

brainlat_comp_table$ci_l <- c(brainlat_comp_table$beta - (1.96*brainlat_comp_table$se))
brainlat_comp_table$ci_h <- c(brainlat_comp_table$beta + (1.96*brainlat_comp_table$se))

brainlat_comp_plot <- ggplot(brainlat_comp_table, 
       aes(x=factor(test, levels = c('Dementia (ADNI)', 'AD (BrainLat)', 'FTD (BrainLat)')), 
           y=beta,  color = dataset)) + 
  geom_pointrange(aes(ymin=ci_l, ymax=ci_h), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(y = 'SDs above healthy controls') +
  ylim(-.1, 1.05)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = 'none') +
  scale_color_manual(labels = c('ADNI', 'BrainLat'),
                     values = c("#075149FF", "#FD7446FF"))


ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/brainlat_fig.png", 
       brainlat_pacni_dx_plot + brainlat_comp_plot + moca_overlay, 
       width = 13, height = 4, units = 'in', dpi = 500)




# no GWR

brainlat_df_nogwr <- LoadFreeSurferStats(fsdir = '/Users/ew198/Documents/brainpace/data/brainlat/freesurfer',
                                   sublist = '/Users/ew198/Documents/brainpace/data/brainlat/gid.csv',
                                   missing_gwr = TRUE)

brainlat_pacni_nogwr <- ExportDunedinPACNI(data = brainlat_df_nogwr,
                                     outdir = '/Users/ew198/Documents/brainpace/data/brainlat',
                                     missing_gwr = TRUE)

colnames(brainlat_pacni_nogwr) <- c('ID', 'pacni_nogwr')

brainlat_pacni_bag_dem_nogwr <- left_join(brainlat_pacni_bag_dem, brainlat_pacni_nogwr, by = 'ID')


# GWR and no GWR correlation
cor.test(brainlat_pacni_bag_dem_nogwr$pacni, brainlat_pacni_bag_dem_nogwr$pacni_nogwr)

bl_gwr_nogwr_scatter <- ggplot(data=brainlat_pacni_bag_dem_nogwr, 
                            aes(x=scale(pacni), y=scale(pacni_nogwr))) +
  geom_point() +
  geom_smooth(method='lm', color = 'red')+
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ylim(-3.5, 5)+
  xlim(-3.5, 5)+
  labs(x = "DunedinPACNI with GWR",
       y = "DunedinPACNI without GWR") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15))

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/figs_bl_gwr_nogwr_scatter.png", 
       bl_gwr_nogwr_scatter,
       width = 4, height = 4, units = 'in', dpi = 500)



