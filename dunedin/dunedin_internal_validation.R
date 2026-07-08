# PACNI vs brainAGE
# ethan whitman
# 8/16/23

library(haven)
library(dplyr)
library(superheat)
library(ggplot2)

# load in all iterations
file_list <- list.files('/Users/ew198/Documents/brainpace/results/alpha_gridsearch/predictions')
file_dir <- '/Users/ew198/Documents/brainpace/results/alpha_gridsearch/predictions/'

# load one to get snum
load("/Users/ew198/Documents/brainpace/results/random_splits_test/predictions/predictions_Fri_Aug_11_16_28_08_2023_918.Rdata")

pred_values_df <- data.frame(snum=df$snum)
for (f in file_list){
  load(paste0(file_dir, f))
  print(head(df))
  pred_values_df <- data.frame(pred_values_df, df[,2:3])
}

ridge_preds <- pred_values_df[,c(1,seq(from=2, to=ncol(pred_values_df), by=2))]
enet_preds <- pred_values_df[,c(1,seq(from=3, to=ncol(pred_values_df), by=2))]


ridge_enet_cors <- rep(NA, 100)
for (i in 1:100){
  ridge_enet_cors[i] <- cor.test(ridge_preds[,i+1], enet_preds[,i+1])$estimate
}

unique(ridge_enet_cors)

# load PACNI model predictions

# load in model with DK+GWR post grid search
load("/Users/ew198/Documents/brainpace/results/full_model_train/dk_gwr/predictions_Wed_Sep_20_13_44_09_2023_533.Rdata")
final_pred_dk_gwr <- df


colnames(final_pred_dk_gwr) <- c("snum", "dk_gwr_pred")


# load validation phenotypes
phenotypes_elife <- as.data.frame(read_sav("/Users/ew198/Documents/brainpace/data/internal_validation/phenotypes_elife_2023.sav"))
load('/Users/ew198/Documents/methylation/DunedinClocks45.rdata')
dunedinpace <- DunedinClocks45

# load pace of aging
load('/Users/ew198/Documents/brainpace/data/dunedin/Dunedin_wholebodyaging_Nov2023.rdata')


pacni_behavior <- left_join(final_pred_dk_gwr, phenotypes_elife, by = 'snum')
pacni_behavior <- left_join(pacni_behavior, Dunedin_wholebodyaging_Nov2023[,c('snum', 'PaceOfAgingP45')], by = 'snum')




# cor.test(median_pred_behavior$prediction_ridge, median_pred_behavior$PaceOfAgingP45)
# cor.test(median_pred_behavior$prediction_ridge, median_pred_behavior$brainAgeGap45_Liem17_ctrd)

models <- c('dk_gwr_pred', 'PaceOfAgingP45')
behavvars <- c('Health45', 'ZFacialAge45', 'balClsMax45', 'ChairStands45', 'StepPlace45', 'GripRelUnit45', 'GripDom45', 'GpDomTim45', 'PhyLimts45',
               'Velocity_avg45', 'vci45a', 'pri45a', 'wmi45a', 'psi45a', 'fsiq45a', 'IQ_ResChange.x')

brainage_comp_results <- data.frame(behavvar=character(), comb=numeric(), model=character(), beta=numeric(), error=numeric(), rqsuare=numeric(), p=numeric())

for(model in models){
  for(behavvar in behavvars){
    pacni_model <- summary(lm(scale(get(behavvar))~scale(get(model))+sex, data=pacni_behavior ) )
    brainage_comp_results <- rbind(brainage_comp_results, 
                                   data.frame(behavvar=rep(behavvar, 4), model=c(model),
                                              comb=c(0,0,1,1),
                                              beta=c(pacni_beta=pacni_model$coefficients[2,1]),
                                              error=c(pacni_beta=pacni_model$coefficients[2,2]),
                                              rsquare=c(pacni_beta=pacni_model$adj.r.squared),
                                              p=c(pacni_beta=pacni_model$coefficients[2,4]))
                                )
  }
}

brainage_comp_results$ci_l <- brainage_comp_results$beta - (1.96*brainage_comp_results$error)
brainage_comp_results$ci_h <- brainage_comp_results$beta + (1.96*brainage_comp_results$error)


# comparison with PoA
cor.test(pacni_behavior$PaceOfAgingP45, pacni_behavior$dk_gwr_pred)

colfunc <- colorRampPalette(c('#197EC0FF','#71D0F5FF', '#FD8CC1FF', '#FD8CC1FF', '#C80813FF', '#C80813FF'))

poa_pacni_cor <- ggplot(data=pacni_behavior, 
       aes(x=PaceOfAgingP45, 
           y=((dk_gwr_pred - mean(dk_gwr_pred)) / sd(dk_gwr_pred) * 0.3 + 1)) ) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'black')+
  geom_point(aes(color = pacni_behavior$PaceOfAgingP45), size = 2) +
  geom_smooth(method='lm', color = 'black')+
  xlim(.35, 2.5)+
  ylim(.35, 2.5)+
  labs(y = "DunedinPACNI",
       x = "Pace of Aging") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20, hjust = .5),
        title = element_text(size = 20),
        legend.position = 'none') +
  scale_color_gradientn(colors = c(colfunc(5)))



### PACNI vs PoA bar chart


brainage_comp_results_pacni_poa <- bag_comp_results[bag_comp_results$model %in% c('dk_gwr_pred', 'PaceOfAgingP45'),]
brainage_comp_results_pacni_poa <- brainage_comp_results_pacni_poa[brainage_comp_results_pacni_poa$behavvar %in% 
                                                           c('balClsMax45', 'Velocity_avg45', 'StepPlace45', 'ChairStands45', 'GripRelUnit45', 'GpDomTim45', 'PhyLimts45', 'fsiq45a', 'vci45a', 'pri45a', 'wmi45a', 'psi45a','IQ_ResChange.x', 'Health45', 'ZFacialAge45'),]

custom_colors <- c("dk_gwr_pred" = '#1A9993FF', "PaceOfAgingP45" = '#8A9197FF')
custom_behav_order <- rev(c('balClsMax45', 'Velocity_avg45', 'StepPlace45', 'ChairStands45', 'GripRelUnit45', 'GpDomTim45', 'PhyLimts45',  'Health45',
                            'fsiq45a', 'vci45a', 'pri45a', 'wmi45a', 'psi45a','IQ_ResChange.x', 
                            'ZFacialAge45'))

behav_pretty_order <- rev(c('Balance', 'Gait speed', 'Step in place', 'Chair stands', 'Grip strength', 'Motor coordination', 'Physical limitations', 'Self-reported health', 
                        'Full-scale IQ', 'Verbal reasoning', 'Perceptual reasoning', 'Working memory', 'Processing speed','IQ decline', 
                        'Facial aging'))

brainage_comp_results_pacni_poa$behavvar <- factor(
  brainage_comp_results_pacni_poa$behavvar,
  levels = custom_behav_order
)

pacni_poa_bar <- ggplot(data = brainage_comp_results_pacni_poa, 
       aes(x = behavvar, y = abs(beta), group = model, color = model,
           ymin = abs(beta) - 1.96 * error, ymax = abs(beta) + 1.96 * error)) +
  geom_pointrange(position = position_dodge(width = 0.8), stat = "identity", size =.7) +
  ylab(expression(italic(β)))+
  geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5),
             color = c('#8A9197FF', rep('white', 5), '#8A9197FF',  rep('white', 7)), 
             alpha = c(rep(.5, 14)))+
  scale_color_manual(values = custom_colors, labels = c("DunedinPACNI", "Pace of Aging")) +
  scale_x_discrete(labels = behav_pretty_order) +
  geom_hline(yintercept=0)+
  coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 0),
        legend.position='none')

brainage_comp_results_pacni_poa$ci_l <- brainage_comp_results_pacni_poa$beta - (1.96*brainage_comp_results_pacni_poa$error)
brainage_comp_results_pacni_poa$ci_h <- brainage_comp_results_pacni_poa$beta + (1.96*brainage_comp_results_pacni_poa$error)

# write.csv(brainage_comp_results_pacni_poa, file = '/Users/ew198/Documents/brainpace/results/dunedin_validation/brainage_comp_results_pacni_poa.csv')

# 1700 x 760
ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_1.png", 
       poa_pacni_cor + pacni_poa_bar, 
       width = 17, height = 7.6, units = 'in', dpi = 500)

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_1_l.png", pacni_poa_bar, width = 8.5, height = 7.6, units = 'in', dpi = 500)


