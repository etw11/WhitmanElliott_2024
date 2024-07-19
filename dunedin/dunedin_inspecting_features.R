# inspecting PACNI weights using DK parcellation
# ethan whitman
# 11/20/23

library(ggplot2)
library(freesurferformats)
library(dplyr)
library(gridExtra)
library(grid)
library(ciftiTools)
ciftiTools.setOption("wb_path", "/Applications/workbench/bin_macosx64/wb_command")
library(ggseg)
library(ggpubr)

# load surfaces
lh_annot <- read.fs.annot('/Users/ew198/freesurfer/subjects/fsaverage5/label/lh.aparc.annot')
rh_annot <- read.fs.annot('/Users/ew198/freesurfer/subjects/fsaverage5/label/rh.aparc.annot')
lh_inflated <- read_surf('/Users/ew198/freesurfer/subjects/fsaverage5/surf/lh.inflated.gii')
rh_inflated <- read_surf('/Users/ew198/freesurfer/subjects/fsaverage5/surf/rh.inflated.gii')

# load haufe transformed feature importance
load('/Users/ew198/Documents/brainpace/results/full_model_train/dk_gwr/predictions_Wed_Sep_20_13_44_09_2023_533.Rdata')
dk_gwr_df <- data.frame(ROI=ROIs_full,
                    fi = enet_coefs_haufe)

# filter ROIs set to 0 by elastic net
load('/Users/ew198/Documents/brainpace/results/full_model_train/dk_gwr/dk_gwr_model_Wed_Sep_20_13_44_09_2023_533.Rdata')
dk_gwr_model <- enet
coefs_dk_gwr <- data.frame(ROI=rownames(coef(dk_gwr_model$finalModel, dk_gwr_model$bestTune$lambda)),
                          coef=as.numeric(coef(dk_gwr_model$finalModel, dk_gwr_model$bestTune$lambda)))

dk_gwr_df_full <- dplyr::left_join(dk_gwr_df, coefs_dk_gwr, by = "ROI")

dk_gwr_df_full$fi_filter <- rep(NA, nrow(dk_gwr_df_full))
dk_gwr_df_full$fi_filter[dk_gwr_df_full$coef != 0] <- dk_gwr_df_full$fi[dk_gwr_df_full$coef != 0]

################ DK+GWR ###############

# ridge bar graph
dk_gwr_bar_table <- dk_gwr_df_full[order(dk_gwr_df_full$fi_filter),]
colfunc <- colorRampPalette(c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'))
dk_gwr_bar_table$pcolor <- colfunc(nrow(dk_gwr_bar_table))
dk_gwr_bar_table$ROI <- factor(dk_gwr_bar_table$ROI, levels = dk_gwr_bar_table[order(dk_gwr_bar_table$fi_filter),]$ROI)
dk_gwr_bar_table$cat <- rep(NA, nrow(dk_gwr_bar_table))
dk_gwr_bar_table$cat[grepl('GMV', dk_gwr_bar_table$ROI)] <-  'GMV'
dk_gwr_bar_table$cat[grepl('SA_', dk_gwr_bar_table$ROI)] <-  'SA'
dk_gwr_bar_table$cat[grepl('CT_', dk_gwr_bar_table$ROI)] <-  'CT'
dk_gwr_bar_table$cat[grepl('GWR_', dk_gwr_bar_table$ROI)] <-  'GWR'
dk_gwr_bar_table$cat[is.na(dk_gwr_bar_table$cat)] <-  'ASEG'

custom_colors <- c("ASEG" = "#F0E442", "CT" = "#E69F00", "GMV" = "#56B4E9", 'GWR'='#009E73', 'SA'='#999999')
ggplot(data=dk_gwr_bar_table[!is.na(dk_gwr_bar_table$fi_filter),], 
       aes(x=ROI, y=fi_filter, fill=cat)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = custom_colors) +  # Set custom colors
  ylim(-.2, .2)+
  theme_classic()  +
  labs(y = "Feature Importance", title = "Ridge") +
  theme(axis.text.x = element_text(size = 0),
        axis.text.y = element_text(vjust = 0.5, size = 15),
        axis.title.y = element_text(size = 20, vjust = 0),
        axis.title.x = element_text(size = 0),
        legend.text = element_text(size=15),
        title = element_text(size = 0))


# cortical thickness

ctgwr_ourlabs_rl <- data.frame(labs=as.character(dk_gwr_bar_table$ROI[dk_gwr_bar_table$cat=='CT']),
                               fi_filter=dk_gwr_bar_table$fi_filter[dk_gwr_bar_table$cat=='CT'],
                            pcolor=dk_gwr_bar_table$pcolor[dk_gwr_bar_table$cat=='CT'])

# left
ctgwr_ourlabs_l <- ctgwr_ourlabs_rl[!grepl('_right', ctgwr_ourlabs_rl$labs),]
ctgwr_ourlabs_l$labs <- gsub("CT_", "", ctgwr_ourlabs_l$labs)
ctgwr_ourlabs_l$labs <- gsub("G.S", "G_and_S", ctgwr_ourlabs_l$labs)
ctgwr_ourlabs_l <- data.frame(labs=c(ctgwr_ourlabs_l$labs[1:(nrow(ctgwr_ourlabs_l))], 'Medial_wall'),
                              fi_filter=c(ctgwr_ourlabs_l$fi_filter[1:(nrow(ctgwr_ourlabs_l))],0),
                           pcolor=c(ctgwr_ourlabs_l$pcolor[1:(nrow(ctgwr_ourlabs_l))],0))

# right
ctgwr_ourlabs_r <- ctgwr_ourlabs_rl[!grepl('_left', ctgwr_ourlabs_rl$labs),]
ctgwr_ourlabs_r$labs <- gsub("CT_", "", ctgwr_ourlabs_r$labs)
ctgwr_ourlabs_r$labs <- gsub("G.S", "G_and_S", ctgwr_ourlabs_r$labs)
ctgwr_ourlabs_r <- data.frame(labs=c(ctgwr_ourlabs_r$labs[1:(nrow(ctgwr_ourlabs_r))], 'Medial_wall'),
                              fi_filter=c(ctgwr_ourlabs_r$fi_filter[1:(nrow(ctgwr_ourlabs_r))],0),
                           pcolor=c(ctgwr_ourlabs_r$pcolor[1:(nrow(ctgwr_ourlabs_r))],0))


# left
ct_ourlabs_l_test <- ctgwr_ourlabs_l
colnames(ct_ourlabs_l_test) <- c('region', 'fi_filter', 'pcolor')
ct_ourlabs_l_test$region <- sub("_left", "", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub(".1", "", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("rostral", "rostral ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("superior", "superior ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("inferior", "inferior ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("anterior", "anterior ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("posterior", "posterior ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("medial", "medial ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("lateral", "lateral ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("pars", "pars ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("caudal", "caudal ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("middle", "middle ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("isthmus", "isthmus ", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("frontalpole", "frontal pole", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("temporalpole", "temporal pole", ct_ourlabs_l_test$region)
ct_ourlabs_l_test$region <- sub("transversetemporal", "transverse temporal", ct_ourlabs_l_test$region)

ct_left <- ggseg(ct_ourlabs_l_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'left',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Cortical Thickness')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

# right
ct_ourlabs_r_test <- ctgwr_ourlabs_r
colnames(ct_ourlabs_r_test) <- c('region', 'fi_filter', 'pcolor')
ct_ourlabs_r_test$region <- sub("_right", "", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub(".1", "", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("rostral", "rostral ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("superior", "superior ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("inferior", "inferior ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("anterior", "anterior ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("posterior", "posterior ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("medial", "medial ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("lateral", "lateral ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("pars", "pars ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("caudal", "caudal ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("middle", "middle ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("isthmus", "isthmus ", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("frontalpole", "frontal pole", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("temporalpole", "temporal pole", ct_ourlabs_r_test$region)
ct_ourlabs_r_test$region <- sub("transversetemporal", "transverse temporal", ct_ourlabs_r_test$region)

ct_right <- ggseg(ct_ourlabs_r_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'right',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Cortical Thickness')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')


# cortical surface area

sagwr_ourlabs_rl <- data.frame(labs=as.character(dk_gwr_bar_table$ROI[dk_gwr_bar_table$cat=='SA']),
                            fi_filter=dk_gwr_bar_table$fi_filter[dk_gwr_bar_table$cat=='SA'],
                            pcolor=dk_gwr_bar_table$pcolor[dk_gwr_bar_table$cat=='SA'])

# left
sagwr_ourlabs_l <- sagwr_ourlabs_rl[!grepl('_right', sagwr_ourlabs_rl$labs),]
sagwr_ourlabs_l$labs <- gsub("SA_", "", sagwr_ourlabs_l$labs)
sagwr_ourlabs_l$labs <- gsub("G.S", "G_and_S", sagwr_ourlabs_l$labs)
sagwr_ourlabs_l <- data.frame(labs=c(sagwr_ourlabs_l$labs[1:(nrow(sagwr_ourlabs_l))], 'Medial_wall'),
                              fi_filter=c(sagwr_ourlabs_l$fi_filter[1:(nrow(sagwr_ourlabs_l))],0),
                           pcolor=c(sagwr_ourlabs_l$pcolor[1:(nrow(sagwr_ourlabs_l))],0))

# right
sagwr_ourlabs_r <- sagwr_ourlabs_rl[!grepl('_left', sagwr_ourlabs_rl$labs),]
sagwr_ourlabs_r$labs <- gsub("SA_", "", sagwr_ourlabs_r$labs)
sagwr_ourlabs_r$labs <- gsub("G.S", "G_and_S", sagwr_ourlabs_r$labs)
sagwr_ourlabs_r <- data.frame(labs=c(sagwr_ourlabs_r$labs[1:(nrow(sagwr_ourlabs_r))], 'Medial_wall'),
                              fi_filter=c(sagwr_ourlabs_r$fi_filter[1:(nrow(sagwr_ourlabs_r))],0),
                           pcolor=c(sagwr_ourlabs_r$pcolor[1:(nrow(sagwr_ourlabs_r))],0))

# left
sa_ourlabs_l_test <- sagwr_ourlabs_l
colnames(sa_ourlabs_l_test) <- c('region', 'fi_filter', 'pcolor')
sa_ourlabs_l_test$region <- sub("_left", "", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub(".1", "", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("rostral", "rostral ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("superior", "superior ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("inferior", "inferior ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("anterior", "anterior ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("posterior", "posterior ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("medial", "medial ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("lateral", "lateral ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("pars", "pars ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("caudal", "caudal ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("middle", "middle ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("isthmus", "isthmus ", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("frontalpole", "frontal pole", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("temporalpole", "temporal pole", sa_ourlabs_l_test$region)
sa_ourlabs_l_test$region <- sub("transversetemporal", "transverse temporal", sa_ourlabs_l_test$region)

sa_left <- ggseg(sa_ourlabs_l_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'left',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Surface Area')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

# right
sa_ourlabs_r_test <- sagwr_ourlabs_r
colnames(sa_ourlabs_r_test) <- c('region', 'fi_filter', 'pcolor')
sa_ourlabs_r_test$region <- sub("_right", "", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub(".1", "", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("rostral", "rostral ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("superior", "superior ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("inferior", "inferior ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("anterior", "anterior ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("posterior", "posterior ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("medial", "medial ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("lateral", "lateral ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("pars", "pars ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("caudal", "caudal ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("middle", "middle ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("isthmus", "isthmus ", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("frontalpole", "frontal pole", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("temporalpole", "temporal pole", sa_ourlabs_r_test$region)
sa_ourlabs_r_test$region <- sub("transversetemporal", "transverse temporal", sa_ourlabs_r_test$region)

sa_right <- ggseg(sa_ourlabs_r_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'right',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Surface Area')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')



# cortical GMV

gmvgwr_ourlabs_rl <- data.frame(labs=as.character(dk_gwr_bar_table$ROI[dk_gwr_bar_table$cat=='GMV']),
                                fi_filter=dk_gwr_bar_table$fi_filter[dk_gwr_bar_table$cat=='GMV'],
                             pcolor=dk_gwr_bar_table$pcolor[dk_gwr_bar_table$cat=='GMV'])

# left
gmvgwr_ourlabs_l <- gmvgwr_ourlabs_rl[!grepl('_right', gmvgwr_ourlabs_rl$labs),]
gmvgwr_ourlabs_l$labs <- gsub("GMV_", "", gmvgwr_ourlabs_l$labs)
gmvgwr_ourlabs_l$labs <- gsub("G.S", "G_and_S", gmvgwr_ourlabs_l$labs)
gmvgwr_ourlabs_l <- data.frame(labs=c(gmvgwr_ourlabs_l$labs[1:(nrow(gmvgwr_ourlabs_l))], 'Medial_wall'),
                               fi_filter=c(gmvgwr_ourlabs_l$fi_filter[1:(nrow(gmvgwr_ourlabs_l))],0),
                            pcolor=c(gmvgwr_ourlabs_l$pcolor[1:(nrow(gmvgwr_ourlabs_l))],0))

# right
gmvgwr_ourlabs_r <- gmvgwr_ourlabs_rl[!grepl('_left', gmvgwr_ourlabs_rl$labs),]
gmvgwr_ourlabs_r$labs <- gsub("GMV_", "", gmvgwr_ourlabs_r$labs)
gmvgwr_ourlabs_r$labs <- gsub("G.S", "G_and_S", gmvgwr_ourlabs_r$labs)
gmvgwr_ourlabs_r <- data.frame(labs=c(gmvgwr_ourlabs_r$labs[1:(nrow(gmvgwr_ourlabs_r))], 'Medial_wall'),
                               fi_filter=c(gmvgwr_ourlabs_r$fi_filter[1:(nrow(gmvgwr_ourlabs_r))],0),
                            pcolor=c(gmvgwr_ourlabs_r$pcolor[1:(nrow(gmvgwr_ourlabs_r))],0))

# left
gmv_ourlabs_l_test <- gmvgwr_ourlabs_l
colnames(gmv_ourlabs_l_test) <- c('region', 'fi_filter', 'pcolor')
gmv_ourlabs_l_test$region <- sub("_left", "", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub(".1", "", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("rostral", "rostral ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("superior", "superior ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("inferior", "inferior ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("anterior", "anterior ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("posterior", "posterior ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("medial", "medial ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("lateral", "lateral ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("pars", "pars ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("caudal", "caudal ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("middle", "middle ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("isthmus", "isthmus ", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("frontalpole", "frontal pole", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("temporalpole", "temporal pole", gmv_ourlabs_l_test$region)
gmv_ourlabs_l_test$region <- sub("transversetemporal", "transverse temporal", gmv_ourlabs_l_test$region)

gmv_left <- ggseg(gmv_ourlabs_l_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'left',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Gray Matter Volume')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

# right
gmv_ourlabs_r_test <- gmvgwr_ourlabs_r
colnames(gmv_ourlabs_r_test) <- c('region', 'fi_filter', 'pcolor')
gmv_ourlabs_r_test$region <- sub("_right", "", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub(".1", "", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("rostral", "rostral ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("superior", "superior ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("inferior", "inferior ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("anterior", "anterior ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("posterior", "posterior ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("medial", "medial ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("lateral", "lateral ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("pars", "pars ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("caudal", "caudal ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("middle", "middle ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("isthmus", "isthmus ", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("frontalpole", "frontal pole", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("temporalpole", "temporal pole", gmv_ourlabs_r_test$region)
gmv_ourlabs_r_test$region <- sub("transversetemporal", "transverse temporal", gmv_ourlabs_r_test$region)

gmv_right <- ggseg(gmv_ourlabs_r_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'right',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Gray Matter Volume')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

# cortical GWR

gwr_ourlabs_rl <- data.frame(labs=as.character(dk_gwr_bar_table$ROI[dk_gwr_bar_table$cat=='GWR']),
                             fi_filter=dk_gwr_bar_table$fi_filter[dk_gwr_bar_table$cat=='GWR'],
                             pcolor=dk_gwr_bar_table$pcolor[dk_gwr_bar_table$cat=='GWR'])

# left
gwr_ourlabs_l <- gwr_ourlabs_rl[!grepl('_right', gwr_ourlabs_rl$labs),]
gwr_ourlabs_l$labs <- gsub("GWR_", "", gwr_ourlabs_l$labs)
gwr_ourlabs_l$labs <- gsub("G.S", "G_and_S", gwr_ourlabs_l$labs)
gwr_ourlabs_l <- data.frame(labs=c(gwr_ourlabs_l$labs[1:(nrow(gwr_ourlabs_l))], 'Medial_wall'),
                            fi_filter=c(gwr_ourlabs_l$fi_filter[1:(nrow(gwr_ourlabs_l))],0),
                            pcolor=c(gwr_ourlabs_l$pcolor[1:(nrow(gwr_ourlabs_l))],0))

# right
gwr_ourlabs_r <- gwr_ourlabs_rl[!grepl('_left', gwr_ourlabs_rl$labs),]
gwr_ourlabs_r$labs <- gsub("GWR_", "", gwr_ourlabs_r$labs)
gwr_ourlabs_r$labs <- gsub("G.S", "G_and_S", gwr_ourlabs_r$labs)
gwr_ourlabs_r <- data.frame(labs=c(gwr_ourlabs_r$labs[1:(nrow(gwr_ourlabs_r))], 'Medial_wall'),
                            fi_filter=c(gwr_ourlabs_r$fi_filter[1:(nrow(gwr_ourlabs_r))],0),
                            pcolor=c(gwr_ourlabs_r$pcolor[1:(nrow(gwr_ourlabs_r))],0))

# left
gwr_ourlabs_l_test <- gwr_ourlabs_l
colnames(gwr_ourlabs_l_test) <- c('region', 'fi_filter', 'pcolor')
gwr_ourlabs_l_test$region <- sub("_left", "", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub(".1", "", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("rostral", "rostral ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("superior", "superior ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("inferior", "inferior ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("anterior", "anterior ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("posterior", "posterior ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("medial", "medial ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("lateral", "lateral ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("pars", "pars ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("caudal", "caudal ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("middle", "middle ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("isthmus", "isthmus ", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("frontalpole", "frontal pole", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("temporalpole", "temporal pole", gwr_ourlabs_l_test$region)
gwr_ourlabs_l_test$region <- sub("transversetemporal", "transverse temporal", gwr_ourlabs_l_test$region)

gwr_left <- ggseg(gwr_ourlabs_l_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'left',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Gray-White Signal Intensity Ratio')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

# right
gwr_ourlabs_r_test <- gwr_ourlabs_r
colnames(gwr_ourlabs_r_test) <- c('region', 'fi_filter', 'pcolor')
gwr_ourlabs_r_test$region <- sub("_right", "", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub(".1", "", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("rostral", "rostral ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("superior", "superior ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("inferior", "inferior ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("anterior", "anterior ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("posterior", "posterior ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("medial", "medial ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("lateral", "lateral ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("pars", "pars ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("caudal", "caudal ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("middle", "middle ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("isthmus", "isthmus ", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("frontalpole", "frontal pole", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("temporalpole", "temporal pole", gwr_ourlabs_r_test$region)
gwr_ourlabs_r_test$region <- sub("transversetemporal", "transverse temporal", gwr_ourlabs_r_test$region)

gwr_right <- ggseg(gwr_ourlabs_r_test,
      atlas = dk, 
      mapping = aes(fill = fi_filter),
      color = 'black',
      hemisphere = 'right',
      size = .3,
      position = 'stacked') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), na.value = 'gray', limits = c(-.18, .18)) +
  theme_void()+
  ggtitle('Gray-White Signal Intensity Ratio')+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none')

             
            
grid_top <- arrangeGrob(textGrob("Cortical Thickness", vjust = 1.5, gp = gpar(fontsize = 25)), 
                        textGrob("Cortical Surface Area", vjust = 1.5, gp = gpar(fontsize = 25)),
                        ct_left, sa_left,
                        ct_right, sa_right, 
                        ncol=2, 
                        heights = c(.1, 2.5, 2.5))

space <- unit(.5, 'cm')
grid_bottom <- arrangeGrob(textGrob("Gray Matter Volume", vjust = 1.5, gp = gpar(fontsize = 25)), 
                           textGrob("Gray-White Signal Intensity Ratio", vjust = 1.5, gp = gpar(fontsize = 25)),
                           gmv_left, gwr_left, 
                           gmv_right, gwr_right, 
                           ncol = 2, 
                           heights = c(.1, 2.5, 2.5))

# 1000 x 1000
grid.arrange(grid_top,
             textGrob(""),
             grid_bottom,
             heights = c(5, 0.1, 5))

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_2.png",
       grid.arrange(grid_top,
                    textGrob(""),
                    grid_bottom,
                    heights = c(5, 0.1, 5)),
       width = 10.5, height = 10.5, units = 'in', dpi = 500)


# aseg

color_pal <- colfunc(200)
weight_range <- max(dk_gwr_bar_table$fi, na.rm=T)-min(dk_gwr_bar_table$fi, na.rm=T)
interval <- weight_range/200
spaced_vals <- seq(from=-.12, to=.12, by=interval)
dk_gwr_bar_table$color_position <- rep(NA, nrow(dk_gwr_bar_table))
dk_gwr_bar_table$color_from_pos <- rep(NA, nrow(dk_gwr_bar_table))
for (i in (1:nrow(dk_gwr_bar_table))){
  if(!is.na(dk_gwr_bar_table$fi[i])){
    dk_gwr_bar_table$color_position[i] <- sum(dk_gwr_bar_table$fi[i] > spaced_vals)
    dk_gwr_bar_table$color_from_pos[i] <- color_pal[c(sum(dk_gwr_bar_table$fi[i] >= spaced_vals)+1)]
  }
}

aseg_dk_gwr_bar_table <- dk_gwr_bar_table[dk_gwr_bar_table$cat=='ASEG',]

# left
aseg_dk_gwr_bar_table_l <- aseg_dk_gwr_bar_table
aseg_dk_gwr_bar_table_l$region <- aseg_dk_gwr_bar_table$ROI
aseg_dk_gwr_bar_table_l$region <- tolower(aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("left.", "", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("\\.", " ", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("cc", "CC", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("aCCumbens", "accumbens", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("mid_posterior", "mid posterior", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("mid_anterior", "mid anterior", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("ventraldc", "ventral DC", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("white\\.matter", "white matter", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("lat\\.vent", "inferior lateral vent", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("inf", "", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("x5th", "5th", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("x4th", "4th", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("x3rd", "3rd", aseg_dk_gwr_bar_table_l$region)
aseg_dk_gwr_bar_table_l$region <- sub("_", " ", aseg_dk_gwr_bar_table_l$region)

# coronal slice

ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_3_cor_l.png",
       aseg_dk_gwr_bar_table_l %>%
         ggplot() +
         geom_brain(atlas = aseg, 
                    mapping = aes(fill = fi_filter),
                    hemi = 'left') +
         scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), limits = c(-.18, .18)) +
         theme_void()+
         theme(plot.title = element_text(hjust = 0.5, size = 0),
               plot.margin = margin(5, 0, 5, 0)),
       width = 5, height = 10, units = 'in', dpi = 500)


# sagittal slice
ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_3_sag_l.png",
aseg_dk_gwr_bar_table_l %>%
  ggplot() +
  geom_brain(atlas = aseg, 
             mapping = aes(fill = fi_filter)) +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), limits = c(-.18, .18)) +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0)),
width = 5, height = 10, units = 'in', dpi = 500)

# right
aseg_dk_gwr_bar_table_r <- aseg_dk_gwr_bar_table
aseg_dk_gwr_bar_table_r$region <- aseg_dk_gwr_bar_table$ROI
aseg_dk_gwr_bar_table_r$region <- tolower(aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("right.", "", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("\\.", " ", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("cc", "CC", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("aCCumbens", "accumbens", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("mid_posterior", "mid posterior", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("mid_anterior", "mid anterior", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("ventraldc", "ventral DC", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("white\\.matter", "white matter", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("lat\\.vent", "inferior lateral vent", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("inf", "", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("x5th", "5th", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("x4th", "4th", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("x3rd", "3rd", aseg_dk_gwr_bar_table_r$region)
aseg_dk_gwr_bar_table_r$region <- sub("_", " ", aseg_dk_gwr_bar_table_r$region)

# coronal slice
ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_3_cor_r.png",
aseg_dk_gwr_bar_table_r %>%
  ggplot() +
  geom_brain(atlas = aseg, 
             mapping = aes(fill = fi_filter),
             hemi = 'right') +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), limits = c(-.18, .18)) +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none'),
width = 5, height = 10, units = 'in', dpi = 500)


# sagittal slice
ggsave("/Users/ew198/Documents/brainpace/drafts/figure_panels/fig2_3_sag_r.png",
aseg_dk_gwr_bar_table_r %>%
  ggplot() +
  geom_brain(atlas = aseg, 
             mapping = aes(fill = fi_filter)) +
  scale_fill_gradientn(colors = c('#197EC0FF','#71D0F5FF', 'white', '#FD8CC1FF', '#C80813FF'), limits = c(-.18, .18)) +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, size = 0),
        plot.margin = margin(5, 0, 5, 0),
        legend.position = 'none'),
width = 5, height = 10, units = 'in', dpi = 500)



