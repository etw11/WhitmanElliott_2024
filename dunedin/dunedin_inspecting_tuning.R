# inspecting parameter tuning results
# 9/11/23
# ethan whitman

library(superheat)

# load in all iterations
file_list <- list.files('/Users/ew198/Documents/brainpace/results/tuning_parameters_dk_gwr')
file_dir <- '/Users/ew198/Documents/brainpace/results/tuning_parameters_dk_gwr/'
pred_alpha_all <- data.frame(alpha=numeric(), lambda=numeric(), RMSE=numeric(), Rsquared=numeric(), MAE=numeric(), 
                             RMSED=numeric(), RsquaredSD=numeric(), MAESD=numeric(), r=numeric(), items=numeric())
best_tune_all <- data.frame(alpha=numeric(), lambda=numeric(), items=numeric())

for (f in file_list){
  load(paste0(file_dir, f))
  
  pred_alpha_all <- rbind(pred_alpha_all, cbind(as.data.frame(enet$results),  
                                                      items=rep(sum(as.numeric(coef(enet$finalModel, enet$bestTune$lambda)) !=0 ), nrow(as.data.frame(enet$results))) ))
  best_tune_all <- rbind(best_tune_all, data.frame(alpha=enet$bestTune$alpha, 
                                                   lambda=enet$bestTune$lambda,
                                                   items=sum(as.numeric(coef(enet$finalModel, enet$bestTune$lambda)) !=0 ) ) )
}

# matrix of lambda x alpha values

lambdas <- 10^seq(-3, 3, length = 25)
alphas <- seq(0,1,length = 15)

# tuning_parameters_dkt_nogwr_r2_mat <- matrix(NA,25,15)
# tuning_parameters_dkt_nogwr_mae_mat <- matrix(NA,25,15)
tuning_parameters_dk_gwr_r2_mat <- matrix(NA,25,15)
tuning_parameters_dk_gwr_mae_mat <- matrix(NA,25,15)

i <- 1
for(l in lambdas){
  j <- 1
  for (a in alphas){
    temp_grid <- pred_alpha_all[as.logical(c(pred_alpha_all$alpha==a)*c(pred_alpha_all$lambda==l)), ]
    tuning_parameters_dk_gwr_r2_mat[i,j] <- mean(temp_grid$Rsquared)
    tuning_parameters_dk_gwr_mae_mat[i,j] <- mean(temp_grid$MAE)
    j <- j+1
  }
  i <- i+1
}

colnames(tuning_parameters_dk_gwr_r2_mat) <- round(alphas, 3)
rownames(tuning_parameters_dk_gwr_r2_mat) <- round(lambdas, 3)
colnames(tuning_parameters_dk_gwr_mae_mat) <- round(alphas, 3)
rownames(tuning_parameters_dk_gwr_mae_mat) <- round(lambdas, 3)

tuning_parameters_dk_gwr_r2_mat[is.na(tuning_parameters_dk_gwr_r2_mat)] <- 0

superheat(sqrt(tuning_parameters_dk_gwr_r2_mat),
          heat.pal = c("red", "white", "lightblue"),
          heat.pal.values = c(0,.9, 1),
          bottom.label.text.size = 4,
          bottom.label.text.angle = 45,
          left.label.text.size = 4,
          legend.height = 0.075,
          legend.text.size = 10,
          left.label.col = "white", bottom.label.col = "white",
          X.text = round(sqrt(tuning_parameters_dk_gwr_r2_mat), 2)
)


superheat(tuning_parameters_dk_gwr_mae_mat,
          heat.pal = c("lightblue", "white", "red"),
          heat.pal.values = c(0, 0.5, 1),
          bottom.label.text.size = 4,
          bottom.label.text.angle = 45,
          left.label.text.size = 4,
          legend.height = 0.075,
          legend.text.size = 10,
          left.label.col = "white", bottom.label.col = "white",
          X.text = round(tuning_parameters_dk_gwr_mae_mat, 2)
)



# looking at most common selected pairs

selected <- unique(best_tune_all[!duplicated(best_tune_all),1:2])

sums <- c(rep(0,nrow(selected)))
for (i in 1:nrow(selected)){
  s <- 0
  for (b in 1:nrow(best_tune_all)){
    if (sum(selected[i,] == best_tune_all[b,1:2]) == 2){
      s <- s + 1
    }
  }
  sums[i] <- s
}
selected_sums <- cbind(selected, sums)

sum(as.character(best_tune_all$alpha) == "0.785714285714286")

Mode(best_tune_all$alpha)
Mode(best_tune_all$lambda)
Mode(best_tune_all$items)

# distribution of lambda selections
best_lambdas <- ggplot(data = best_tune_all, aes(x = lambda)) + 
  geom_histogram() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15))

# distribution of alpha selections
best_alphas <- ggplot(data = best_tune_all, aes(x = alpha)) + 
  geom_histogram() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15))

# distribution of N items
best_items <- ggplot(data = best_tune_all, aes(x = items)) + 
  geom_histogram() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 20),
        legend.text = element_text(size = 15))

grid.arrange(best_lambdas, best_alphas, best_items, nrow = 1)

