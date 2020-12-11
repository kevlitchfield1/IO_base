library(pheatmap)
library("CombinePValue")
library(ggpubr)
library(ggplot2)
library(cowplot)
library(nlme)
library(plyr)
require(reshape2)
library(gridExtra)
library(cowplot)
library(scales)
library('dplyr')
library(ggbeeswarm)
library(meta)
library("glmnet")
library("ROCR")
library("xgboost")
library("DiagrammeR")
library("pROC")

# N=1000, 500 discovery, 500 validation, frequency of severe tox=25%, non sever tox=75%

# Proteomics data based on Luminex platform preliminary for 25 markers
# Difference between groups based on an effect size (minimum 1.15) seen in preliminary data from Sophie Papa lab
proteomics_tox<-data.frame(matrix(1.15*sample(10:2000,(250*25),replace=T),nrow=250))
proteomics_no_tox<-data.frame(matrix(sample(10:2000,(750*25),replace=T),nrow=750))
proteomics<-rbind(proteomics_tox,proteomics_no_tox)
proteomics$proteomics_signature<-rowMeans(proteomics)

# Standardised PRS, scale -3 to +3
# PRS difference between groups estimated from a clinical trial of 930 patients where we had information on any grade 3/grade 4 toxicity following treatment with Capecitabine +/- Bevacizumab
genotypes_tox<-data.frame(matrix(rnorm(250)+0.05, nrow=250));names(genotypes_tox)[1]<-"genetic_PRS"
genotypes_no_tox<-data.frame(matrix(rnorm(750), nrow=750));names(genotypes_no_tox)[1]<-"genetic_PRS"
genotypes<-rbind(genotypes_tox,genotypes_no_tox)

# Assume panel of 500 autoantibodies, scale 0 to 100
# Difference in groups based on minimum effect size of 1.03, taken from Gowen et al. 2018
serology_tox<-data.frame(matrix(1.03*sample(1:100,(250*500),replace=T),nrow=250))
serology_no_tox<-data.frame(matrix(sample(1:100,(750*500),replace=T),nrow=750))
serology<-rbind(serology_tox,serology_no_tox)
serology$serology_signature<-rowMeans(serology)

# labels
labels<-data.frame(c(rep(1,250),rep(0,750)));names(labels)[1]<-"label"

# Merge data
all_dat<-cbind(proteomics$proteomics_signature,genotypes,serology$serology_signature,labels)
train_index <- sample(seq_len(nrow(all_dat)), size = (nrow(all_dat)*0.5))
train_data<-all_dat[train_index,];test_data<-all_dat[-train_index,]

# Hyper parameters from grid search using caret package
mva_mod <- xgboost(data = as.matrix(train_data[,1:3]), label = train_data$label,objective = "binary:logistic",verbose=F,max.depth = 4,nrounds = 20,eta=0.2,"eval_metric" = "logloss")
importance_matrix <- xgb.importance(model = mva_mod)
xgb.plot.importance(importance_matrix = importance_matrix,left_margin = 18)

# Plot ROC curve
pred <- predict(mva_mod, as.matrix(test_data[,1:3]))
pr <- prediction(as.numeric(pred),test_data$label)
auc <- performance(pr, measure = "auc")
auc_val<-auc@y.values[[1]]
auc <- performance(pr,"tpr","fpr");auc2 <- performance(pr2,"tpr","fpr")
plot(auc,col="darkblue",lwd=3,xaxt="n")+axis(2, at = seq(0.0,1.0, by = 0.2), las=2)+mtext(paste0("Multivariate AUC: ",round(auc_val,2)),col="darkblue",line=-8,cex=0.8)+ abline(coef = c(0,1))

