rm(list=ls())
library(gtsummary)
library(knitr)
library(kableExtra)
library(tidyverse)
library(focus)
library(withr)
library(gridExtra)
library(factoextra)
library(fastDummies)
library(survminer) 
library(survival) 
library(cowplot)
library(ggpubr)
library(survival)
library(ranger)
library(ggplot2)
library(ggfortify)
library(glmnet)



#load data and change path to folder in pdf() folders ### change family!!

x <- data.matrix(xpred) 
y <- data.matrix(Surv(df2$surv_days,df2$incident_case))

#Turn x into dummies

#head(x %>%  select_if(is.factor))
#dum <- dummy_cols(x %>%  select(covid), select_columns = colnames(x %>%  select(covid)))%>%select(-covid,-covid_No)
#x<-cbind(x%>%select(-covid,-incident_case),dum)

x_train <- model.matrix( ~ .-1, x)

fold <- sample(rep(seq(5), length = nrow(x)))
glm<-cv.glmnet (x_train, y,family="cox",foldid = fold ,alpha = 1)


pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_II_Analysis/plots_nested_analysis/glm_lasso.pdf"))
plot(glm)
dev.off()


beta_lasso = coef(glm, s="lambda.min")[2:(ncol(x_train)+1),]
selected_lasso = names(beta_lasso)[which(beta_lasso!=0)]
print(paste0(length(selected_lasso), " selected variables"))
print(selected_lasso)

coefs_lass=exp(coef(glm, s="lambda.min")[2:(ncol(x_train)+1),][which(beta_lasso!=0)])

# See if there are any filters applied messing with results?!!

out = VariableSelection(x = x_train,seed=1,y = y,family = "cox")

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_II_Analysis/plots_nested_analysis/Cali_plot_1.pdf"))
CalibrationPlot(out)
dev.off()

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

par(mar=c(10,5,1,1))
pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_II_Analysis/plots_nested_analysis/Selection_proportion.pdf"))

plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)

abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}
dev.off()



out = GraphicalModel(xdata=x_train, verbose=FALSE)

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_II_Analysis/plots_nested_analysis/Cali_plot_2.pdf"))
CalibrationPlot(out)
dev.off()