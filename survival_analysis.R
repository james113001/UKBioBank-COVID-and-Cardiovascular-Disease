library(survival)
library(glmnet)
library(dplyr)
library(data.table)
library(skimr)
library(focus)
library(ggplot2)


############# Dataset should be cleaned - only thing to consider is some of the deaths as covid don't have a test date !!!! -- specdate should be chosen for them ###

setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Data")
factors <- readRDS("CVD_ukbiobank_surv_anal.rds")%>% mutate(eid=as.character(eid)) %>% distinct()


skim(factors)
# setwd("/rds/general/project/hda_21-22/live/TDS/Group_8/Group-Presentation/extraction_and_recording/outcome_definitions_COVID")
# covid_deaths<-readRDS("./Outputs/output_final.rds")%>% mutate(eid=as.character(eid))

setwd("/rds/general/project/hda_21-22/live/TDS/General/Data")
deathcause<-fread("death_cause.txt")
deathcause
deathcause1 <-deathcause %>% filter(cause_icd10=="U071")
skim(deathcause1)

deathdates<-fread("death.txt")
test<- left_join(data.frame(deathcause), data.frame(deathdates), by="eid") 

deathdeets<-test %>%
  #filter(cause_icd10=="U071")%>%
  dplyr::select(-ins_index.x, -arr_index, -level, -ins_index.y, -dsource, -source)%>%
  mutate(eid=as.character(eid))


#covid_deaths %>% filter(!is.na(date_death),incident_case==1)

#used deathdeets instead of covid_deaths
factors %>% select(eid)
class(data.frame(deathdeets %>% select(eid))[1])

covid_survival0<- right_join( data.frame(deathdeets)%>%filter(cause_icd10=="U071"),data.frame(factors), by="eid") %>% distinct()



covid_survival<- covid_survival0 %>%
  mutate(
    cause_icd10=ifelse(is.na(cause_icd10) ,0,cause_icd10),
    cause_icd10=ifelse(cause_icd10=="U071" ,1,0)) %>%
  dplyr::rename(covid_case=cause_icd10)


######lets try to censor all deaths as 1
# covid_survival<- covid_survival0 %>%
#   mutate(cause_icd10=ifelse(cause_icd10!=0,1,0)) %>%
#   dplyr::rename(covid_case=cause_icd10)
######

#remove multiple causes of death in same indiv
covid_survival = covid_survival[order(covid_survival[,'eid'],-covid_survival[,'covid_case']),]
covid_survival = covid_survival[!duplicated(covid_survival$eid),]


#censor living people to latest covid death date, calculate event times
covid_survival_clean<- covid_survival %>%
  dplyr::select(-result, -case) %>%
  mutate(date_of_death=ifelse(is.na(date_of_death)==TRUE, "05/11/2021", date_of_death),
         date_of_death=as.Date(date_of_death, format="%d/%m/%Y"),
         specdate=as.Date(specdate, format="%Y-%m-%d"),
         surv_months = as.numeric(difftime(date_of_death, specdate, units='days'))/30
         )
##three ppl got covid result after dying
covid_survival_clean$surv_months <- ifelse(covid_survival_clean$surv_months < 0, 0, covid_survival_clean$surv_months)
covid_survival_clean[covid_survival_clean$surv_months<0, ]#$surv_months

######ignore
# hmm<-deathdeets%>%
#   filter(cause_icd10=="U071")
# hmm<-hmm[rev(order(as.Date(hmm$date_of_death, format="%d/%m/%Y"))),] #finding the latest covid death
# head(hmm,10) ##latest covid death was nov. 05 2021
# 
# hmm2<-covid_survival_clean[rev(order(as.Date(covid_survival_clean$specdate, format="%Y-%m-%d"))),]
# head(hmm2$specdate,5) ##latest positive test was Oct. 17 2021
             

view<-covid_survival_clean%>% filter(surv_months!=0)
head(view[order(view$surv_months),],10)

finaldata<-covid_survival_clean%>%
  select(-date_of_death, -eid, -specdate) %>%
  filter(surv_months!=0)

covariates<-finaldata %>%
  select(surv_months,covid_case)

factors<-finaldata %>%
  select(-surv_months,-covid_case)




factors<-as.matrix(factors)

#KM plot
DLBCL.surv <- Surv(covariates$surv_months, covariates$covid_case)
DLBCL.fit <- survfit(DLBCL.surv ~ 1)
plot(DLBCL.fit, main = "Kaplan-Meier estimate of the survival function",
     xlab = "Time (months)", ylab = "Survival function", xlim = c(0, 20))


set.seed(1551)
fold <- sample(rep(seq(5), length = nrow(factors)))

cox_lasso <- cv.glmnet(factors, DLBCL.surv, family = "cox",
                       foldid = fold, standardize = TRUE)


best_lambda <- cox_lasso$lambda.min

#Coeffs
coef(cox_lasso, s = "lambda.min")

#Variable selection
out = VariableSelection(x = factors,seed=1,y = DLBCL.surv, family = "cox")
pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_I_Analysis/Cali_plot_1.pdf"))

CalibrationPlot(out)
dev.off()

selprop=SelectionProportions(out)
print(selprop)

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)

par(mar=c(10,5,1,1))
pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_I_Analysis/Selection_proportion.pdf"))


plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)

abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}
dev.off()

out = GraphicalModel(xdata=factors, verbose=FALSE)

pdf(file = paste0("/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_I_Analysis/Caliplot2.pdf"))
CalibrationPlot(out)
dev.off()


#166.5
pdf(file="/rds/general/project/hda_21-22/live/TDS/Group_8/Comp_Epi/Part_I_Analysis/plots/lambda_plots.pdf")
plot(cox_lasso$glmnet.fit, xvar = "lambda")
dev.off()

summary(cox_lasso)
plot(cox_lasso)
