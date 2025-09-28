medicare_loader = function(resp)
{
library('haven')
setwd("/Users/cgmeixide/Projects/implied_interventions/oregon_Data")

descriptive <- read_dta("oregonhie_descriptive_vars.dta")
#inperson    <- read_dta("oregonhie_inperson_vars.dta")
survey0m    <- read_dta("oregonhie_survey0m_vars.dta")
survey6m    <- read_dta("oregonhie_survey6m_vars.dta")
survey12m   <- read_dta("oregonhie_survey12m_vars.dta")
#ed          <- read_dta("oregonhie_ed_vars.dta")
#patterns    <- read_dta("oregonhie_patterns_vars.dta")
#stateprogs  <- read_dta("oregonhie_stateprograms_vars.dta")

if(FALSE) {
str(descriptive)
str(inperson)
str(survey0m)
str(survey6m)
str(survey12m)
str(ed)
str(patterns)
str(stateprogs)

names(descriptive)
names(inperson)
names(survey0m)
names(survey6m)
names(survey12m)
names(ed)
names(patterns)
names(stateprogs)

names(descriptive)

# descriptve, survey0m, survey6m

names(descriptive)
names(survey0m)
names(survey6m)

}
responses=c('needmet_med_12m','cost_tot_oop_12m', 'happiness_12m', 'doc_any_12m', 'er_any_12m','health_gen_12m','dep_dx_12m')         # Any doctor or ER visit at 6 months     # Depression score or self-reported happiness

covariatedes=c('zip_msa_list','birthyear_list','female_list', 'english_list')

covariatesur=c('dia_dx_0m','dep_dx_0m')#c('edu_0m', 'race_white_0m', 'race_black_0m', 'race_hisp_0m', 'hhinc_cat_0m', 'employ_0m')

treatment='approved_app'

instrument='treatment'

W <- as.matrix(descriptive[ , covariatedes])
W = cbind(W, as.matrix(survey0m[ , covariatesur]))

sum(complete.cases(W))

A=as.matrix(descriptive[ , treatment])
Z=as.matrix(descriptive[ , instrument])
A[is.na(A)]=0
for(i in 1:length(responses)) {
Y=as.matrix(survey12m[ , responses[i]])
#print(sum(complete.cases(Y)))
}

Y=as.matrix(survey12m[ , responses[resp]])

all=data.frame(Y=Y,A=A,Z=Z,W=W)

data=all[complete.cases(all),]

miny=min(data$W.birthyear_list)

maxy=max(data$W.birthyear_list)

data$W.birthyear_list = (data$W.birthyear_list - miny)/(maxy - miny)
sum(is.na(A))
return(data)

}




