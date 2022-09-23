### This (mech shop) sets up for sensitivity analysis for lower VE against severe outcomes in older adults including.
### This reduction in VE against severe outcomes includes: (1) a faster speed of waning, and (2) lower strength of inital protection.
### Creates:SA_VE_older_muted_SO


require(ggpubr); require(readr); require(gridExtra); require(ggplot2); require(tidyverse);


### SETUP
rm(list=ls())
plotting_standard =  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = 'black'))


### (1) SPEED OF WANING ##############################################################################################################################################
raw <- read.csv(file = '1_inputs/VE_severe_outcomes.csv',header=TRUE)

predicted_distribution = data.frame()
subplot_list = list()
for (i in 1:length(unique(raw$age_group))){
  for (d in 1: length(unique(raw$dose))){
    this_age = unique(raw$age_group)[i]
    this_dose = unique(raw$dose)[d] 
    
    workshop_real = raw[raw$dose == this_dose & raw$age_group == this_age,]

    attach(workshop_real)
    #model = lm(log(VE)~days)
    model = lm(VE~days)
    detach(workshop_real)
    
    time <- seq(0, 365)
    #workshop_predicted <- exp(predict(model,list(days=time)))
    workshop_predicted <- predict(model,list(days=time))
    workshop_predicted = data.frame(cbind(days = time, VE = workshop_predicted))
    
    workshop_predicted = workshop_predicted %>%
      mutate(dose = unique(raw$dose)[d],
             age_group = unique(raw$age_group)[i],
             rsquared = model_rsquared )
    
    predicted_distribution = rbind(predicted_distribution,workshop_predicted)
  }
}

plot_dose2 = ggplot() + 
  geom_line(data=predicted_distribution[predicted_distribution$dose == 2,],aes(x=days,y=VE,color=as.factor(age_group))) +
  geom_point(data = raw[raw$dose == 2,], aes(x=days,y=VE,color=as.factor(age_group))) +
  ylim(0,1)+
  labs(color='age group') +
  xlab('days since vaccination') +
  ggtitle('primary schedule') +
  plotting_standard
plot_dose3 = ggplot() + 
  geom_line(data=predicted_distribution[predicted_distribution$dose == 3,],aes(x=days,y=VE,color=as.factor(age_group))) +
  geom_point(data = raw[raw$dose == 3,], aes(x=days,y=VE,color=as.factor(age_group))) +
  ylim(0,1)+
  labs(color='age group') +
  xlab('days since vaccination') +
  ggtitle('booster dose') +
  plotting_standard
grid.arrange(plot_dose2,plot_dose3, nrow=2)


### Internal VE
apply_distribution <- predicted_distribution %>%
  group_by(dose,age_group) %>%
  mutate(VE_internal = VE / max(VE)) %>%
  ungroup() %>%
  select(dose,age_group,days,VE_internal,VE)

plot_dose2 = ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$dose == 2,], aes(x=days,y=VE_internal,color=as.factor(age_group))) +
  ylim(0,1)  +
  ylab('% of max protection') +
  xlab('days since vaccination') +
  labs(color='age group') +
  ggtitle('primary schedule') +
  plotting_standard
plot_dose3 = ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$dose == 3,], aes(x=days,y=VE_internal,color=as.factor(age_group))) +
  ylim(0,1) +
  ylab('% of max protection') +
  xlab('days since vaccination') +
  labs(color='age group') +
  ggtitle('booster dose') +
  plotting_standard
grid.arrange(plot_dose2,plot_dose3, nrow=2)
#_________________________________________________________________________________________________________________________________________



### (2) RATIO ##############################################################################################################################################
workshop_overall = raw %>% 
  filter(age_group == 'overall') %>% 
  rename(VE_overall = VE) %>%
  select(dose,days,VE_overall)
workshop_age = raw %>% 
  filter(! age_group == 'overall') %>%
  left_join(workshop_overall, by = c('dose','days')) %>%
  mutate(VE_overall = VE/VE_overall) 

apply_ratio = workshop_age %>% filter(days == 22) %>% rename(VE_ratio = VE_overall,agegroup_RAW = age_group) %>% select(dose,agegroup_RAW,VE_ratio)
#_________________________________________________________________________________________________________________________________________



### COVERT TO MODEL AGE GROUP ##################################################################################################################################
CS_age_groupings = c(0,59,79,110) #age groupings in VE estimate data
pop_RAW <- pop_setting_orig %>%
  mutate(agegroup_RAW = cut(age,breaks = CS_age_groupings, include.lowest = T, labels = unique(apply_ratio$agegroup_RAW)),
         agegroup_MODEL = cut(age,breaks = age_groups_num, include.lowest = T, labels = age_group_labels)) %>%
  ungroup() %>%
  group_by(agegroup_MODEL) %>%
  mutate(model_group_percent = population/sum(population))

#ratio
workshop = pop_RAW %>% left_join(apply_ratio) %>% 
  mutate(interim = model_group_percent * VE_ratio)
workshop = aggregate(workshop$interim, by=list(category = workshop$dose,workshop$agegroup_MODEL),FUN=sum)
colnames(workshop) = c('dose','age_group','VE_ratio')

apply_ratio_MODEL = workshop %>% arrange(dose) %>%
  mutate(schedule = case_when(
    dose > 2 ~ 'booster',
    TRUE ~ 'primary'
  ))  %>%  
  select(-dose)

#waning
workshop = apply_distribution %>% rename(agegroup_RAW = age_group)

workshop = pop_RAW %>% left_join(workshop) %>% 
  mutate(interim = model_group_percent * VE_internal)
workshop = aggregate(workshop$interim, by=list(category = workshop$dose,workshop$agegroup_MODEL,workshop$days),FUN=sum)
colnames(workshop) = c('dose','age_group','days','VE_internal')

apply_distribution_MODEL = workshop %>% arrange(dose)  %>%
  mutate(schedule = case_when(
    dose > 2 ~ 'booster',
    TRUE ~ 'primary'
  )) %>%  
  select(-dose)

plot_dose2 = ggplot() + 
  geom_line(data = apply_distribution_MODEL[apply_distribution_MODEL$schedule == 'primary',], aes(x=days,y=VE_internal,color=as.factor(age_group))) +
  ylim(0,1)  +
  ylab('% of max protection') +
  xlab('days since vaccination') +
  labs(color='age group') +
  ggtitle('primary schedule') +
  plotting_standard
plot_dose3 = ggplot() + 
  geom_line(data = apply_distribution_MODEL[apply_distribution_MODEL$schedule == 'booster',], aes(x=days,y=VE_internal,color=as.factor(age_group))) +
  ylim(0,1) +
  ylab('% of max protection') +
  xlab('days since vaccination') +
  labs(color='age group') +
  ggtitle('booster dose') +
  plotting_standard
grid.arrange(plot_dose2,plot_dose3, nrow=2)
#_________________________________________________________________________________________________________________________________________



### APPLY TO POINT ESTIMATES ##################################################################################################################################
load(file = "1_inputs/VE_estimates_imputed.Rdata")
load(file = "1_inputs/VE_booster_estimates.Rdata")

point_estimates = VE_estimates_imputed %>% 
  filter(outcome_family == 'severe_outcome' & !(vaccine_type == 'Pfizer' & dose == 3)) %>%
  select(strain,vaccine_type,dose,outcome,outcome_family,VE) %>%
  mutate(schedule = case_when(
    dose > 2 ~ 'booster',
    dose == 2 & vaccine_type == "Johnson & Johnson" ~ 'booster',
    TRUE ~ 'primary'
  ))

point_estimates_booster = VE_booster_estimates %>% 
  filter(outcome_family == 'severe_outcome') %>%
  select(strain,vaccine_type,primary_if_booster,dose,outcome,outcome_family,VE) %>%
  mutate(schedule = 'booster')

point_estimates = bind_rows(point_estimates,point_estimates_booster)

together = point_estimates %>% 
  left_join(apply_ratio_MODEL,by='schedule') %>%
  left_join(apply_distribution_MODEL, by = c('schedule','age_group')) %>%
  rename(VE_days = VE) %>%
  mutate(VE_days = VE_days*VE_internal*VE_ratio/100)
#_________________________________________________________________________________________________________________________________________



### PLOT ##################################################################################################################################
if (exists("vax_type_list") == FALSE){  vax_type_list = c("AstraZeneca","Johnson & Johnson", "Pfizer", "Sinopharm" ) }

waning_to_plot = together %>%
  filter(vaccine_type %in% vax_type_list) %>%
  mutate(immunity = paste(vaccine_type,dose))

strain_test = 'omicron'
outcome_test = 'severe_disease'
vaccine_type_test = 'Johnson & Johnson'
vaccine_type_test = 'Pfizer'

ggplot() +
  geom_point(data=waning_to_plot[waning_to_plot$strain == strain_test  & waning_to_plot$outcome == outcome_test & waning_to_plot$vaccine_type == vaccine_type_test,],
            aes(x=days,y=VE_days,color=as.factor(age_group),shape=as.factor(dose)),na.rm=TRUE) +
  labs(title=(paste("Waning of VE against","(",strain_test,")"))) +
  xlab("days since vaccination") +
  ylab("% max protection") +
  ylim(0,1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

waning = together %>% mutate(waning = TRUE)
no_waning = together %>% mutate(waning = FALSE) %>%
  group_by(strain,vaccine_type,primary_if_booster,dose,age_group,outcome) %>%
  mutate(VE_days = max(VE_days))

### SAVE ##################################################################################################################################
SA_VE_older_muted_SO = rbind(waning,no_waning) %>% select(strain, vaccine_type,primary_if_booster, dose, outcome,days,age_group,VE_days,waning)
save(SA_VE_older_muted_SO, file = '1_inputs/SA_VE_older_muted_SO.Rdata')
#_________________________________________________________________________________________________________________________________________