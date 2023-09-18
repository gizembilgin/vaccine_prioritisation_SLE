### This (mech shop) applies a waning distribution to the point estimates of VE against severe outcomes.
### Creates: VE_waning_distribution_SO

require(ggpubr); require(readr); require(gridExtra); require(ggplot2); require(tidyverse);
#rm(list=ls())


###(1/3) Load distribution
raw <- read.csv(file = '1_inputs/VE_severe_outcomes.csv',header=TRUE)

### INITIAL VISULALISATION 
plot1 = ggplot(data = raw[raw$dose == 2,]) + 
  ggtitle('Primary schedule') +
  geom_point(aes(x=days,y=VE,color=as.factor(age_group))) + 
  geom_errorbar(aes(x=days, ymin=LB, ymax=UB))+
  ylim(0,1)
plot2 = ggplot(data = raw[raw$dose == 3,]) + 
  ggtitle('Booster') +
  geom_point(aes(x=days,y=VE,color=as.factor(age_group))) + 
  geom_errorbar(aes(x=days, ymin=LB, ymax=UB)) +
  ylim(0,1)
grid.arrange(plot1, plot2, nrow=2)


predicted_distribution = data.frame()
plot_list = list()
for (j in 1:length(unique(raw$age_group))){
  subplot_list = list()
  for (i in 1: length(unique(raw$dose))){
    workshop_real = raw[raw$dose == unique(raw$dose)[i] & 
                                     raw$age_group == unique(raw$age_group)[j],]
    attach(workshop_real)
    model = lm(VE~days)
    
    #summary(model)
    model_rsquared = summary(model)$adj.r.squared
    detach(workshop_real)
    
    time <- seq(0, 365)
    workshop_predicted <- predict(model,list(days=time))
    workshop_predicted = data.frame(cbind(days = time, VE = workshop_predicted))
    
    subplot_list[[i]] = ggplot() + 
      geom_point(data=workshop_real,aes(x=days,y=VE)) +
      geom_line(data = workshop_predicted, aes(x=days,y=VE)) +
      ylim(0,1)+
      ggtitle(paste(unique(raw$dose)[i],unique(raw$age_group)[j],round(model_rsquared,2)))
    
    workshop_predicted = workshop_predicted %>%
      mutate(age_group = unique(raw$age_group)[j],
             dose = unique(raw$dose)[i],
             rsquared = model_rsquared )
    
    predicted_distribution = rbind(predicted_distribution,workshop_predicted)
  }
  plot_list[[j]] = grid.arrange(subplot_list[[1]], subplot_list[[2]], nrow=2)
}

predicted_distribution = predicted_distribution %>% mutate(plot_label = paste(age_group,"(R squared",round(rsquared,digits=3)))

ggplot() + 
  geom_point(data = raw[raw$dose == 2,], aes(x=days,y=VE,color=as.factor(age_group))) +
  geom_line(data = predicted_distribution[predicted_distribution$dose == 2,], aes(x=days,y=VE,color=as.factor(plot_label))) +
  ylim(0,1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'))
ggplot() + 
  geom_point(data = raw[raw$dose == 3,], aes(x=days,y=VE,color=as.factor(age_group))) +
  geom_line(data = predicted_distribution[predicted_distribution$dose == 3,], aes(x=days,y=VE,color=as.factor(age_group))) +
  ylim(0,1) +
  theme_bw() + 
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'))

plot1 = ggplot(data = raw[raw$dose == 2 & raw$age_group == 'overall',]) + 
  ggtitle('Primary schedule') +
  geom_point(aes(x=days,y=VE)) + 
  geom_errorbar(aes(x=days, ymin=LB, ymax=UB))+
  geom_line(data = predicted_distribution[predicted_distribution$dose == 2 & predicted_distribution$age_group == 'overall',], 
            aes(x=days,y=VE))+
  ylim(0,1)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'))
plot2 = ggplot(data = raw[raw$dose == 3 & raw$age_group == 'overall',]) + 
  ggtitle('Booster') +
  geom_point(aes(x=days,y=VE)) + 
  geom_errorbar(aes(x=days, ymin=LB, ymax=UB)) +
  geom_line(data = predicted_distribution[predicted_distribution$dose == 3 & predicted_distribution$age_group == 'overall',], 
            aes(x=days,y=VE)) +
  ylim(0,1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'))
grid.arrange(plot1, plot2, nrow=2)


ggplot(data = raw[raw$age_group == 'overall',]) + 
  geom_point(aes(x=days,y=VE,color=as.factor(dose))) + 
  geom_errorbar(aes(x=days, ymin=LB, ymax=UB)) +
  geom_line(data = predicted_distribution[ predicted_distribution$age_group == 'overall',], 
            aes(x=days,y=VE,color=as.factor(dose))) +
  ylim(0,1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'))

predicted_distribution_SO = predicted_distribution
save(predicted_distribution_SO, file = '1_inputs/VE_predicted_distribution_SO.Rdata')



###(2/3) Apply predicted distributions
#(A) calculate internal
apply_distribution <- predicted_distribution_SO %>%
  filter(age_group == 'overall') %>% #COMEBACK SENSITIVITY
  mutate(outcome_family =  'severe_outcome') %>%
  group_by(outcome_family,dose) %>%
  mutate(VE_internal = VE / max(VE)) %>%
  ungroup() %>%
  select(dose,outcome_family,days,VE_internal)
ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$dose == 2,], aes(x=days,y=VE_internal)) +
  ylim(0,1)
ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$dose == 3,], aes(x=days,y=VE_internal)) +
  ylim(0,1)

#copy dose 2 for dose 1
# dose_1 = apply_distribution %>% filter(dose == 2) %>% mutate(dose = 1)
# apply_distribution = rbind(apply_distribution,dose_1)
apply_distribution = apply_distribution %>%
  mutate(schedule = case_when(
    dose > 2  ~ 'booster',
    TRUE ~ 'primary'
  )) %>% 
  select(-dose)
#columns = outcome_family, days, VE_Internal, schedule


#(B) apply to point estimates
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
  left_join(apply_distribution, by = c('outcome_family','schedule')) %>%
  rename(VE_days = VE) %>%
  mutate(VE_days = VE_days*VE_internal/100)



###(3/3) Plot distributions and save VE_waning_distribution
#(A/B) Plot
if (exists("vax_type_list") == FALSE){  vax_type_list = c("AstraZeneca","Johnson & Johnson", "Pfizer", "Sinopharm" ) }

waning_to_plot = together %>%
  filter(vaccine_type %in% vax_type_list) %>%
  mutate(immunity = paste(vaccine_type,dose))

strain_test = 'omicron'
outcome_test = 'severe_disease'

ggplot() +
  geom_line(data=waning_to_plot[waning_to_plot$strain == strain_test  & waning_to_plot$outcome == outcome_test,],
            aes(x=days,y=VE_days,color=as.factor(immunity)),na.rm=TRUE) +
  labs(title=(paste("Waning of VE against","(",strain_test,")"))) +
  xlab("days since vaccination") +
  ylab("% max protection") +
  ylim(0,1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

waning = together %>% mutate(waning = TRUE)
no_waning = together %>% mutate(waning = FALSE) %>%
  group_by(strain,vaccine_type,primary_if_booster,dose,outcome) %>%
  mutate(VE_days = max(VE_days))


VE_waning_distribution_SO = rbind(waning,no_waning) %>% select(strain, vaccine_type,primary_if_booster, dose, outcome,days,VE_days,waning)
save(VE_waning_distribution_SO, file = '1_inputs/VE_waning_distribution_SO.Rdata')


