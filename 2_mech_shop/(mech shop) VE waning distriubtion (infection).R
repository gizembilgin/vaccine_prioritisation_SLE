### This (mech shop) applies a waning distribution to the point estimates of VE against acquisition
### (VE_estimates_imputed created in `(mech shop) VE primary doses point estimates.R`).
### Creates: VE_waning_distribution


require(ggpubr); require(readr); require(gridExtra); require(ggplot2); require(tidyverse);
rm(list=ls())


###(1/3) Load distribution
andrews <- read.csv(file = '1_inputs/VE_Andrews_shape.csv',header=TRUE)

### INITIAL VISULALISATION 
dose_number = 2
plot1 = ggplot(data = andrews[andrews$strain == 'delta' & andrews$dose == dose_number,]) + 
  ggtitle(paste('Delta')) +
  geom_point(aes(x=week,y=VE,color=as.factor(vaccine_type))) + 
  geom_errorbar(aes(x=week, ymin=VE_lower, ymax=VE_upper))
plot2 = ggplot(data = andrews[andrews$strain == 'omicron' & andrews$dose == dose_number,]) + 
  ggtitle(paste('Omicron')) +
  geom_point(aes(x=week,y=VE,color=as.factor(vaccine_type))) + 
  geom_errorbar(aes(x=week, ymin=VE_lower, ymax=VE_upper))
grid.arrange(plot1, plot2, nrow=2)


andrews <- andrews %>%  filter(dose<3)
andrews_to_fit = andrews
andrews_to_fit$days = andrews_to_fit$week * 7

predicted_distribution = data.frame()
plot_list = list()
for (j in 1:length(unique(andrews_to_fit$vaccine_type))){
  subplot_list = list()
  for (i in 1: length(unique(andrews_to_fit$strain))){
    workshop_real = andrews_to_fit[andrews_to_fit$strain == unique(andrews_to_fit$strain)[i] & 
                                     andrews_to_fit$vaccine_type == unique(andrews_to_fit$vaccine_type)[j],]
    attach(workshop_real)
    model = lm(log(VE)~days)
    model_rsquared = summary(model)$adj.r.squared
    detach(workshop_real)
    
    time <- seq(0, 365)
    workshop_predicted <- exp(predict(model,list(days=time)))
    workshop_predicted = data.frame(cbind(days = time, VE = workshop_predicted))
    
    subplot_list[[i]] = ggplot() + 
      geom_point(data=workshop_real,aes(x=days,y=VE)) +
      geom_line(data = workshop_predicted, aes(x=days,y=VE)) +
      ylim(0,100)+
      ggtitle(paste(unique(andrews_to_fit$strain)[i],unique(andrews_to_fit$vaccine_type)[j],round(model_rsquared,2)))
    
    workshop_predicted = workshop_predicted %>%
      mutate(vaccine_type = unique(andrews_to_fit$vaccine_type)[j],
             strain = unique(andrews_to_fit$strain)[i],
             rsquared = model_rsquared )
    
    predicted_distribution = rbind(predicted_distribution,workshop_predicted)
  }
  plot_list[[j]] = grid.arrange(subplot_list[[1]], subplot_list[[2]], nrow=2)
}

predicted_distribution = predicted_distribution %>% mutate(plot_label = paste(vaccine_type,"(R squared",round(rsquared,digits=3)))

ggplot() + 
  geom_point(data = andrews_to_fit[andrews_to_fit$strain == 'delta',], aes(x=days,y=VE,color=as.factor(vaccine_type))) +
  geom_line(data = predicted_distribution[predicted_distribution$strain == 'delta',], aes(x=days,y=VE,color=as.factor(plot_label))) +
  ylim(0,100) +
  labs(color='vaccine type') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))
ggplot() + 
  geom_point(data = andrews_to_fit[andrews_to_fit$strain == 'omicron',], aes(x=days,y=VE,color=as.factor(vaccine_type))) +
  geom_line(data = predicted_distribution[predicted_distribution$strain == 'omicron',], aes(x=days,y=VE,color=as.factor(vaccine_type))) +
  ylim(0,100) +
  theme_bw() + 
  labs(color='vaccine type') +
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))

save(predicted_distribution, file = '1_inputs/VE_predicted_distribution.Rdata')



###(2/3) Apply predicted distributions
#(A) calculate internal
apply_distribution <- predicted_distribution %>%
  mutate(vaccine_mode = case_when(
    vaccine_type %in% c('Pfizer','Moderna') ~ 'mRNA',
    vaccine_type == 'AstraZeneca' ~ 'viral'),
    outcome_family =  'acquisition') %>%
  group_by(outcome_family,strain,vaccine_type) %>%
  mutate(VE_internal = VE / max(VE)) %>%
  ungroup() %>%
  select(strain,vaccine_type,vaccine_mode,days,VE_internal,VE)
ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$strain == 'delta',], aes(x=days,y=VE_internal,color=as.factor(vaccine_type))) +
  ylim(0,1)
ggplot() + 
  geom_line(data = apply_distribution[apply_distribution$strain == 'omicron',], aes(x=days,y=VE_internal,color=as.factor(vaccine_type))) +
  ylim(0,1)


#(B) take direct where available
direct = apply_distribution %>%
  mutate(dose = 2) %>%
  select(-VE_internal)

#(C) ratio apply distribution
# > unique(VE_estimates_imputed$vaccine_type)
# [1] "AstraZeneca"       "Johnson & Johnson" "Moderna"           "Pfizer"            "Sinopharm"         "Sinovac"     
# impute:
# (1) dose 1 from existing dose 2 (AZ, Moderna, Pfizer)
#         Note: this may not be true, but also not relevant for a long period
# (2) dose 1 and 2 from 'closest' vaccine (mode = viral) i.e. AZ for J&J, Sinopharm and Sinovac
#         Note: possible issue with J&J and later boosters
D = 2
load(file = "1_inputs/VE_estimates_imputed.Rdata")

imputed = data.frame()
for (s in c('omicron','delta')){
  for (d in 1:D){
    for (t in 1:length(unique(VE_estimates_imputed$vaccine_type))){
      this_vax = unique(VE_estimates_imputed$vaccine_type)[t]
      if (! this_vax %in% c("AstraZeneca","Johnson & Johnson","Moderna","Pfizer","Sinopharm","Sinovac" )){
        stop('need to rethink mech shop VE waning implementation line 111')
      }
      
      if (!(d %in% direct$dose & this_vax %in% direct$vaccine_type)){ #avoid if direct exists
        if (this_vax %in% direct$vaccine_type){ # use own distribution if exists
          workshop = apply_distribution %>% filter(vaccine_type == this_vax & strain == s) %>% mutate(dose = d)
          
          ratio = VE_estimates_imputed %>% filter(vaccine_type == this_vax & strain == s & outcome == 'any_infection')
          ratio = min(ratio$VE)/max(ratio$VE)
          
          workshop$VE = workshop$VE * ratio
          
          imputed = rbind(imputed,workshop)
          
        } else { #use AZ for all others (since mode = viral)
          workshop = apply_distribution %>% filter(vaccine_type == "AstraZeneca" & strain == s) %>% 
            mutate(dose = d, vaccine_type = this_vax)
          
          ratio_top = VE_estimates_imputed %>% filter(vaccine_type == this_vax & strain == s & outcome == 'any_infection' & dose == d)
          
          ratio_bottom = VE_estimates_imputed %>% filter(vaccine_type == "AstraZeneca" & strain == s & outcome == 'any_infection' & dose == d)
          
          if (this_vax == "Johnson & Johnson" & d == 1){ 
            #correction otherwise J&J d=1 outtrumps AZ d=1 so entirely that J&J d=1>>d=2
            ratio_bottom = VE_estimates_imputed %>% filter(vaccine_type == "AstraZeneca" & strain == s & outcome == 'any_infection' & dose == 2)
          }
          
          ratio = ratio_top$VE/ratio_bottom$VE
          
          workshop$VE = workshop$VE * ratio
          workshop$VE[workshop$VE>100] = 100
          
          imputed = rbind(imputed,workshop)
          
        } 
      }
    }
  }
}
imputed = imputed %>% select(-VE_internal)

together = rbind(imputed,direct) %>%
  select(strain,vaccine_type,dose,days,VE) %>%
  rename(VE_days = VE) %>%
  mutate(VE_days = VE_days/100)



###(3/3) Plot distributions and save VE_waning_distribution
#(A/B) Plot
waning_to_plot = together %>%
  filter(dose < 3) %>%
  #filter(vaccine_type %in% vax_type_list) %>%
  mutate(immunity = paste(vaccine_type,dose))

strain_test = 'omicron'
ggplot() +
  geom_line(data=waning_to_plot[waning_to_plot$strain == strain_test,],
            aes(x=days,y=VE_days,color=as.factor(immunity)),na.rm=TRUE) +
  labs(title=(paste("Waning of VE against infection","(",strain_test,")"))) +
  xlab("days since vaccination") +
  ylab("% max protection") +
  ylim(0,1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


#(B/B) Save
waning = together %>% mutate(waning = TRUE)
no_waning = together %>% mutate(waning = FALSE) %>%
  group_by(strain,vaccine_type,dose) %>%
  mutate(VE_days = max(VE_days))

VE_waning_distribution = rbind(waning,no_waning) %>% select(strain,vaccine_type,dose,days,VE_days,waning)
save(VE_waning_distribution, file = '1_inputs/VE_waning_distribution.Rdata')