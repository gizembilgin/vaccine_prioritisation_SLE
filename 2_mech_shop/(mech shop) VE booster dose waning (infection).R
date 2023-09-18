### This (mech shop) applies a waning distribution to the point estimates of booster VE against infection
### Creates: VE_waning_distribution (updated from `(mech shop) VE waning distribution (infection)` to include booster doses)

require(ggpubr); require(readr); require(gridExtra); require(ggplot2); require(tidyverse);


### Showing that Pfizer booster to CoronaVac EQ to Pfizer dose 2 waning
#rm(list=ls())
raw <- read.csv(file = '1_inputs/VE_acq.csv',header=TRUE)
raw = raw %>% filter(dose == 3 & age_group == 'overall') %>% select(dose,days,VE) %>% mutate(label = 'Pfizer booster with CoronaVac primary schedule')

raw_2 <- read.csv(file = '1_inputs/VE_Andrews_shape.csv',header=TRUE)
raw_2 = raw_2 %>% filter(dose == 2 & vaccine_type == 'Pfizer' & strain == 'omicron') %>% mutate(days = week * 7, VE = VE/100,label = 'Pfizer second dose') %>% select(dose,days,VE,label)

raw = rbind(raw,raw_2)

predicted_distribution = data.frame()

  subplot_list = list()
  for (i in 1: length(unique(raw$dose))){
    workshop_real = raw[raw$dose == unique(raw$dose)[i],]
    this_label = unique(workshop_real$label)
    attach(workshop_real)
    model = lm(log(VE)~days)
    
    #summary(model)
    model_rsquared = summary(model)$adj.r.squared
    detach(workshop_real)
    
    time <- seq(0, 365)
    workshop_predicted <- exp(predict(model,list(days=time)))
    workshop_predicted = data.frame(cbind(days = time, VE = workshop_predicted))
    
    subplot_list[[i]] = ggplot() + 
      geom_point(data=workshop_real,aes(x=days,y=VE)) +
      geom_line(data = workshop_predicted, aes(x=days,y=VE)) +
      #ylim(0,1)+
      #ggtitle(paste(unique(raw$dose)[i],round(model_rsquared,2)))
      ggtitle(paste(this_label)) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(color = 'black'))
    
    workshop_predicted = workshop_predicted %>%
      mutate(dose = unique(raw$dose)[i],
             rsquared = model_rsquared )
    
    predicted_distribution = rbind(predicted_distribution,workshop_predicted)
  }
  grid.arrange(subplot_list[[1]], subplot_list[[2]], nrow=2)
  #_________________________________________________________________________________________________________________________________________

  
  ### Internal VE
  apply_distribution <- predicted_distribution %>%
    group_by(dose) %>%
    mutate(VE_internal = VE / max(VE)) %>%
    ungroup() %>%
    select(dose,days,VE_internal,VE)
  ggplot() + 
    geom_line(data = apply_distribution, aes(x=days,y=VE_internal,color=as.factor(dose))) +
    ylim(0,1)
  apply_distribution = apply_distribution %>% filter(dose == 2) %>% select(-dose) %>% mutate(outcome_family = 'acquisition')
  #_________________________________________________________________________________________________________________________________________
  
  
  ### Apply distribution
  load(file = "1_inputs/VE_booster_estimates.Rdata")
  
  imputed = data.frame()
    for (t in 1:length(unique(VE_booster_estimates$primary_if_booster))){
         this_vax = unique(VE_booster_estimates$primary_if_booster)[t]
      
          workshop = apply_distribution  %>% 
            mutate(primary_if_booster = this_vax)
          
          ratio_top = VE_booster_estimates %>% filter(primary_if_booster == this_vax & outcome == 'any_infection')
          ratio_bottom = raw %>% filter(label == 'Pfizer second dose' & days == min(days))
          
          ratio = (ratio_top$VE/100)/ratio_bottom$VE
          
          workshop$VE = workshop$VE * ratio
          workshop$VE[workshop$VE>100] = 100
          
          imputed = rbind(imputed,workshop)
    }

  imputed = imputed %>% 
    select(-VE_internal) %>%
    rename(VE_days=VE)
  
  plot_1 = ggplot() + 
    geom_line(data = imputed, aes(x=days,y=VE_days,color=as.factor(primary_if_booster))) +
    theme_bw() + 
    xlab("") + 
    ylim(0,1) +
    labs(title=(paste("VE against acqusition")),color='primary schedule') +
    ylab("") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(color = 'black'))
 
  ###load booster from severe outcome `(mech shop) VE waning distribution (severe outcomes).R`
  check = waning_to_plot %>% filter(is.na(primary_if_booster) == FALSE)
  plot_2 = ggplot() +
    geom_line(data=check[check$strain == 'omicron' & check$outcome  == 'death',],
              aes(x=days,y=VE_days,color=as.factor(primary_if_booster)),na.rm=TRUE) +
    labs(title=(paste("VE against death")),color='primary schedule') +
    xlab("") +
    ylim(0,1) +
    ylab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  plot_3 = ggplot() +
    geom_line(data=check[check$strain == 'omicron' & check$outcome  == 'severe_disease',],
              aes(x=days,y=VE_days,color=as.factor(primary_if_booster)),na.rm=TRUE) +
    labs(title=(paste("VE against severe_disease")),color='primary schedule') +
    xlab("days since vaccination") +
    ylim(0,1) +
    ylab("") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  grid.arrange(plot_1,plot_2,plot_3, nrow =3 )
  #_________________________________________________________________________________________________________________________________________
  

  
  
  ### Apply distribution
  waning = imputed %>% mutate(waning = TRUE)
  no_waning = imputed %>% mutate(waning = FALSE) %>%
    group_by(primary_if_booster) %>%
    mutate(VE_days = max(VE_days))
  
  together = rbind(waning,no_waning) %>% 
    mutate(strain = 'omicron', dose = 3, vaccine_type = "Pfizer") %>%
    select(strain,vaccine_type,dose,days,VE_days,waning,primary_if_booster)
  
  load(file = '1_inputs/VE_waning_distribution.Rdata')
  VE_waning_distribution = bind_rows(VE_waning_distribution,together)
  
  save(VE_waning_distribution, file = '1_inputs/VE_waning_distribution.Rdata')
  
  
  