### This script calculates the incidence of severe outcomes by age_group and vaccination status per time step of the model
### Incidences of severe outcomes vary across time due to the waning of vaccine- and infection-derived immunity

# VE_loop set to 0 when no sensitivity analysis of reduced VE in older adults or adults with comorbidities is conducted
if (exists("VE_loop") == FALSE){ VE_loop = 0} 
if (VE_loop == 0){
  if (length(unique(VE_waning_distribution$outcome)) == 1){ #'any_infection'
    save_VE_waning_distribution = VE_waning_distribution
  } else{
    VE_waning_distribution = save_VE_waning_distribution
  }

  load( file = '1_inputs/VE_waning_distribution_SO.Rdata')
  VE_waning_distribution_SO = VE_waning_distribution_SO %>% filter(waning == waning_toggle_severe_outcome )
  VE_waning_distribution = bind_rows(VE_waning_distribution,VE_waning_distribution_SO)
  
  workshop = VE_waning_distribution %>% 
    filter(dose == 3 & strain == strain_now) %>%
    group_by(strain,outcome,vaccine_type,dose,days,waning,.add = TRUE) %>%
    summarise(VE_days = mean(VE_days),.groups = "keep")
  #ggplot(workshop) +  geom_line(data=workshop[workshop$waning == TRUE,],aes(x=days,y=VE_days,linetype = as.factor(outcome) ))
  VE_waning_distribution = VE_waning_distribution %>% filter(! dose == 3) %>% select(-primary_if_booster)
  VE_waning_distribution = rbind(VE_waning_distribution,workshop)
  
} else if (VE_loop == 1 &'VE_older_adults' %in% names(sensitivity_analysis_toggles)){
  #Note: VE_loop == 2 (comorb) will use this same dn
  load( file = '1_inputs/SA_VE_older_muted_SO.Rdata')
  VE_waning_distribution_SO = SA_VE_older_muted_SO %>% filter(waning == waning_toggle_severe_outcome )
  VE_waning_distribution = bind_rows(save_VE_waning_distribution,VE_waning_distribution_SO)
  VE_waning_distribution = VE_waning_distribution %>% group_by(age_group)
  
  workshop = VE_waning_distribution %>% 
    filter(dose == 3 & strain == strain_now) %>%
    group_by(strain,outcome,vaccine_type,dose,days,waning,.add = TRUE) %>%
    summarise(VE_days = mean(VE_days),.groups = "keep")
  #ggplot(workshop) +  geom_line(data=workshop[workshop$waning == TRUE,],aes(x=days,y=VE_days,linetype = as.factor(outcome) ))
  VE_waning_distribution = VE_waning_distribution %>% filter(! dose == 3) %>% select(-primary_if_booster)
  VE_waning_distribution = rbind(VE_waning_distribution,workshop)
}  



#(A/C) calculate VE against severe outcomes by day
if (waning_toggle_severe_outcome == TRUE){
  if (ticket == 1 | vax_strategy_toggle == "on"){ #only have to run when vax strategy changing
    VE_tracker = data.frame()
    for (outcome in c('death','severe_disease')){
      for (day in 0:num_time_steps){
        workshop = VE_time_step(strain_now,date_start+day,outcome)
        workshop = workshop %>% mutate(date=day,
                                       outcome_VE=outcome)
        VE_tracker = rbind(VE_tracker,workshop)
      }
    }
    VE_tracker$date = date_start + VE_tracker$date
  }
} else{
  VE_tracker= crossing(risk_group = risk_group_labels,
                       dose = seq(1,D),
                       vaccine_type = vax_type_list,
                       age_group = age_group_labels,
                       date = unique(incidence_log$date),
                       outcome_VE = c('death','severe_disease'))
  workshop = VE_estimates_imputed %>% 
    filter(strain == strain_now) %>% 
    select(vaccine_type,dose,outcome,VE) %>% 
    mutate(VE=VE/100) %>%
    rename(outcome_VE = outcome)
  VE_tracker = VE_tracker %>% left_join(workshop, by = c("dose", "vaccine_type", "outcome_VE"))
}

workshop = crossing(risk_group = risk_group_labels,
                    dose = 0,
                    vaccine_type = "unvaccinated",
                    age_group = age_group_labels,
                    date = unique(incidence_log$date),
                    outcome_VE = c('death','severe_disease'),
                    VE = 0)
VE_tracker = rbind(VE_tracker,workshop)

if ('VE_adults_comorb' %in% names(sensitivity_analysis_toggles)){
  VE_tracker =  VE_tracker %>%
    mutate(VE = case_when(
      risk_group == risk_group_name & age_group %in% c("30 to 44","45 to 59") ~ VE * sensitivity_analysis_toggles$VE_adults_comorb,
      TRUE ~ VE
    ))
}



#(B/C) 
load(file = '1_inputs/severe_outcome_FINAL.Rdata')
if (risk_group_toggle == "on"){
  #if risk-group included, then adjust general population incidence rate so pop-level estimates stay the same

  severe_outcomes_list = unique(severe_outcome_FINAL$outcome)
  severe_outcome_FINAL_wRisk = data.frame()
  
  for (o in 1:length(severe_outcomes_list)){
    this_outcome = severe_outcomes_list[o]
    for (i in 1:num_age_groups){
      this_age = age_group_labels[i]
      row = severe_outcome_FINAL %>% filter(outcome == this_outcome & age_group == this_age)
      
      IR = row$percentage
      P = pop_setting$pop[pop_setting$age_group == this_age]
      P_general = pop_risk_group_dn$pop[pop_risk_group_dn$risk_group == 'general_public' & pop_risk_group_dn$age_group == this_age]
      P_risk = pop_risk_group_dn$pop[pop_risk_group_dn$risk_group == risk_group_name & pop_risk_group_dn$age_group == this_age]
      if (P != (P_general+P_risk)){stop('Line 15 P_general + P_risk != P_overall')}
      
      IR_gen = ((IR*P)/((RR_estimate*P_risk)/P_general +1))/P_general
      IR_risk = ((IR*P)/(P_general/(RR_estimate*P_risk) +1))/P_risk
      
      row_gen = row %>% mutate(percentage = IR_gen,risk_group = 'general_public')
      row_risk = row %>% mutate(percentage = IR_risk,risk_group = risk_group_name)
      severe_outcome_FINAL_wRisk = rbind(severe_outcome_FINAL_wRisk,row_gen,row_risk)
    }
  }
  
  severe_outcome_FINAL_wRisk$percentage[is.na(severe_outcome_FINAL_wRisk$percentage)]=0
  severe_outcome_FINAL = severe_outcome_FINAL_wRisk
}



#(C/C) calculate severe outcome incidence by vax_status
if (risk_group_toggle == "on"){

  if (risk_group_name == "pregnant_women"){  # Add excess neonatal deaths due to increased risk of stillbirth or preterm delivery
    #ASSUMPTION - no dependence on pregnant woman's age
    stillbirth_prev = 15.6/1000
    stillbirth_risk = 1.80
    preterm_prev = 6.2/1000 
    preterm_risk = 1.47
    
    row = crossing(outcome = 'neonatal_deaths',
                              outcome_long = 'adverse pregnancy outcome resulting in neonatal death',
                              age_group  = age_group_labels,
                              percentage = ((stillbirth_risk-1)*stillbirth_prev + (preterm_risk-1)*preterm_prev),
                              outcome_VE = 'severe_disease',
                              risk_group = 'pregnant_women') 
    severe_outcome_FINAL = rbind(severe_outcome_FINAL,row)
  }
  
  severe_outcome_this_run = severe_outcome_FINAL %>% 
    left_join(VE_tracker, by = c("age_group", "outcome_VE", "risk_group")) %>%
    mutate(percentage = percentage*(1-VE)) %>%
    select(date,outcome,outcome_long,age_group,risk_group,vaccine_type,dose,percentage)
  
} else{
  severe_outcome_this_run = severe_outcome_FINAL %>% 
    left_join(VE_tracker, by = c("age_group", "outcome_VE")) %>%
    mutate(percentage = percentage*(1-VE)) %>%
    select(date,outcome,outcome_long,age_group,risk_group,vaccine_type,dose,percentage)
}
#_______________________________________________________________________________
