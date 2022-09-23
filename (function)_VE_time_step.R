### This function calculates the population-level vaccine effectiveness at any given time by vaccine_type, dose and age_group

VE_time_step <- function(strain_now,date_now,outcome){
  
  #(1) load VE_distribution
  VE_distribution <- VE_waning_distribution[VE_waning_distribution$outcome == outcome &
                                              VE_waning_distribution$strain == strain_now,] 
  if (strain_now == 'WT'){
    VE_distribution <- VE_waning_distribution[VE_waning_distribution$outcome == outcome &
                                                VE_waning_distribution$strain == 'delta',] 
  }
  if ('VE_older_adults' %in% names(sensitivity_analysis_toggles)){
    if ('age_group' %in% colnames(VE_distribution)){
      if (length(unique(VE_distribution$age_group)) == 1){
        VE_distribution = VE_distribution %>% ungroup() %>% select(-age_group)
      }
    }
  }
  
  
  #(2) doses delivered to this date
  vax_to_this_date <- vaccination_history_FINAL[vaccination_history_FINAL$date <= date_now,] 
  if (nrow(vaccination_history_FINAL[vaccination_history_FINAL$dose == 8,])>0){
    vax_to_this_date = vax_to_this_date %>% mutate(
    dose = case_when(
      dose == 8 ~ booster_dose_number,
      TRUE ~ dose
      ))
  }
  
  vax_to_this_date <- vax_to_this_date %>% 
    select(risk_group,vaccine_type,dose,date,age_group,doses_delivered_this_date) %>%
    rename(doses = doses_delivered_this_date)
  
  total_doses_up_to_this_date <- vax_to_this_date %>%
    group_by(risk_group,vaccine_type,dose,age_group) %>%
    summarise(total_delivered = sum(doses),.groups = "keep")
    
  vax_to_this_date <- vax_to_this_date %>%
    left_join(total_doses_up_to_this_date, by = c("risk_group", "vaccine_type", "dose", "age_group")) %>%
    mutate(prop = case_when(
      total_delivered >0 ~ doses/total_delivered,
      total_delivered == 0 ~ 0
      ),
      days = as.numeric(date_now - date ))
  
  #<interlude> to add together all days >365 to 365
  meddling <- vax_to_this_date[vax_to_this_date$days > 364,]
  if(length(unique(meddling$days))>1){
    meddling <- meddling %>%
      group_by(risk_group,vaccine_type,dose,age_group) %>%
      summarise(prop = sum(prop),.groups = "keep") %>% 
      mutate(days=365)
    
    vax_to_this_date <- rbind(vax_to_this_date[vax_to_this_date$days<365,c(colnames(meddling))],
                          meddling)
  }

  #(3) Bring VE d'n and AIR history together
  if ('age_group' %in% colnames(VE_distribution) ){
    workshop <- vax_to_this_date %>%
      left_join(VE_distribution, by = c("vaccine_type", "dose", "days",'age_group')) %>%
      select(risk_group,vaccine_type,dose,days,age_group,VE_days,prop) %>%
      mutate(VE_weighted = VE_days*prop)
  } else{
    workshop <- vax_to_this_date %>%
      left_join(VE_distribution, by = c("vaccine_type", "dose", "days")) %>%
      select(risk_group,vaccine_type,dose,days,age_group,VE_days,prop) %>%
      mutate(VE_weighted = VE_days*prop)
  }
  
  
  #(4) Aggregate to estimate population VE for doses
  workshop <- workshop %>%
    group_by(risk_group,dose,vaccine_type,age_group) %>%
    summarise(VE = sum(VE_weighted),.groups = "keep")
  
  if(nrow(workshop[round(workshop$VE,digits=2)>1,])){stop('VE > 1!')}
  
  #<interim> add none covered vaccines
  for (i in 1:J){
    for (t in 1:num_vax_types){
      for (d in 1:D){
        for (r in 1:num_risk_groups){
          this_vax = vax_type_list[t]
          if (!( this_vax %in% unique(workshop$vaccine_type[workshop$risk_group == risk_group_labels[r] & workshop$dose == d & workshop$age_group == age_group_labels[i]]))){
            workshop2 = crossing(risk_group = risk_group_labels[r],
                                 dose = d,
                                 vaccine_type = this_vax,
                                 age_group = age_group_labels[i],
                                 VE =0) 
            workshop = rbind(workshop,workshop2)
          } 
        }
      }
    }
  }

  workshop[is.na(workshop)] <-0 
  
  VE_tidy = workshop
  
  return(VE_tidy)
  
}