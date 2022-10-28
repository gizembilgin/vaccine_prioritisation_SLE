### This script sets the initial state of the model.
### If the model is running from the fitted date this will be the last next_state.


###### (1/5) Vaccination
#(A/B) Coverage 
#(i/iii) Vaccine coverage at end of known history
vaccine_coverage_end_history = vaccination_history_TRUE %>% 
  filter(date == max(vaccination_history_TRUE$date)) %>%
  select(dose,vaccine_type,age_group,risk_group,coverage_this_date)
#_________________________________________________

#(ii/iii) Add hypothetical campaign (if 'on') ____
if (vax_strategy_toggle == "on" & vax_risk_strategy_toggle == "off"){
  vaccination_history_FINAL = 
    vax_strategy(vax_strategy_start_date        = vax_strategy_toggles$vax_strategy_start_date,
                 vax_strategy_num_doses         = vax_strategy_toggles$vax_strategy_num_doses,
                 vax_strategy_roll_out_speed    = vax_strategy_toggles$vax_strategy_roll_out_speed,
                 vax_age_strategy               = vax_strategy_toggles$vax_age_strategy,  
                 vax_delivery_group             = 'universal',
                 vax_dose_strategy              = vax_strategy_toggles$vax_dose_strategy,            
                 vax_strategy_vaccine_type      = vax_strategy_toggles$vax_strategy_vaccine_type,            
                 vax_strategy_vaccine_interval  = vax_strategy_toggles$vax_strategy_vaccine_interval,            
                 vax_strategy_max_expected_cov  = vax_strategy_toggles$vax_strategy_max_expected_cov
    )
  
  #update attributes!
  list_doses = unique(vaccination_history_FINAL$dose)
  list_doses = list_doses[! list_doses %in% c(8)]
  num_vax_doses = D = length(list_doses)
  num_vax_doses = max(num_vax_doses,vax_strategy_toggles$vax_dose_strategy)
  vax_type_list = sort(unique(vaccination_history_FINAL$vaccine_type))
  num_vax_types = T = length(unique(vaccination_history_FINAL$vaccine_type))
  num_vax_classes = num_vax_doses*num_vax_types + 1                 # + 1 for unvaccinated
  parameters$num_vax_types = num_vax_types
  parameters$num_vax_doses = num_vax_doses
  
} else if (vax_strategy_toggle == "on" & vax_risk_strategy_toggle == "on"){
  
  if (!'risk_group_accessibility' %in% names(apply_risk_strategy_toggles)){
    apply_risk_strategy_toggles$risk_group_accessibility = FALSE
  }
  if (!'risk_group_acceptability' %in% names(apply_risk_strategy_toggles)){
    apply_risk_strategy_toggles$risk_group_acceptability = vax_strategy_toggles$vax_strategy_max_expected_cov
  }
  if (!'risk_group_age_broaden' %in% names(apply_risk_strategy_toggles)){
    apply_risk_strategy_toggles$risk_group_age_broaden = FALSE
  }
  
  vaccination_history_FINAL = 
   apply_risk_strategy(vax_risk_strategy     = apply_risk_strategy_toggles$vax_risk_strategy,            
                       vax_risk_proportion   = apply_risk_strategy_toggles$vax_risk_proportion,      
                       vax_doses_general     = apply_risk_strategy_toggles$vax_doses_general,      
                       vax_doses_risk        = apply_risk_strategy_toggles$vax_doses_risk,
                       risk_group_acceptability = apply_risk_strategy_toggles$risk_group_acceptability,
                       risk_group_accessibility = apply_risk_strategy_toggles$risk_group_accessibility,
                       risk_group_age_broaden   = apply_risk_strategy_toggles$risk_group_age_broaden
   )

  
  #update attributes!
  list_doses = unique(vaccination_history_FINAL$dose)
  list_doses = list_doses[! list_doses %in% c(8)]
  num_vax_doses = D = length(list_doses)
  num_vax_doses = max(num_vax_doses,vax_strategy_toggles$vax_dose_strategy)
  vax_type_list = sort(unique(vaccination_history_FINAL$vaccine_type))
  num_vax_types = T = length(unique(vaccination_history_FINAL$vaccine_type))
  num_vax_classes = num_vax_doses*num_vax_types + 1                 # + 1 for unvaccinated
  parameters$num_vax_types = num_vax_types
  parameters$num_vax_doses = num_vax_doses
  
  date_complete_at_risk_group = vaccination_history_FINAL %>% 
    filter(risk_group == risk_group_name &
             doses_delivered_this_date > 0)
  date_complete_at_risk_group = max(date_complete_at_risk_group$date)
  
  if (is.na(risk_group_lower_cov_ratio) & 
      max(vaccination_history_FINAL$date[vaccination_history_FINAL$dose == 1 & vaccination_history_FINAL$risk_group == risk_group_name & vaccination_history_FINAL$doses_delivered_this_date>0]) >
      max(vaccination_history_FINAL$date[vaccination_history_FINAL$dose == 1 & vaccination_history_FINAL$risk_group != risk_group_name & vaccination_history_FINAL$doses_delivered_this_date>0])){
    stop('high risk group finish dose 1 after lower risk group')
  }
  
  #sum across day in case date fitted < date_now
  if ('FROM_vaccine_type' %in% names(vaccination_history_FINAL)){
    vaccination_history_FINAL = vaccination_history_FINAL %>%
      group_by(date,vaccine_type,vaccine_mode,dose,age_group,risk_group,FROM_dose,FROM_vaccine_type) %>%
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date), .groups = 'keep')
  } else{
    vaccination_history_FINAL = vaccination_history_FINAL %>%
      group_by(date,vaccine_type,vaccine_mode,dose,age_group,risk_group) %>%
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date), .groups = 'keep')
  }

  
} else {
  vaccination_history_FINAL = vaccination_history_TRUE
}


### sensitivity analysis - booster doses in 2023
if (length(booster_toggles)>1){
  source(paste(getwd(),"/(function)_booster_dose_delivery.R",sep=""))

  if (length(booster_prioritised_strategies)>1){
    source(paste(getwd(),"/(function)_prioritised_booster_dose_delivery.R",sep=""))
    
    vaccination_history_FINAL =
      booster_strategy_prioritised(
        booster_risk_strategy = booster_prioritised_strategies$strategy,
        booster_risk_proportion = booster_prioritised_strategies$risk_proportion
      )
    rm(booster_strategy_prioritised)
    
  } else{
    vaccination_history_FINAL =
      booster_strategy( booster_strategy_start_date = booster_toggles$start_date,       # start of hypothetical vaccination program
                        booster_dose_allocation     = booster_toggles$dose_allocation,  # num of doses avaliable
                        booster_rollout_speed       = booster_toggles$rollout_speed,    # doses delivered per day
                        booster_delivery_risk_group = booster_toggles$delivery_risk_group,
                        booster_delivery_includes_previously_boosted = booster_toggles$delivery_includes_previously_boosted,
                        booster_age_strategy        = booster_toggles$age_strategy,     # options: "oldest", "youngest","50_down","uniform"
                        booster_strategy_vaccine_type = booster_toggles$vaccine_type,   # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"
                        booster_strategy_vaccine_interval = booster_toggles$vaccine_interval
      )
    rm(booster_strategy)
  }


  #update attributes!
  list_doses = unique(vaccination_history_FINAL$dose)
  list_doses = list_doses[! list_doses %in% c(8)]
  num_vax_doses = D = length(list_doses)
  vax_type_list = sort(unique(vaccination_history_FINAL$vaccine_type))
  num_vax_types = T = length(vax_type_list)
  num_vax_classes = num_vax_doses*num_vax_types + 1                 # + 1 for unvaccinated
  parameters$num_vax_types = num_vax_types
  parameters$num_vax_doses = num_vax_doses

  #sum across day in case date fitted < date_now
  if ('FROM_vaccine_type' %in% names(vaccination_history_FINAL)){
    vaccination_history_FINAL = vaccination_history_FINAL %>%
      group_by(date,vaccine_type,vaccine_mode,dose,age_group,risk_group,FROM_dose,FROM_vaccine_type) %>%
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date), .groups = 'keep')
  } else{
    vaccination_history_FINAL = vaccination_history_FINAL %>%
      group_by(date,vaccine_type,vaccine_mode,dose,age_group,risk_group) %>%
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date), .groups = 'keep')
  }
}



#UPDATE: Delay & Interval 
vaxCovDelay = crossing(dose = seq(1,num_vax_doses),delay = 0)
vaxCovDelay = vaxCovDelay %>%
  mutate(delay = case_when(
    dose == 1 ~ 21,
    TRUE ~ 14 #all other doses
  ))
#_________________________________________________


#extract other attributes
if (nrow(vaccination_history_FINAL[vaccination_history_FINAL$dose == 8,])>0){
  booster_type = unique(vaccination_history_FINAL$vaccine_type[vaccination_history_FINAL$dose == 8])
  if (booster_type == 'Johnson & Johnson'){
    booster_dose_number = 2
  } else{
    booster_dose_number = 3
  }
}


fitted_seed_dates =  c(as.Date('2021-04-25'),c(as.Date('2021-09-01'))) #COMEBACK - this is hard coded :(

if (fitting == "on"){ #seed date specified in Command Deck
  if (length(fitted_seed_dates[! fitted_seed_dates %in% seed_date])>0 ){stop('seed dates for fitting have changed')}
} else if (outbreak_timing == "after"){ seed_date = max(vaccination_history_FINAL$date)  #outbreak after vaccine rollout
} else if(outbreak_timing == "during"){ seed_date = date_start + 7                       #outbreak during vaccine rollout
} else if(outbreak_timing == "off")   { seed_date = as.Date('1900-01-01')                #no outbreak
}

variant_change_date = c(fitted_seed_dates,seed_date)
variant_change_date = unique(sort(variant_change_date))
variant_change_date = variant_change_date[variant_change_date>as.Date('2021-06-01')] #because delta is not immune escape

date_now = date_start
#_________________________________________________



#(iii/iii)  Initial coverage _______________________
#Including coverage_this_date for projected doses
vaccination_history_FINAL = vaccination_history_FINAL %>% 
  left_join(pop_risk_group_dn, by = c("age_group", "risk_group")) %>%
  group_by(risk_group,age_group,vaccine_type,dose) %>%
  mutate(coverage_this_date = 100*cumsum(doses_delivered_this_date)/pop) %>%
  select(-pop)
vaccination_history_FINAL$coverage_this_date[is.na(vaccination_history_FINAL$coverage_this_date)] = NA

#calculate vaccine coverage of inital state
vaccine_coverage = crossing(risk_group = risk_group_labels,
                            dose = c(1:num_vax_doses),
                            vaccine_type = unique(vaccination_history_FINAL$vaccine_type),
                            age_group = age_group_labels,
                            cov = 0) 
for (r in 1:num_risk_groups){ # risk group
  for (i in 1:J){             # age
    for (t in 1:T){           # vaccine type
      for (d in 1:D){         # vaccine dose
        workshop_type =  unique(vaccination_history_FINAL$vaccine_type)[t]
        workshop_age  = age_group_labels[i]
        workshop_risk = risk_group_labels[r]
        this_vax_max_date = max(vaccination_history_FINAL$date[vaccination_history_FINAL$vaccine_type == workshop_type])
        
          if ((date_start - vaxCovDelay$delay[vaxCovDelay$dose == d])<= this_vax_max_date &
              (date_start - vaxCovDelay$delay[vaxCovDelay$dose == d])>= min(vaccination_history_FINAL$date)){
            
            workshop_value =  vaccination_history_FINAL$coverage_this_date[
              vaccination_history_FINAL$date == date_start - vaxCovDelay$delay[vaxCovDelay$dose == d]
              & vaccination_history_FINAL$age_group == workshop_age
              & vaccination_history_FINAL$risk_group == workshop_risk
              & vaccination_history_FINAL$dose == d
              & vaccination_history_FINAL$vaccine_type == workshop_type] / 100
            
            vaccine_coverage$cov[
                vaccine_coverage$dose == d &
                vaccine_coverage$vaccine_type == workshop_type &
                vaccine_coverage$age_group == workshop_age &
                vaccine_coverage$risk_group == workshop_risk 
            ] = max(workshop_value,0)
             
          } else if ((date_start -vaxCovDelay$delay[vaxCovDelay$dose == d])> this_vax_max_date){
            workshop_value =
              vaccination_history_FINAL$coverage_this_date[
                vaccination_history_FINAL$date == this_vax_max_date
                & vaccination_history_FINAL$dose == d
                & vaccination_history_FINAL$age_group == workshop_age
                & vaccination_history_FINAL$risk_group == workshop_risk
                & vaccination_history_FINAL$vaccine_type == workshop_type]/100 
            
            vaccine_coverage$cov[
              vaccine_coverage$dose == d &
                vaccine_coverage$vaccine_type == workshop_type &
                vaccine_coverage$age_group == workshop_age &
                vaccine_coverage$risk_group == workshop_risk 
              ] =  max(workshop_value,0)
          } 
      }
    }
  }
}
vaccine_coverage$cov[is.na(vaccine_coverage$cov)] = 0
#CHECK
check =  vaccine_coverage %>% group_by(dose,age_group,risk_group) %>% summarise(check = sum(cov),.groups = "keep") %>% filter(check>1)
if(nrow(check)>1){stop('inital vaccine coverage > 100%')}


#(B/B) Vaccine Effectiveness (VE)
load( file = '1_inputs/VE_waning_distribution.Rdata')
VE_waning_distribution = VE_waning_distribution %>%
  filter(waning == waning_toggle_acqusition) %>%
  mutate(outcome = 'any_infection')

#ASSUMPTION - averaging across heterogeneous combinations of VE dose 3
workshop = VE_waning_distribution %>% 
  filter(dose == 3 & strain == strain_now) %>%
  group_by(strain,outcome,vaccine_type,dose,days,waning,.add = TRUE) %>%
  summarise(VE_days = mean(VE_days),.groups = "keep") 
VE_waning_distribution = VE_waning_distribution %>% 
  filter(! dose == 3) %>%
  select(-primary_if_booster)
VE_waning_distribution = rbind(VE_waning_distribution,workshop)

if ((date_start - vaxCovDelay$delay[vaxCovDelay$dose == d])>= min(vaccination_history_POP$date)){
  VE = VE_inital = VE_time_step(strain_inital,date_start,'any_infection')
  #VE_onwards_inital <- VE_time_step(strain_inital,date_start,'transmission')
}

if (fitting == "off"){
  waning_shape_plot_list = list()
  waning_to_plot = VE_waning_distribution %>% 
    filter(dose < 3 & 
             vaccine_type %in% vax_type_list &
             strain == strain_now) %>%
    mutate(immunity = paste(vaccine_type,dose))
  
  for (o in 1:length(unique(waning_to_plot$outcome))){
    waning_shape_plot_list[[o]] = ggplot() +
      geom_line(data=waning_to_plot[waning_to_plot$outcome == unique(waning_to_plot$outcome)[o],],
                aes(x=days,y=VE_days,color=as.factor(immunity)),na.rm=TRUE) +
      labs(title=(paste("Waning of VE against",paste(unique(waning_to_plot$outcome)[o]),"(",strain_now,")"))) +
      xlab("days since vaccination") +
      ylab("% max protection") +
      ylim(0,1)+
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  }
}
#___________________________________________________________________



###### (2/5) Seroprevalence
load(file = "1_inputs/seroprev.Rdata")
seroprev = seroprev %>%
  filter(setting == setting & 
           year == as.numeric(format(date_start, format="%Y")))
if (as.numeric(format(date_start, format="%Y")) > 2022){
  load(file = "1_inputs/seroprev.Rdata")
  seroprev = seroprev[seroprev$setting == setting & seroprev$year ==  2022,]
}
#___________________________________________________________________



###### (3/5) NPI
if(date_start <= max(NPI_estimates$date)){ #if date within known NPI range
  NPI_inital = NPI_estimates$NPI[NPI_estimates$date==date_start]
} else {  #take average of last month on record
  NPI_inital = NPI_estimates %>% 
    filter(date>(max(NPI_estimates$date)-30))
  NPI_inital = mean(NPI_inital$NPI)
}
NPI = NPI_inital = 1 - as.numeric(NPI_inital)/100 #effective NPI is (1-NPI), calculate here so we don't have to repeatedly calculate
#________________________________________________________________



###### (4/5) Hence, initial state
J=num_age_groups
T=num_vax_types
D=num_vax_doses
RISK=num_risk_groups
count=J*(T*D+1)*RISK # +1 is unvax

seed = 0.001*sum(pop) #seed of outbreak at seed_dates (introduction of delta and omicron)

initialInfected = seed*AverageSymptomaticPeriod/(AverageSymptomaticPeriod+AverageLatentPeriod) 
initialExposed  = seed*AverageLatentPeriod/(AverageSymptomaticPeriod+AverageLatentPeriod) 
initialRecovered = seroprev %>% 
  left_join(pop_setting, by = "age_group") %>% 
  mutate(R = seroprev*pop/100) %>%
  select(age_group,R)

date = seq(1,lengthInfectionDerivedImmunity)
date = date_start - date
workshop = as.data.frame(date)
workshop$daily_cases = sum(initialRecovered$R)/lengthInfectionDerivedImmunity
hist_cases = workshop

if (fitting == "on"){

  S_inital=E_inital=I_inital=R_inital=(rep(0,count)) 
  
  #number of active infected/recovered cases
  if (risk_group_toggle == "on"){
    recovered_risk = initialRecovered %>% 
      left_join(risk_dn, by = "age_group") %>%
      mutate(risk_group = risk_group_name,
             R = R*prop) %>% 
      select(risk_group,age_group,R)
    recovered_general_public   = initialRecovered %>% 
      left_join(risk_dn, by = "age_group") %>%
      mutate(risk_group = 'general_public',
             R = R*(1-prop)) %>% 
      select(risk_group,age_group,R)
    initialRecovered = rbind(recovered_general_public,recovered_risk)
  } else{
    initialRecovered = initialRecovered %>%
      mutate(risk_group = 'general_public')
  }
  
  #age distribution of cases
  #LIMITATION: no data so assuming uniform across age groups
  initialClasses = pop_risk_group_dn %>% ungroup() %>% 
    mutate(I = initialInfected*pop/sum(pop),
           E = initialExposed*pop/sum(pop)) %>% 
    left_join(initialRecovered, by = c("risk_group", "age_group")) %>%
    mutate(S = pop - E - I - R) %>%
    select(-pop) %>%
    pivot_longer(
      cols = I:S,
      names_to = 'class',
      values_to = 'state_inital'
    ) 
  if (round(sum(initialClasses$state_inital)) != sum(pop)){stop('size of inital state != size of population')}  
  
  
  # distribute across vaccine classes 
  # LIMITATION: no data so assuming infections are spread equally across vax classes
  empty_unvaccinated = crossing(class = disease_class_list,
                                risk_group = risk_group_labels,
                                dose = 0,
                                vaccine_type = 'unvaccinated',
                                age_group = age_group_labels,
                                state_inital = 0)
  empty_vaccinated = crossing(class = disease_class_list,
                              risk_group = risk_group_labels,
                              dose = seq(1,num_vax_doses),
                              vaccine_type = vax_type_list,
                              age_group = age_group_labels,
                              state_inital = 0)
  state_tidy = rbind(empty_unvaccinated,empty_vaccinated)
  
  for (num in 1:num_disease_classes){
    
    workshop = initialClasses %>% filter(class ==  disease_class_list[num] )
    
    for (r in 1:RISK){
      for (i in 1:J){ # age
        #pop*(1-cov1A-cov1B-cov1C)
        #unvaccinated
        state_tidy$state_inital[state_tidy$dose == 0 & 
                                  state_tidy$age_group == age_group_labels[i] &
                                  state_tidy$risk_group == risk_group_labels[r] &
                                  state_tidy$class == disease_class_list[num]] = 
          workshop$state_inital[workshop$risk_group ==  risk_group_labels[r] & workshop$age_group == age_group_labels[i]]*
          (1-sum(vaccine_coverage$cov[vaccine_coverage$dose == 1 & 
                                        vaccine_coverage$age_group == age_group_labels[i] & 
                                        vaccine_coverage$risk_group == risk_group_labels[r]]))
        
        for (t in 1:T){  # vaccine type
          for (d in 1:D){ # vaccine dose
    
            if (d != D){
              #pop*(cov1A-cov2A)
              state_tidy$state_inital[state_tidy$dose == d  & state_tidy$vaccine_type == vax_type_list[t]& 
                                        state_tidy$age_group == age_group_labels[i] &
                                        state_tidy$risk_group == risk_group_labels[r] &
                                        state_tidy$class == disease_class_list[num]] = 
                workshop$state_inital[workshop$risk_group ==  risk_group_labels[r] & workshop$age_group == age_group_labels[i]]*
                (vaccine_coverage$cov[vaccine_coverage$dose == d &
                                        vaccine_coverage$vaccine_type == vax_type_list[t] &
                                        vaccine_coverage$age_group == age_group_labels[i] & 
                                        vaccine_coverage$risk_group == risk_group_labels[r]] -
                   vaccine_coverage$cov[vaccine_coverage$dose == d+1 &
                                          vaccine_coverage$vaccine_type == vax_type_list[t] &
                                          vaccine_coverage$age_group == age_group_labels[i] & 
                                          vaccine_coverage$risk_group == risk_group_labels[r]])
            }
            if (d == D){
              #pop*cov2A
              state_tidy$state_inital[state_tidy$dose == d & state_tidy$vaccine_type == vax_type_list[t] & 
                                        state_tidy$age_group == age_group_labels[i] & 
                                        state_tidy$risk_group == risk_group_labels[r] &
                                        state_tidy$class == disease_class_list[num]] = 
                workshop$state_inital[workshop$risk_group ==  risk_group_labels[r] & workshop$age_group == age_group_labels[i]]*
                (vaccine_coverage$cov[vaccine_coverage$dose == d &
                                        vaccine_coverage$vaccine_type == vax_type_list[t]&
                                        vaccine_coverage$age_group == age_group_labels[i]& 
                                        vaccine_coverage$risk_group == risk_group_labels[r]])
            }
          }
        }
      }
    }
  }    
  
  state_tidy$class      = factor(state_tidy$class, levels = disease_class_list)
  state_tidy$risk_group = factor(state_tidy$risk_group, levels = risk_group_labels)
  state_tidy$age_group  = factor(state_tidy$age_group, levels = age_group_labels)
  
  state_tidy = state_tidy %>% arrange(class,risk_group,dose,vaccine_type,age_group)
  state_tidy$state_inital[is.na(state_tidy$state_inital)] = 0
  
  
  #construct long array that ODE solver requires
  Incidence_inital=(rep(0,J*(T*D+1)*RISK)) 
  Exposed_incidence_inital = rep(0,J*2)
  state=c(state_tidy$state_inital,Incidence_inital,Exposed_incidence_inital) 
  
  rm (empty_unvaccinated,empty_vaccinated,state_tidy)
  if (round(sum(state)) != sum(pop)){stop('inital state doesnt align with population size!')}
  
} else if (! fitting == "on"){ #load last fitted state of the model
  
  incidence_log = fitted_incidence_log #for rho_inital

  #include additional vaccine types if don't exist
  if (! length(unique(fitted_next_state$vaccine_type)) == length(unique(vaccination_history_FINAL$vaccine_type))+1){ #+1 for unvaccinated
    this_vax = unique(vaccination_history_FINAL$vaccine_type)[!  unique(vaccination_history_FINAL$vaccine_type) %in% unique(fitted_next_state$vaccine_type)]
    copy_vax = unique(vaccination_history_FINAL$vaccine_type)[unique(vaccination_history_FINAL$vaccine_type) %in% unique(fitted_next_state$vaccine_type)][1]
    
    filler = fitted_next_state %>% 
      filter(vaccine_type == copy_vax) %>%
      mutate(pop = 0, 
             vaccine_type = this_vax)
    
    fitted_next_state = rbind(fitted_next_state,filler)
    
    fitted_next_state$class = factor(fitted_next_state$class, levels = disease_class_list)
    fitted_next_state$risk_group = factor(fitted_next_state$risk_group, levels = risk_group_labels)
    fitted_next_state$age_group = factor(fitted_next_state$age_group, levels = age_group_labels)
    fitted_next_state$vaccine_type = factor(fitted_next_state$vaccine_type, levels = vax_type_list)
    
    fitted_next_state = fitted_next_state %>% arrange(class,risk_group,dose,vaccine_type,age_group)
  }
  
  #include additional doses if don't exist
  if (! length(unique(fitted_next_state$dose)) == length(unique(vaccination_history_FINAL$dose[!vaccination_history_FINAL$dose == 8]))+1){  #+1 for unvaccinated
    this_dose = max(fitted_next_state$dose) +1
    copy_dose = max(fitted_next_state$dose)
    
    filler = fitted_next_state %>% 
      filter(dose == copy_dose) %>%
      mutate(pop = 0, 
             dose = this_dose)
    
    fitted_next_state = rbind(fitted_next_state,filler)
    
    fitted_next_state$class        = factor(fitted_next_state$class, levels = disease_class_list)
    fitted_next_state$risk_group   = factor(fitted_next_state$risk_group, levels = risk_group_labels)
    fitted_next_state$age_group    = factor(fitted_next_state$age_group, levels = age_group_labels)
    fitted_next_state$vaccine_type = factor(fitted_next_state$vaccine_type, levels = vax_type_list)
    
    fitted_next_state = fitted_next_state %>% arrange(class,risk_group,dose,vaccine_type,age_group)
  }
  
  S_next=fitted_next_state$pop[fitted_next_state$class == 'S']
  E_next=fitted_next_state$pop[fitted_next_state$class == 'E']
  I_next=fitted_next_state$pop[fitted_next_state$class == 'I']
  R_next=fitted_next_state$pop[fitted_next_state$class == 'R']

  Incidence_inital=(rep(0,J*(T*D+1)*RISK)) 
  Exposed_incidence_inital = rep(0,J*2)
  
  next_state_FINAL=as.numeric(c(S_next,E_next,I_next,R_next,
                                Incidence_inital,Exposed_incidence_inital)) #setting Incid to repeated 0s
  rm(S_next,E_next,I_next,R_next)
  
  if (round(sum(next_state_FINAL)) != sum(pop)){stop('inital state doesnt align with population size!')}
}
#________________________________________________________________



###### (5/5) calculate inital infection-derived protection (rho) & beta
if (waning_toggle_rho_acqusition == TRUE ){
  rho_inital = rho_time_step(date_start)
} else{
  rho_inital = 0.95 #Chemaitelly et al. 2 week estimate
}
if (rho_inital > 1){stop('rho is > 1')}

#LIMITATION WARNING: no age-specific susceptibility to infection is included (no delta data available)
source(paste(getwd(),"/(function)_calculate_R0_Reff.R",sep=""))
beta = rep(beta_fitted_values$beta_optimised[beta_fitted_values$strain == strain_inital],num_age_groups)
#________________________________________________________________

