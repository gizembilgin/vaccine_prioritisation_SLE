### This function allocates booster doses by age group, risk group and dose
### It has been designed for:
### (1) sensitivity analysis S4.3 in the vaccine allocation paper - allocation of additional doses
### (2) antiviral model - allocation of booster doses in sequential years regardless of how many doses they have previously had (all doses after the 3rd assumed roughly equal to VE(3rd dose))


###### Coding vaccine prioritisation strategies
#options(scipen = 100) #removes scientific notation
#options(warn = 0) # = 2 when debugging
#NOTE: limitation that hypoth vaccine_type must be consistent with prev hypoth vaccine type!

booster_strategy <- function(
                         booster_strategy_start_date,       # start of hypothetical vaccination program
                         booster_dose_allocation,        # num of doses avaliable
                         booster_rollout_speed,   # doses delivered per day
                         booster_delivery_risk_group = c(risk_group_name,'general_public'),
                         booster_delivery_includes_previously_boosted = 'Y',
                         booster_age_strategy,              # options: "oldest", "youngest","50_down","uniform"
                         booster_strategy_vaccine_type,     # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"  
                         booster_strategy_vaccine_interval, # (days) since primary schedule
                         booster_prioritised = 'N',
                         vaccination_history_FINAL_local = vaccination_history_FINAL
                
){
  
  ##### SETUP ###################################################################
  ### WARNINGS 
  if (booster_strategy_start_date <= max(vaccination_history_TRUE$date)){ 
    warning ('Your hypothetical vaccine campaign start date needs to be in the future!')
  }
  if (!(booster_strategy_vaccine_type %in% c("Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"))){
    stop('pick a valid vaccine type, or check your spelling!')
  }


  
  ### IMPORTS
  booster_strategy_vaccine_interval = round(booster_strategy_vaccine_interval)
  if (booster_strategy_vaccine_interval < 60){warning('Are you sure that there is less than 2 months between primary and booster doses?')}
  prioritisation_csv <- read.csv("1_inputs/prioritisation_strategies.csv",header=TRUE)
  
  if (booster_strategy_vaccine_type == "Johnson & Johnson" ){booster_dose_number = 2
  } else{  booster_dose_number = 3}
  #_______________________________________________________________________________
  
  
 
  ##### BOOSTER DOSES ###################################################################################################################################
  ### create booster_speed_modifier
  
  booster_speed_modifier = data.frame()
  if (booster_strategy_start_date <= max(vaccination_history_FINAL_local$date)){
    booster_speed_modifier = vaccination_history_FINAL_local %>% 
      filter(date >= max(vaccination_history_TRUE$date) &
               date >= booster_strategy_start_date) %>%
      group_by(date) %>% 
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep")  %>% 
      ungroup() %>% 
      mutate(doses_avaliable = booster_rollout_speed - doses_delivered_this_date) %>%
      select(date,doses_avaliable)
  }
  #additional dates to make up to two years
  additional_dates = seq(booster_strategy_start_date, booster_strategy_start_date + 52*7*2, by = 'days')
  additional_dates = additional_dates[! additional_dates %in% booster_speed_modifier$date]
  if (length(additional_dates) > 0 ){
    additional_dates = data.frame(date = additional_dates)
    additional_dates = additional_dates %>% mutate(doses_avaliable = booster_rollout_speed)
    booster_speed_modifier = rbind(booster_speed_modifier,additional_dates)
  }
  booster_speed_modifier = booster_speed_modifier %>%
    mutate(day = as.numeric(date - booster_strategy_start_date + 1),
           cumsum = cumsum(doses_avaliable))
  
  
  ###<intermission check_elig>###
  check_elig = vaccination_history_FINAL_local  %>%
    filter(risk_group %in% booster_delivery_risk_group)  %>%
    mutate(date = date + booster_strategy_vaccine_interval) 
  if (booster_delivery_includes_previously_boosted == "N"){
    check_elig = check_elig %>%
      mutate(primary_schedule_complete = case_when( #ASSUMPTION - this flag isn't used, so assuming boosted whether or not completed primary schedule
        vaccine_type == "Johnson & Johnson" ~  'Y',
        dose %in% c(2,3) ~ 'Y',
        TRUE ~ 'N'), 
        boosted = case_when(
          vaccine_type == "Johnson & Johnson" & dose == 2 ~ 'Y',
          dose == 3 ~ 'Y',
          TRUE ~ 'N')) %>% 
      filter(boosted == booster_delivery_includes_previously_boosted)
  }
  check_elig = check_elig %>% 
    ungroup() %>% group_by(date,dose) %>% summarise(total_doses = sum(doses_delivered_this_date),.groups = "keep") %>%
    ungroup() %>% mutate(eligible_individuals = cumsum(total_doses)) %>%
    select(-total_doses)
  
  #remove double counted tidy
  for (d in 2:num_vax_doses){
    remove = check_elig %>% 
      ungroup() %>%
      filter(dose == d) %>%
      rename(complete_vax = eligible_individuals) %>%
      select(date,complete_vax)
    
    if (nrow(remove)>0){
      check_elig = check_elig %>% 
        left_join(remove, by = c('date')) %>%
        mutate(eligible_individuals = case_when(
          dose == (d-1) & complete_vax > eligible_individuals ~ 0, #this shouldn't be triggered
          dose == (d-1) ~ eligible_individuals - complete_vax,
          TRUE ~ eligible_individuals,
        )) %>%
        select(-complete_vax)
    }
  }
  if (nrow(check_elig[check_elig$dose == 8,])>0){
    remove = check_elig  %>%
      ungroup() %>%
      filter(dose == 8) %>%
      select(-vaccine_type,-dose) %>%
      rename(boosted_vax = eligible_individuals,
             dose = FROM_dose, 
             vaccine_type = FROM_vaccine_type)
    
    if (nrow(remove)>0){
      check_elig = check_elig %>% 
        left_join(remove, by = c('age_group','vaccine_type','risk_group','dose')) %>%
        mutate(eligible_individuals = case_when(
          is.na(boosted_vax) ~ eligible_individuals, #this shouldn't be triggered
          TRUE ~ eligible_individuals - boosted_vax,
        )) %>%
        select(-boosted_vax)
    }
  }
  check_elig = check_elig %>% 
    ungroup() %>% group_by(date) %>% summarise(eligible_individuals = sum(eligible_individuals))

  booster_speed_modifier = booster_speed_modifier %>%
    left_join(check_elig, by = "date")
  
  #make sure that individuals are eligible when boosted
  check_less_than = booster_speed_modifier  %>%
    filter(cumsum > eligible_individuals) 
    
  while (nrow(check_less_than) > 0) {
    booster_speed_modifier = booster_speed_modifier %>%
      mutate(difference = cumsum - eligible_individuals)
    this_date = check_less_than$date[1]
    booster_speed_modifier$doses_avaliable[booster_speed_modifier$date == this_date] =
      booster_speed_modifier$doses_avaliable[booster_speed_modifier$date == this_date] - booster_speed_modifier$difference[booster_speed_modifier$date == this_date]
    booster_speed_modifier = booster_speed_modifier %>%
      mutate(cumsum = cumsum(doses_avaliable))
    check_less_than = booster_speed_modifier %>% filter(cumsum > eligible_individuals)
  }
  ###
  
  
  #####(1/4) Calculate the eligible population ###################################
  if (8 %in% unique(vaccination_history_FINAL_local$dose)){
    eligible_pop =  vaccination_history_FINAL_local %>% 
      group_by(dose,vaccine_type,risk_group,age_group,FROM_vaccine_type,FROM_dose) %>%
      summarise(eligible_individuals = sum(doses_delivered_this_date) , .groups = 'keep')
  } else{
    eligible_pop =  vaccination_history_FINAL_local %>% 
      group_by(dose,vaccine_type,risk_group,age_group) %>%
      summarise(eligible_individuals = sum(doses_delivered_this_date), .groups = 'keep')
  }

 
  
  #remove double counted tidy
  for (d in 2:num_vax_doses){
    remove = eligible_pop %>% 
      ungroup() %>%
      filter(dose == d) %>%
      rename(complete_vax = eligible_individuals) %>%
      select(vaccine_type,risk_group,age_group,complete_vax)
    
    if (nrow(remove)>0){
      eligible_pop = eligible_pop %>% 
        left_join(remove, by = c('age_group','vaccine_type','risk_group')) %>%
        mutate(eligible_individuals = case_when(
          dose == (d-1) & complete_vax > eligible_individuals ~ 0, #this shouldn't be triggered
          dose == (d-1) ~ eligible_individuals - complete_vax,
          TRUE ~ eligible_individuals,
        )) %>%
        select(-complete_vax)
    }
  }
  if (nrow(eligible_pop[eligible_pop$dose == 8,])>0){
    remove = eligible_pop  %>%
      ungroup() %>%
      filter(dose == 8) %>%
      select(-vaccine_type,-dose) %>%
      rename(boosted_vax = eligible_individuals,
             dose = FROM_dose, 
             vaccine_type = FROM_vaccine_type)
    
    if (nrow(remove)>0){
      eligible_pop = eligible_pop %>% 
        left_join(remove, by = c('age_group','vaccine_type','risk_group','dose')) %>%
        mutate(eligible_individuals = case_when(
          is.na(boosted_vax) ~ eligible_individuals, #this shouldn't be triggered
          TRUE ~ eligible_individuals - boosted_vax,
        )) %>%
        select(-boosted_vax)
    }
  }
  
  eligible_pop = eligible_pop %>%
    mutate(primary_schedule_complete = case_when(
      vaccine_type == "Johnson & Johnson" ~  'Y',
      dose %in% c(2,3,8) ~ 'Y',
      TRUE ~ 'N'),
      boosted = case_when(
        vaccine_type == "Johnson & Johnson" & dose == 2 ~ 'Y',
        dose %in% c(3,8) ~ 'Y',
        TRUE ~ 'N')) %>%
    filter(
      eligible_individuals > 0,
      risk_group %in% booster_delivery_risk_group &
        # primary_schedule_complete == "Y" & #CAN CHANGE THIS TOGGLE #ASSUMPTION / COMEBACK - including those with only 1 primary dose
        boosted == booster_delivery_includes_previously_boosted
    )
  
  #collapse dose == 8
  eligible_pop = eligible_pop %>%
    group_by(dose,vaccine_type,risk_group,age_group) %>%
    summarise(eligible_individuals = sum(eligible_individuals), .groups = 'keep')
  
  if (nrow(eligible_pop) == 0){
    warning('no one is eligible for this additional booster')
    return()
  }
  #_______________________________________________________________________________
  
  
  
  #####(2/4) Place priority # on age group by strategy ############################
  if (booster_age_strategy %in% unique(prioritisation_csv$strategy)) {
    elected_strategy = prioritisation_csv[prioritisation_csv$strategy == booster_age_strategy,c('age_group','priority')]
    eligible_pop <- eligible_pop %>% 
      left_join(elected_strategy, by = "age_group")
  } else if (booster_age_strategy == "manual_overwrite"){
    eligible_pop <- eligible_pop %>%
      mutate(priority= case_when(
        age_group == '70 to 100'~ 7,
        age_group == '60 to 69' ~ 6,
        age_group == '45 to 59' ~ 5,
        age_group == '30 to 44' ~ 4,
        age_group == '18 to 29' ~ 3,
        age_group == '10 to 17' ~ 2,
        age_group == '5 to 9' ~ 1,
        age_group == '0 to 4' ~ 99
      ))
  }
  #Note: 99 used as a 'never' indicator
  #_______________________________________________________________________________
  
  
  
  #####(3/4)  Distribute doses by priority #######################################
  eligible_pop <- eligible_pop %>% mutate(doses_delivered = 0)
  doses_to_deliver = booster_dose_allocation
  priority_num = 1
  n = length(unique(eligible_pop$priority))
  highest_priority = 1
  if (length(unique(eligible_pop$priority))>1){  highest_priority = sort(unique(eligible_pop$priority), partial = n - 1)[n - 1] }#highest valid priority

  while (doses_to_deliver > 0 &
         priority_num <= (highest_priority)) {
    priority_group = eligible_pop[eligible_pop$priority == priority_num, ]
    
    workshop = doses_to_deliver - sum(priority_group$eligible_individuals) #check enough for all priority groups
    
    if (workshop < 0) {
      prop_to_deliver = doses_to_deliver / sum(priority_group$eligible_individuals)
      priority_group$doses_delivered = prop_to_deliver * priority_group$eligible_individuals
      doses_to_deliver = 0
    } else if (workshop > 0) {
      priority_group$doses_delivered =  priority_group$eligible_individuals
      doses_to_deliver = doses_to_deliver - sum(priority_group$eligible_individuals)
    }
    
    eligible_pop$doses_delivered[eligible_pop$priority == priority_num] = priority_group$doses_delivered
    
    priority_num = priority_num + 1
  }
  if (nrow(eligible_pop[eligible_pop$doses_delivered > eligible_pop$eligible_individuals, ]) >
      0) {
    stop('more doses delivered than eligible individuals!')
  }
    #_______________________________________________________________________________
    
    
    
    #####(4/4) Distribute between days #############################################
    # we should use:
    # (1) booster_dose_allocation - doses to deliver
    # (2) booster_rollout_speed - max doses delivered per day
    # (3) booster_strategy_start_date - first day of doses delivered
    
    VA =  eligible_pop %>% mutate(doses_left = doses_delivered)
    priority_num = 1
    priority_group  = as.character(unique(VA$age_group[VA$priority == priority_num]))
    
    ceiling = min(sum(eligible_pop$doses_delivered),booster_dose_allocation)
    timeframe = ceiling/booster_rollout_speed
    daily_per_dose = booster_rollout_speed
    
    #modify if available doses restricted
    if (nrow(booster_speed_modifier[booster_speed_modifier$doses_avaliable != booster_rollout_speed,])>0){
      booster_speed_modifier$day = booster_speed_modifier$date - booster_strategy_start_date + 1
      
      exception_list = unique(booster_speed_modifier$day)
      workshop = booster_speed_modifier %>% filter(cumsum<ceiling)
      workshop = nrow(workshop) + 1
      
      if (max(booster_speed_modifier$cumsum)+daily_per_dose>ceiling){
        timeframe = workshop
      } else{
        extension_time = (ceiling - max(booster_speed_modifier$cumsum))/daily_per_dose
        timeframe = timeframe + extension_time
      }
    } else{
      exception_list = c(-1)
    }
    
    timeframe = ceiling(timeframe)
    
    booster_delivery_outline <- crossing(day = c(1:timeframe),
                                         dose = unique(eligible_pop$dose), 
                                         vaccine_type = unique(eligible_pop$vaccine_type),
                                         age_group = age_group_labels,
                                         risk_group = booster_delivery_risk_group,
                                         doses_delivered = c(0))
    
    for (day in 1:timeframe){
      avaliable = daily_per_dose
      daily_per_dose_here = daily_per_dose
      
      #apply limiting
      if (day %in% exception_list){
        avaliable = booster_speed_modifier$doses_avaliable[booster_speed_modifier$day == day]
        daily_per_dose_here = avaliable 
      }
      
      #ensuring that we don't overshoot available doses
      if (day == timeframe){
        total_delivered = sum(booster_delivery_outline$doses_delivered)
        avaliable = min(booster_dose_allocation-total_delivered,daily_per_dose,ceiling-total_delivered)
        
        #CHECK
        total_delivered = total_delivered + avaliable
        if(! total_delivered == ceiling){stop('delivered booster doses do not align with expected')}
        
      }
      if (avaliable > sum(VA$doses_left)){avaliable = sum(VA$doses_left)}
      
      while(avaliable>0 & priority_num <= highest_priority){
        
        if(sum(VA$doses_left[VA$priority == priority_num])>0){ 
          #i.e., while we still have doses to deliver in this priority group
          
          if(0 %in% VA$doses_left[VA$priority == priority_num & VA$dose == 1]){
            age_complete = VA$age_group[VA$doses_left == 0 & VA$priority == priority_num & VA$dose == 1]
            VA$priority[VA$age_group %in% age_complete] = VA$priority[VA$age_group %in% age_complete] + 100
            
            priority_group = as.character(unique(VA$age_group[VA$priority == priority_num]))
          } #FIX - when one age group in the priority group runs out first
          
          workshop_doses = min(sum(VA$doses_left[VA$priority == priority_num]),
                               daily_per_dose_here,
                               avaliable)
          #either deliver max capacity or number left in this group, whichever is fewer
          
          leftover=0
          VA_pt = VA #snapshot
          
          for (i in length(priority_group):1){
            for (d in unique(VA$dose)){
              for (t in 1:length(unique(VA$vaccine_type))){
                for (r in 1:length(unique(VA$risk_group))){
                  workshop_age = priority_group[i]
                  workshop_type = unique(VA$vaccine_type)[t]
                  workshop_risk_group = unique(VA$risk_group)[r]
                  
                  if (nrow(VA_pt[VA_pt$age_group == workshop_age & VA_pt$dose == d & VA_pt$vaccine_type == workshop_type & VA_pt$risk_group == workshop_risk_group,])>0){
                    workshop_prop = VA_pt$doses_left[VA_pt$age_group == workshop_age & VA_pt$dose == d & VA_pt$vaccine_type == workshop_type & VA_pt$risk_group == workshop_risk_group]/
                      sum(VA_pt$doses_left[VA_pt$priority == priority_num])
                    workshop_calc = workshop_doses * workshop_prop + leftover
                    
                    if (workshop_calc > VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type & VA$risk_group == workshop_risk_group]){
                      leftover = abs(workshop_calc - VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type & VA$risk_group == workshop_risk_group])
                      workshop_calc = VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type & VA$risk_group == workshop_risk_group]
                      
                    } else{
                      leftover = 0
                    }
                    booster_delivery_outline$doses_delivered[booster_delivery_outline$day == day &
                                                               booster_delivery_outline$dose == d &
                                                               booster_delivery_outline$vaccine_type == workshop_type &
                                                               booster_delivery_outline$age_group == workshop_age &
                                                               booster_delivery_outline$risk_group == workshop_risk_group] = workshop_calc
                    VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type & VA$risk_group == workshop_risk_group] = 
                      VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type & VA$risk_group == workshop_risk_group] - workshop_calc
                  }
                }
              }
            }
          }
          
          avaliable = avaliable - workshop_doses
          
          
        } else if (sum(VA$doses_left[VA$priority == priority_num])==0){
          priority_num = priority_num+1
          priority_group = as.character(unique(VA$age_group[VA$priority == priority_num]))
        } else{
          stop('negative doses left, reconsider!')
        }
      } #<end while loop>
      
    }
    
    ### formating booster_delivery_outline to align with vaccination_history_TRUE
    booster_delivery_outline = booster_delivery_outline %>%
      mutate(date = booster_strategy_start_date + (day-1)) %>%
      rename(doses_delivered_this_date = doses_delivered,
             FROM_vaccine_type = vaccine_type) %>%
      mutate(FROM_dose = case_when(dose != 8 ~ dose, 
                                   dose == 8 ~ booster_dose_number), 
             vaccine_type = booster_strategy_vaccine_type,
             dose = 8,
             vaccine_mode = case_when(
               vaccine_type %in% c("Moderna","Pfizer") ~ 'mRNA',
               vaccine_type %in% c("AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac") ~ 'viral'),
             coverage_this_date = NA #shouldn't be used anyway
      ) %>% 
      select(date,vaccine_type,vaccine_mode,dose,coverage_this_date,doses_delivered_this_date,age_group,FROM_dose,FROM_vaccine_type)
    
    #CHECK #########################################
    #Checking all doses delivered
    if (round(sum(booster_delivery_outline$doses_delivered_this_date)) != round(sum(eligible_pop$doses_delivered))) { #if not all doses delivered
      if (is.na(restriction_date) == TRUE){ #if no restriction date -> error
        stop('not all booster doses delivered')
      } else{ #else if restriction date
        if (round(as.numeric(restriction_date - booster_strategy_start_date +1) * booster_rollout_speed) == round(sum(workshop$doses))){
          #if restriction date causing non delivery of doses
        } else{  
          stop('not all booster doses delivered')
        }
      } 
    }
    
    #################################################
    ggplot(booster_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(FROM_vaccine_type)))

  if (booster_prioritised == 'Y'){
    return(booster_delivery_outline)
  } else{
    vaccination_history_MODF = bind_rows(vaccination_history_FINAL_local,booster_delivery_outline)
    return(vaccination_history_MODF)
  }

}
#_______________________________________________________________________________