### This function allocates future doses by age_group and dose strategy


###### Coding vaccine prioritisation strategies
#options(scipen = 100) #removes scientific notation
#options(warn = 0) # = 2 when debugging

vax_strategy <- function(vax_strategy_start_date,       # start of hypothetical vaccination program
                         vax_strategy_num_doses,        # num of doses avaliable
                         vax_strategy_roll_out_speed,   # doses delivered per day
                         vax_delivery_group = 'universal', #options = 'universal','at_risk','general_public'. NOTE: 'universal' assumes no risk groups and sets to 'general_public'
                         vax_age_strategy,              # options: "oldest", "youngest","50_down","uniform"
                         vax_dose_strategy,             # number of doses delivered per person
                         vax_strategy_vaccine_type,     # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"  
                         vax_strategy_vaccine_interval, # (days) c() of interval between doses
                         vax_strategy_max_expected_cov,  # value between 0-1 (equivalent to %) of age group willing to be vaccinated
                         restriction_date = NA,
                         vax_roll_out_speed_modifier = data.frame(),
                         vax_end_hist   = vaccine_coverage_end_history
){
  
  ##### SETUP ###################################################################
  ### WARNINGS 
  if (vax_strategy_start_date <= max(vaccination_history_TRUE$date)){ 
    warning ('Your hypothetical vaccine campaign start date needs to be in the future!')
  }
  if (!(vax_strategy_vaccine_type %in% c("Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"))){
    stop('pick a valid vaccine type, or check your spelling!')
  }
  if (vax_delivery_group != "universal" & num_risk_groups == 1){
    warning("You need a risk group to have a risk strategy! We have overwritten vax_delivery_group = 'universal'")
    vax_delivery_group = 'universal'
  }
  if (vax_dose_strategy - 1 > length(vax_strategy_vaccine_interval)){ 
    stop('Please specify the interval between each vaccine dose!')
  }
  
  ### IMPORTS
  vax_strategy_vaccine_interval = round(vax_strategy_vaccine_interval)
  
  prioritisation_csv <- read.csv("1_inputs/prioritisation_strategies.csv",header=TRUE)
  
  if (vax_delivery_group %in% c('universal','general_public')){ this_risk_group = 'general_public'
  } else if (vax_delivery_group == 'at_risk'){ this_risk_group = risk_group_name}
  
  if (vax_delivery_group == 'universal'){ this_pop = data.frame(pop_setting)
  } else if (vax_delivery_group == 'at_risk'){ this_pop = pop_risk_group_dn[pop_risk_group_dn$risk_group == risk_group_name,] %>% select(age_group,pop)
  } else if (vax_delivery_group == 'general_public'){ this_pop = pop_risk_group_dn[pop_risk_group_dn$risk_group == 'general_public',] %>% select(age_group,pop)}
  
  ### IS THIS A BOOSTER?
  booster_dose = "N"
  vax_proportion_booster = 0
  
  if (vax_dose_strategy == 3){booster_dose = "Y"}
  if (vax_dose_strategy == 2 & vax_strategy_vaccine_type == "Johnson & Johnson"){booster_dose = "Y"}
  
  if (booster_dose == "Y"){
    booster_dose_number = vax_dose_strategy
    num_vax_doses = booster_dose_number
    booster_dose_interval = round(vax_strategy_vaccine_interval[vax_dose_strategy-1])
    if (booster_dose_interval < 60){warning('Are you sure that there is less than 2 months between primary and booster doses?')}
    
    ### IF A BOOSTER WE NEED TO CONSIDER INDIVIDUALS WHO HAVE ALREADY RECIEVED PRIMARY DOSES
    #what proportion of the 'at risk' population are untouched by the existing vaccination program?
    workshop_cov = vaccination_history_TRUE %>% # current vaccine coverage in at risk
      filter(date == max(vaccination_history_TRUE$date) &
               risk_group == this_risk_group)
    
    # pop touched by existing vaccination program
    workshop_touched = workshop_cov %>% 
      filter(dose == 1) %>%
      group_by(age_group) %>%
      summarise(cov = sum(coverage_this_date))

    
    #relevant age groups
    eligible_age_groups = prioritisation_csv %>%
      filter(strategy == vax_age_strategy & priority != 99)
    risk_group_age_groups = pop_risk_group_dn %>%
      filter(risk_group == this_risk_group & pop>0)
    
    #consideration of 'unreachable' % either for vaccine hesitancy or access
    unreachable = 1-vax_strategy_max_expected_cov 
    workshop_pop_dn = workshop_touched %>% 
      filter(age_group %in% eligible_age_groups$age_group & age_group %in% risk_group_age_groups$age_group) %>%
      left_join(this_pop, by = "age_group") %>%
      mutate(pop_touched = pop*cov,
             pop_untouched = pop - pop*unreachable - pop*cov)
    
    # calculate proportion to booster / (booster + primary)
    vax_proportion_booster = (sum(workshop_pop_dn$pop_touched))/(sum(workshop_pop_dn$pop_untouched)*vax_dose_strategy+ sum(workshop_pop_dn$pop_touched))  
    #NOTE: we are prioritising doses NOT individuals to ensure that the risk group finish their primary schedule before or at the same time as the general public
    #The code has also been crafted to inflate pop_untouched for a double-dose vaccine
    
    #CHECK
    check = workshop_pop_dn %>% filter(round(pop*vax_strategy_max_expected_cov) != round(pop_touched + pop_untouched))
    if (nrow(check)>0){
      stop('full population not considered')
    }
  }
  #_______________________________________________________________________________
  
  
  
  ##### PRIMARY DOSES ###################################################################################################################################
  #####(1/4) Calculate the eligible population ###################################
  primary_rollout_speed = vax_strategy_roll_out_speed * (1-vax_proportion_booster)
  primary_dose_allocation = round(vax_strategy_num_doses* (1-vax_proportion_booster*(1/vax_dose_strategy))) #because more doses needed for people of primary schedule
  if (nrow(vax_roll_out_speed_modifier)>0){
    primary_speed_modifier =  vax_roll_out_speed_modifier %>% mutate(doses_avaliable = doses_avaliable * (1-vax_proportion_booster),cumsum = cumsum * (1-vax_proportion_booster))
  }

  eligible_pop = this_pop
  colnames(eligible_pop) = c('age_group','eligible_individuals') 
  
  #make long by dose
  workshop = data.frame()
  for (i in 1:num_vax_doses){
    this_dose = as.data.frame(eligible_pop) %>% mutate(dose = i)
    workshop = rbind(workshop,this_dose)
  }
  eligible_pop= workshop
  
  #remove already vaccinated individuals
  existing_coverage = crossing(dose = c(1:num_vax_doses),
                               age_group = age_group_labels,
                               cov_to_date = 0)
  
  for (d in 1:num_vax_doses){
    for (i in 1:num_age_groups){
      #need to sum across vaccine_coverage (as this is vaccination_history_POP split across age groups)
      existing_coverage$cov_to_date[existing_coverage$dose == d & 
                                      existing_coverage$age_group == age_group_labels[i]] = 
        sum(vax_end_hist$coverage_this_date[vax_end_hist$risk_group == this_risk_group &
                                              vax_end_hist$dose == d & 
                                              vax_end_hist$age_group == age_group_labels[i]])
    }
  } 
  
  if (max(existing_coverage$cov_to_date) > vax_strategy_max_expected_cov ){stop('more vaccines already delivered than max expected coverage')}
  
  
  #now remove vaccinated, and vaccine hesistant
  unreachable = 1-vax_strategy_max_expected_cov
  eligible_pop <- eligible_pop %>% 
    left_join(existing_coverage, by = c("age_group", "dose")) %>%
    mutate(eligible_individuals = eligible_individuals *(1-(cov_to_date+unreachable))) %>%
    select(age_group,dose,eligible_individuals)
  
  for (d in 2:num_vax_doses){
    eligible_pop$eligible_individuals[eligible_pop$dose == d] = eligible_pop$eligible_individuals[eligible_pop$dose == 1] 
  }
  #_______________________________________________________________________________
  
  
  
  #####(2/4) Place priority # on age group by strategy ############################
  if (vax_age_strategy %in% unique(prioritisation_csv$strategy)) {
    elected_strategy = prioritisation_csv[prioritisation_csv$strategy == vax_age_strategy,c('age_group','priority')]
    eligible_pop <- eligible_pop %>% 
      left_join(elected_strategy, by = "age_group")
  } else if (vax_age_strategy == "manual_overwrite"){
    eligible_pop <- eligible_pop %>%
      mutate(priority= case_when(
        age_group == '70 to 100' ~ 7,
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
  doses_to_deliver = primary_dose_allocation
  priority_num = 1
  n=length(unique(eligible_pop$priority))
  highest_priority = sort(unique(eligible_pop$priority),partial=n-1)[n-1] #highest valid priority
  
  while (doses_to_deliver>0 & priority_num <= (highest_priority)){
    priority_group = eligible_pop[eligible_pop$priority == priority_num,]
    
    #check enough for all priority groups
    workshop = doses_to_deliver/vax_dose_strategy - sum(priority_group$eligible_individuals[priority_group$dose == 1]) 

    if (workshop < 0){
      prop_to_deliver = doses_to_deliver/vax_dose_strategy / sum(priority_group$eligible_individuals[priority_group$dose == 1])
      
      priority_group$doses_delivered[priority_group$dose == vax_dose_strategy] = prop_to_deliver * priority_group$eligible_individuals[priority_group$dose == vax_dose_strategy]
      if (vax_dose_strategy > 1){
        for (dose in 1:(vax_dose_strategy-1)){
          priority_group$doses_delivered[priority_group$dose == dose] = prop_to_deliver * priority_group$eligible_individuals[priority_group$dose == dose] 
        }
      }
      doses_to_deliver = 0
      
    } else if (workshop>0){
      priority_group$doses_delivered[priority_group$dose == vax_dose_strategy] =  priority_group$eligible_individuals[priority_group$dose == vax_dose_strategy]
      doses_to_deliver = doses_to_deliver - sum(priority_group$eligible_individuals[priority_group$dose == vax_dose_strategy])
      
      if (vax_dose_strategy > 1){
        for (dose in 1:(vax_dose_strategy-1)){
          workshop2 = doses_to_deliver- sum(priority_group$eligible_individuals[priority_group$dose == dose])
          if(workshop2 < 0){
            prop_to_deliver = doses_to_deliver / sum(priority_group$eligible_individuals[priority_group$dose == dose])
            priority_group$doses_delivered[priority_group$dose == dose] = doses_to_deliver
            doses_to_deliver = 0
          } else{
            priority_group$doses_delivered[priority_group$dose == dose] = priority_group$eligible_individuals[priority_group$dose == dose]
            doses_to_deliver = doses_to_deliver - sum(priority_group$eligible_individuals[priority_group$dose == dose])
          }
        }
      }
      
    }
    eligible_pop$doses_delivered[eligible_pop$priority == priority_num] = priority_group$doses_delivered
    priority_num = priority_num + 1  
  }
  #_______________________________________________________________________________
  
  
  
  #####(4/4) Distribute between days #############################################
  # we should use:
  # (1) primary_dose_allocation - doses to deliver
  # (2) primary_rollout_speed - max doses delivered per day
  # (3) vax_strategy_start_date - first day of doses delivered
  
  VA =  eligible_pop %>% mutate(doses_left = doses_delivered)
  priority_num = 1
  priority_group  = as.character(unique(VA$age_group[VA$priority == priority_num]))
  
  ceiling = min(sum(eligible_pop$doses_delivered),primary_dose_allocation) #max delivery of doses (either limited by eligible individuals, or available doses!)
  
  if (vax_dose_strategy == 1){
    timeframe = ceiling/primary_rollout_speed
    daily_per_dose = primary_rollout_speed
    in_interval_flag = "NA"
  } else if (vax_dose_strategy > 1){
    if((ceiling/vax_dose_strategy)/primary_rollout_speed<vax_strategy_vaccine_interval[1]){ #if interval longer than length taken to rollout doses
      timeframe = (ceiling/vax_dose_strategy)/primary_rollout_speed
      daily_per_dose = primary_rollout_speed
      in_interval_flag = "Y"
    } else{
      timeframe = ceiling/(primary_rollout_speed)
      daily_per_dose = primary_rollout_speed/vax_dose_strategy
      in_interval_flag = "N"
    }
  }
  
  #modify if end in sight of restriction date
  if (is.na(restriction_date) == FALSE){
    timeframe = as.numeric(restriction_date - vax_strategy_start_date) + 1
  }
  
  #modify if available doses restricted
  if (nrow(vax_roll_out_speed_modifier)>0){
    exception_list = unique(primary_speed_modifier$day)
    if (in_interval_flag == "Y"){workshop = primary_speed_modifier %>% filter(cumsum<(ceiling/vax_dose_strategy))
    } else{workshop = primary_speed_modifier %>% filter(cumsum<ceiling)}
    workshop = nrow(workshop) + 1
    if (max(primary_speed_modifier$cumsum)+daily_per_dose>(ceiling)){
      timeframe = workshop
    } else{
      extension_time = (ceiling - max(primary_speed_modifier$cumsum))/(daily_per_dose*vax_dose_strategy)
      timeframe = workshop + extension_time
    }
  } else{
    exception_list = c(-1)
  }
  
  timeframe = ceiling(timeframe)
  
  length_track = timeframe
  if (vax_dose_strategy > 1){
    for (add in 1:(vax_dose_strategy-1)){
      length_track=length_track+vax_strategy_vaccine_interval[add]
    }
  }
  
  vax_delivery_outline <- crossing(day = c(1:length_track),
                                   dose = c(1:num_vax_doses),
                                   age_group = age_group_labels,
                                   doses_delivered = c(0))
  
  for (day in 1:timeframe){

    avaliable = daily_per_dose
    daily_per_dose_here = daily_per_dose
    
    if (day %in% exception_list){
      avaliable = primary_speed_modifier$doses_avaliable[primary_speed_modifier$day == day]
      if (in_interval_flag == "N"){avaliable = avaliable/vax_dose_strategy}
      daily_per_dose_here = avaliable 
    }
    
    #ensuring that we don't overshoot available doses
    if (day == timeframe){
      workshop_leftover = ceiling - sum(vax_delivery_outline$doses_delivered)
      avaliable = min(workshop_leftover/vax_dose_strategy,daily_per_dose)
      if (nrow(vax_roll_out_speed_modifier) == 0){
        if(round(workshop_leftover) != round(ceiling - (timeframe-1)*daily_per_dose*vax_dose_strategy)){stop('doses not delivered as expected!')}
      }
    }
    if (avaliable > sum(VA$doses_left)){avaliable = sum(VA$doses_left)}
    
    while(avaliable > 0 & priority_num <= highest_priority){
      
      if(round(sum(VA$doses_left[VA$priority == priority_num]),digits=4)>0){ 
        #i.e., while we still have doses to deliver in this priority group
        
        if(0 %in% VA$doses_left[VA$priority == priority_num & VA$dose == 1]){
          age_complete = VA$age_group[VA$doses_left == 0 & VA$priority == priority_num & VA$dose == 1]
          VA$priority[VA$age_group %in% age_complete] = VA$priority[VA$age_group %in% age_complete] + 100
          
          priority_group = as.character(unique(VA$age_group[VA$priority == priority_num]))
        } #FIX - when one age group in the priority group runs out first
        
        workshop_doses = min(sum(VA$doses_left[VA$priority == priority_num & VA$dose == 1]),
                             daily_per_dose_here,
                             avaliable)
        #either deliver max capacity or number left in this group, whichever is fewer
        
        leftover=0
        VA_pt = VA #snapshot
        
        for (i in length(priority_group):1){
          workshop_age = priority_group[i]
          workshop_prop = sum(VA_pt$doses_left[VA_pt$age_group == workshop_age])/sum(VA_pt$doses_left[VA_pt$priority == priority_num])
          workshop_calc = workshop_doses * workshop_prop + leftover
          
          if (workshop_calc > VA$doses_left[VA$age_group == workshop_age & VA$dose == 1]){
            leftover = abs(workshop_calc - VA$doses_left[VA$age_group == workshop_age & VA$dose == 1])
            workshop_calc = VA$doses_left[VA$age_group == workshop_age & VA$dose == 1]
            
          } else{
            leftover = 0
          }
          vax_delivery_outline$doses_delivered[vax_delivery_outline$day == day &
                                                 vax_delivery_outline$dose == 1 &
                                                 vax_delivery_outline$age_group == workshop_age] = workshop_calc
          VA$doses_left[VA$age_group == workshop_age & VA$dose == 1] = VA$doses_left[VA$age_group == workshop_age & VA$dose == 1] - workshop_calc
          
          if (vax_dose_strategy > 1){
            for (dose in 2:vax_dose_strategy){
              vax_delivery_outline$doses_delivered[vax_delivery_outline$day == day + vax_strategy_vaccine_interval[dose-1] &
                                                     vax_delivery_outline$dose == dose &
                                                     vax_delivery_outline$age_group == workshop_age] = workshop_calc  
              VA$doses_left[VA$age_group == workshop_age & VA$dose == dose] = VA$doses_left[VA$age_group == workshop_age & VA$dose == dose] - workshop_calc
            }
          }
          
          
        }
        
        avaliable = min(avaliable - workshop_doses,ceiling - sum(vax_delivery_outline$doses_delivered))
        
        
      } else if (round(sum(VA$doses_left[VA$priority == priority_num]),digits=4)==0){
        priority_num = priority_num+1
        priority_group = as.character(unique(VA$age_group[VA$priority == priority_num]))
      } else{
        stop('negative doses left, reconsider!')
      }
    } #<end while loop>
    #if(!round(sum(vax_delivery_outline$doses_delivered[vax_delivery_outline$day == day])) == primary_rollout_speed ){warning(paste('issue on day',day))}
  }
  
  ### formatting vax_delivery_outline to align with vaccination_history_TRUE
  vax_delivery_outline = vax_delivery_outline %>%
    mutate(date = vax_strategy_start_date + (day-1),
           vaccine_type = vax_strategy_vaccine_type,
           vaccine_mode = case_when(
             vaccine_type %in% c("Moderna","Pfizer") ~ 'mRNA',
             vaccine_type %in% c("AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac") ~ 'viral'
           ),
           coverage_this_date = NA #shouldn't be used anyway
           ) %>%
    rename(doses_delivered_this_date = doses_delivered) %>% 
    select(date,vaccine_type,vaccine_mode,dose,coverage_this_date,doses_delivered_this_date,age_group)
  
  #CHECK 
  if (round(sum(vax_delivery_outline$doses_delivered_this_date)) != round(sum(eligible_pop$doses_delivered))) { #if not all doses delivered
    if (is.na(restriction_date) == TRUE){ #if no restriction date -> error
      stop('not all doses delivered!')
    } else{ #else if restriction date
      if (round(as.numeric(restriction_date - vax_strategy_start_date +1) * primary_rollout_speed) == round(sum(workshop$doses))){
        #if restriction date causing non delivery of doses
      } else{  
        stop('not all doses delivered!')
      }
    } 
  }
  
  ggplot(vax_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  #Don't forget the booster doses to completely unvaccinated individuals appear here
  ######################################################################################################################################################
  
  
 
  
  
  ##### BOOSTER DOSES ###################################################################################################################################
  if (booster_dose == "Y"){
    #poached from other primary doses (other boosters are included above)
    
    booster_rollout_speed = vax_strategy_roll_out_speed * vax_proportion_booster
    booster_dose_allocation = round(vax_strategy_num_doses* vax_proportion_booster*(1/vax_dose_strategy))
    if (nrow(vax_roll_out_speed_modifier)>0){
      booster_speed_modifier =  vax_roll_out_speed_modifier %>% mutate(doses_avaliable = doses_avaliable * vax_proportion_booster,cumsum = cumsum * vax_proportion_booster)
    } else{
      booster_speed_modifier = crossing(date = unique(vax_delivery_outline$date),
                                        doses_avaliable = booster_rollout_speed) %>%
        mutate(cumsum = cumsum(doses_avaliable))
    }
    
    check_elig = vaccination_history_TRUE %>% #NB: using TRUE because hypoth proj should already be boosted
      mutate(primary_schedule_complete = case_when(
        vaccine_type == "Johnson & Johnson" ~  'Y',
        dose %in% c(2,3) ~ 'Y',
        TRUE ~ 'N'),
        boosted = case_when(
          vaccine_type == "Johnson & Johnson" & dose == 2 ~ 'Y',
          dose == 3 ~ 'Y',
          TRUE ~ 'N'))  %>%
      filter(risk_group %in% this_risk_group &
               # primary_schedule_complete == "Y" & #CAN CHANGE THIS TOGGLE #ASSUMPTION / COMEBACK - including those with only 1 primary dose
               boosted == "N" &
               coverage_this_date > 0) %>%
      mutate(date = date + booster_dose_interval) %>% 
      ungroup() %>% group_by(date) %>% summarise(total_doses = sum(doses_delivered_this_date),.groups = "keep") %>%
      ungroup() %>% mutate(individuals_eligible = cumsum(total_doses)) %>%
      select(-total_doses)

    booster_speed_modifier = booster_speed_modifier %>% 
      left_join(check_elig, by = "date")
    
    #make sure that individuals are eligible when boosted
    check_less_than = booster_speed_modifier  %>%
      filter(cumsum > individuals_eligible) 
    
    while(nrow(check_less_than)>0){
      booster_speed_modifier = booster_speed_modifier %>% 
         mutate(difference = cumsum - individuals_eligible)
      this_date = check_less_than$date[1]
      booster_speed_modifier$doses_avaliable[booster_speed_modifier$date == this_date] = 
        booster_speed_modifier$doses_avaliable[booster_speed_modifier$date == this_date] - booster_speed_modifier$difference[booster_speed_modifier$date == this_date]
      booster_speed_modifier = booster_speed_modifier %>% 
        mutate(cumsum = cumsum(doses_avaliable))
      check_less_than = booster_speed_modifier %>% filter(cumsum>individuals_eligible)
    }
    ###
    
    
    #####(1/4) Calculate the eligible population ###################################
    eligible_pop = this_pop
    colnames(eligible_pop) = c('age_group','eligible_individuals') 
    
    existing_coverage = crossing(dose = unique(vaccination_history_TRUE$dose),
                                 vaccine_type = unique(vaccination_history_TRUE$vaccine_type),
                                 age_group = age_group_labels,
                                 cov_to_date = 0)
    
    for (t in 1:length(unique(vaccination_history_TRUE$vaccine_type))){
      for (d in 1:length(unique(vaccination_history_TRUE$dose))){
        for (i in 1:num_age_groups){
          #need to sum across vaccine_coverage (as this is vaccination_history_POP split across age groups)
          existing_coverage$cov_to_date[existing_coverage$dose == d & 
                                          existing_coverage$age_group == age_group_labels[i] &
                                          existing_coverage$vaccine_type == unique(vaccination_history_TRUE$vaccine_type)[t]] = 
            vax_end_hist$coverage_this_date[vax_end_hist$risk_group == this_risk_group &
                                              vax_end_hist$dose == d & 
                                              vax_end_hist$age_group == age_group_labels[i]&
                                              vax_end_hist$vaccine_type == unique(vaccination_history_TRUE$vaccine_type)[t]]
        }
      } 
    }
    
    eligible_pop <- eligible_pop %>% 
      left_join(existing_coverage, by = "age_group") %>%
      mutate(eligible_individuals = eligible_individuals * cov_to_date) %>%
      select(age_group,dose,vaccine_type,eligible_individuals)
    
    #remove double counted tidy
    for (d in 2:num_vax_doses){
      remove = eligible_pop %>% 
        filter(dose == d) %>%
        rename(complete_vax = eligible_individuals) %>%
        select(-dose)
      
      if (nrow(remove)>0){
        eligible_pop = eligible_pop %>% 
          left_join(remove, by = c('age_group','vaccine_type')) %>%
          mutate(eligible_individuals = case_when(
            dose == (d-1) & complete_vax > eligible_individuals ~ 0, #this shouldn't be triggered
            dose == (d-1) ~ eligible_individuals - complete_vax,
            TRUE ~ eligible_individuals,
          )) %>%
          select(-complete_vax)
      }
    }
    #_______________________________________________________________________________
    
    
    
    #####(2/4) Place priority # on age group by strategy ############################
    if (vax_age_strategy %in% unique(prioritisation_csv$strategy)) {
      elected_strategy = prioritisation_csv[prioritisation_csv$strategy == vax_age_strategy,c('age_group','priority')]
      eligible_pop <- eligible_pop %>% 
        left_join(elected_strategy, by = "age_group")
    } else if (vax_age_strategy == "manual_overwrite"){
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
    n=length(unique(eligible_pop$priority))
    highest_priority = sort(unique(eligible_pop$priority),partial=n-1)[n-1] #highest valid priority
    
    while (doses_to_deliver>0 & priority_num <= (highest_priority)){
      priority_group = eligible_pop[eligible_pop$priority == priority_num,]
      
      workshop = doses_to_deliver - sum(priority_group$eligible_individuals) #check enough for all priority groups
      
      if (workshop < 0){
        prop_to_deliver = doses_to_deliver / sum(priority_group$eligible_individuals)
        priority_group$doses_delivered = prop_to_deliver * priority_group$eligible_individuals
        doses_to_deliver = 0
      } else if (workshop>0){
        priority_group$doses_delivered =  priority_group$eligible_individuals
        doses_to_deliver = doses_to_deliver - sum(priority_group$eligible_individuals)
      }
      
      eligible_pop$doses_delivered[eligible_pop$priority == priority_num] = priority_group$doses_delivered
      
      priority_num = priority_num + 1  
    }
    if (nrow(eligible_pop[eligible_pop$doses_delivered>eligible_pop$eligible_individuals,])>0){stop('more doses delivered than eligible individuals!')}
    #_______________________________________________________________________________
    
    
    
    #####(4/4) Distribute between days #############################################
    # we should use:
    # (1) booster_dose_allocation - doses to deliver
    # (2) booster_rollout_speed - max doses delivered per day
    # (3) vax_strategy_start_date - first day of doses delivered
    
    VA =  eligible_pop %>% mutate(doses_left = doses_delivered)
    priority_num = 1
    priority_group  = as.character(unique(VA$age_group[VA$priority == priority_num]))
    
    ceiling = min(sum(eligible_pop$doses_delivered),booster_dose_allocation)
    timeframe = ceiling/booster_rollout_speed
    daily_per_dose = booster_rollout_speed
    
    #modify if available doses restricted
    booster_speed_modifier$day = booster_speed_modifier$date - vax_strategy_start_date + 1
    if (nrow(booster_speed_modifier)>0){
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
                
                workshop_age = priority_group[i]
                workshop_type = unique(VA$vaccine_type)[t]
                
                workshop_prop = VA_pt$doses_left[VA_pt$age_group == workshop_age & VA_pt$dose == d & VA_pt$vaccine_type == workshop_type]/
                  sum(VA_pt$doses_left[VA_pt$priority == priority_num])
                workshop_calc = workshop_doses * workshop_prop + leftover
                
                if (workshop_calc > VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type]){
                  leftover = abs(workshop_calc - VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type])
                  workshop_calc = VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type]
                  
                } else{
                  leftover = 0
                }
                booster_delivery_outline$doses_delivered[booster_delivery_outline$day == day &
                                                           booster_delivery_outline$dose == d &
                                                           booster_delivery_outline$vaccine_type == workshop_type &
                                                           booster_delivery_outline$age_group == workshop_age] = workshop_calc
                VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type] = 
                  VA$doses_left[VA$age_group == workshop_age & VA$dose == d & VA$vaccine_type == workshop_type] - workshop_calc
                
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
      mutate(date = vax_strategy_start_date + (day-1)) %>%
      rename(FROM_vaccine_type = vaccine_type,
             FROM_dose = dose,
             doses_delivered_this_date = doses_delivered) %>%
      mutate(vaccine_type = vax_strategy_vaccine_type,
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
        if (round(as.numeric(restriction_date - vax_strategy_start_date +1) * booster_rollout_speed) == round(sum(workshop$doses))){
          #if restriction date causing non delivery of doses
        } else{  
          stop('not all booster doses delivered')
        }
      } 
    }
    
    #Checking poached doses are not before eligibility
    check_elig = vaccination_history_TRUE %>%
      mutate(primary_schedule_complete = case_when(
        vaccine_type == "Johnson & Johnson" ~  'Y',
        dose %in% c(2,3) ~ 'Y',
        TRUE ~ 'N'),
      boosted = case_when(
        vaccine_type == "Johnson & Johnson" & dose == 2 ~ 'Y',
        dose == 3 ~ 'Y',
        TRUE ~ 'N'))  %>%
      filter(risk_group %in% this_risk_group &
               # primary_schedule_complete == "Y" & #CAN CHANGE THIS TOGGLE
               boosted == "N" &
               coverage_this_date > 0) %>%
      mutate(date = date + booster_dose_interval) %>% 
      ungroup() %>% group_by(date) %>% summarise(total_doses = sum(doses_delivered_this_date),.groups = "keep") %>%
      ungroup() %>% mutate(individuals_eligible = cumsum(total_doses)) %>%
      select(-total_doses)
    
    hypoth_delivery = booster_delivery_outline %>%
      ungroup() %>% group_by(date) %>% summarise(total_doses = sum(doses_delivered_this_date),.groups = "keep") %>%
      ungroup() %>% mutate(cum_doses_delivered = cumsum(total_doses)) %>%
      select(-total_doses)
    
    check_less_than = check_elig %>% 
      left_join(hypoth_delivery, by = "date") %>%
      filter(round(cum_doses_delivered) > round(individuals_eligible))
    
    if (nrow(check_less_than) > 0 ){
      ggplot() +
        geom_line(data = check_elig, aes(x=date,y=individuals_eligible)) +
        geom_point(data = hypoth_delivery, aes(x=date,y=cum_doses_delivered)) +
        xlab('') + ylab('number of individuals eligible for a booster dose') +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_line(color = 'black'))
      stop('more booster doses delivered than individuals eligible!')
    }
    #################################################
    ggplot(booster_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(FROM_vaccine_type)))
    
    vax_delivery_outline = bind_rows(vax_delivery_outline,booster_delivery_outline)
  }
  
  
  if(vax_delivery_group != 'universal'){
    return(vax_delivery_outline)
  } else if(vax_delivery_group == 'universal'){
    vax_delivery_outline$risk_group = this_risk_group
    vaccination_history_MODF = bind_rows(vaccination_history_TRUE,vax_delivery_outline)
    return(vaccination_history_MODF)
  }
}
#_______________________________________________________________________________