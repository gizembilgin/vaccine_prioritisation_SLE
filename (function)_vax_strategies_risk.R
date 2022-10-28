### This function allocates future doses by risk_group using repeated calls of the vax_strategy() function

apply_risk_strategy <- function(
    vax_risk_strategy,             # options: 'Y','N'
    vax_risk_proportion,           # value between 0-1 (equivalent to %) of doses prioritised to the at risk group
    vax_doses_general,             # number of doses delivered to general pop
    vax_doses_risk,                # number of doses delivered to risk group
    risk_group_accessibility = FALSE, #TRUE means opportunistic vaccination possible for this risk group
    risk_group_acceptability = vax_strategy_toggles$vax_strategy_max_expected_cov, #1-vaccine hesistancy
    risk_group_age_broaden = FALSE # TRUE means the expansion of eligibilty to high-risk children
){
  
  ### WARNINGS 
  if (!vax_risk_strategy %in% c('Y','N')){stop('Is the vax strategy on or off? Y/N')}
  if (vax_risk_proportion<0 | vax_risk_proportion>1){stop('vax_risk_proportion must be between 0-1 (0%-100%)')}
  if (vax_risk_proportion == 0 & vax_risk_strategy == 'Y' & !('additional_doses' %in% names(sensitivity_analysis_toggles))){
    warning("Giving 0% priority is not priority! I have overwritten vax_risk_strategy to EQ 'N'")
    vax_risk_strategy = "N"
  }
  #_______________________________________________________________________________
  
  
  
  ### BRANCH ONE: Are we prioritising the at risk group at all?
  if (vax_risk_strategy == "N"){
    #if not prioritised, set vax_risk_proportion to population proportion of risk group
    risk_split = age_risk_split %>% group_by(risk_group) %>% summarise(sum = sum(split),.groups='keep')
    vax_risk_proportion = risk_split$sum[risk_split$risk_group == risk_group_name]
    
    real_doses = vaccination_history_TRUE %>% 
      filter(! age_group %in% c('0 to 4','5 to 9','10 to 17')) %>%
      group_by(risk_group,age_group,dose) %>%
      summarise(doses = sum(doses_delivered_this_date),.groups='keep') %>%
      left_join(pop_risk_group_dn, by = c('risk_group','age_group')) %>%
      mutate(cov=doses/pop) %>%
      arrange(dose,age_group)
    if (length(unique(na.omit(round(real_doses$cov[real_doses$dose == 1],digits=2))))>1 | length(unique(na.omit(round(real_doses$cov[real_doses$dose == 2],digits=2))))>1){
      if (!'vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
        stop('real doses not equal across risk groups')
      }
    }
    
    toggle_equal_priority = "doses" # STATIC TOGGLE: "individuals", "doses"
    if (vax_doses_risk>vax_doses_general){
      if (toggle_equal_priority == "individuals"){
        #vax_risk_proportion remains the same  
      } else if (toggle_equal_priority == "doses"){
        #calculate unvaccinated proportion by 1-coverage first dose
        workshop = vaccine_coverage_end_history %>% 
          filter(dose == 1) %>%
          group_by(age_group,risk_group) %>%
          summarise(unvax = sum(coverage_this_date),.groups='keep') %>%
          mutate(unvax = 1 - unvax)
        
        #(see proof '2022_09_29 vax_risk_proportion (vax_doses_gen != vax_doses_risk).pdf')
        workshop =  pop_risk_group_dn %>%
          left_join(workshop) %>%
          mutate(adult_pop = case_when(
            age_group %in% c('0 to 4','5 to 9','10 to 17') ~ 0,
            TRUE ~ pop)) %>%
          mutate(dose_inflated_pop = case_when(
            risk_group == risk_group_name ~ adult_pop * (1+unvax*vax_doses_general),
            TRUE ~ adult_pop * unvax * vax_doses_general))
        vax_risk_proportion = sum(workshop$dose_inflated_pop[workshop$risk_group == risk_group_name])/sum(workshop$dose_inflated_pop)
      }
    }
  }
  #___________________________________________________________
  
  
  
  ### BRANCH TWO: What % priority are the at risk group receiving? 
  if (risk_group_accessibility == FALSE){
    speed_risk_group_rollout =  vax_strategy_toggles$vax_strategy_roll_out_speed*vax_risk_proportion
  } else if (risk_group_accessibility == TRUE){
    speed_risk_group_rollout = round(sum(pop_high_risk$pop)/(365/12*4.2))*vax_doses_risk #median first antenatal care visit
  }
  if (risk_group_age_broaden == TRUE){
    save_age_strategy = vax_strategy_toggles$vax_age_strategy 
    vax_strategy_toggles$vax_age_strategy = 'uniform'
  }
  
  at_risk_delivery_outline = vax_strategy(vax_delivery_group = 'at_risk',
                                            vax_dose_strategy              = vax_doses_risk,       
                                            vax_strategy_roll_out_speed    = speed_risk_group_rollout,            
                                            vax_strategy_max_expected_cov  = risk_group_acceptability,
                                            vax_strategy_start_date        = vax_strategy_toggles$vax_strategy_start_date,
                                            vax_strategy_num_doses         = vax_strategy_toggles$vax_strategy_num_doses,
                                            vax_age_strategy               = vax_strategy_toggles$vax_age_strategy,            
                                            vax_strategy_vaccine_type      = vax_strategy_toggles$vax_strategy_vaccine_type,            
                                            vax_strategy_vaccine_interval  = vax_strategy_toggles$vax_strategy_vaccine_interval
  )
  at_risk_delivery_outline = at_risk_delivery_outline %>% mutate(risk_group = risk_group_name)
   
  ggplot(at_risk_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  sum(at_risk_delivery_outline$doses_delivered_this_date)
  total = at_risk_delivery_outline %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep")
  ggplot(total) + geom_point(aes(x=date,y=doses_delivered_this_date))
  
  if (risk_group_age_broaden == TRUE){ #reset for delivery to general public
    vax_strategy_toggles$vax_age_strategy = save_age_strategy
  }
  #___________________________________________________________
  
  
  
  ### BRANCH THREE: Distribute remaining doses to general population
  #<interim update vaccine_coverage_end_history so eligible pop is correct>
  vaccination_history_MODF = rbind(vaccination_history_TRUE,at_risk_delivery_outline)
  
  vaccination_history_MODF = vaccination_history_MODF %>% 
    left_join(pop_risk_group_dn, by = c("age_group", "risk_group")) %>%
    group_by(risk_group,age_group,vaccine_type,dose) %>%
    mutate(coverage_this_date = cumsum(doses_delivered_this_date)/pop) 
  
  vaccine_coverage_end_history_UPDATED = data.frame()
  for (risk in 1:num_risk_groups){
    for (t in 1:length(unique(vaccination_history_MODF$vaccine_type))){
      this_vax = vaccination_history_MODF %>%
        filter(vaccine_type == unique(vaccination_history_MODF$vaccine_type)[t] &
                 risk_group == risk_group_labels[risk])
      this_vax = this_vax %>% filter(date == max(this_vax$date)) %>%
        select(dose,vaccine_type,age_group,risk_group,coverage_this_date)
      vaccine_coverage_end_history_UPDATED = rbind(vaccine_coverage_end_history_UPDATED,this_vax)
    }
  }
  vaccine_coverage_end_history_UPDATED$coverage_this_date[is.na(vaccine_coverage_end_history_UPDATED$coverage_this_date)] = 0
  #<fin>
  
  #<interim calculate doses available per day>
  limiter = data.frame()
  leftover_doses = vax_strategy_toggles$vax_strategy_num_doses - sum(at_risk_delivery_outline$doses_delivered_this_date)
  if (risk_group_accessibility == FALSE){
    limiter = at_risk_delivery_outline %>% group_by(date) %>% 
      summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep") 
    # filter(doses_delivered_this_date > 0 ) %>%
    limiter = limiter  %>% ungroup() %>% 
      mutate(doses_avaliable = vax_strategy_toggles$vax_strategy_roll_out_speed - doses_delivered_this_date,
                                 day = as.numeric(date - min(limiter$date) + 1),
                                 cumsum = cumsum(doses_avaliable)) #see line 315 in function vax strategy
    if (vax_doses_general>1){
      #check that second dose of general availability is = to first dose of general availability
      for (d in 2:vax_doses_general){
        workshop = limiter %>% mutate(day = day - vax_strategy_toggles$vax_strategy_vaccine_interval[d-1]) %>%
          rename(next_dose_avaliab = doses_avaliable) %>%
          select(day,next_dose_avaliab) 
        limiter = limiter %>% left_join(workshop, by = "day")
        limiter = limiter %>% ungroup() %>% 
          mutate(doses_avaliable = case_when(
          doses_avaliable>next_dose_avaliab ~ next_dose_avaliab,
          TRUE ~ doses_avaliable)) %>%
          mutate( cumsum = cumsum(doses_avaliable))
      }
    }
  } 
  #<fin>
  
  generalPublic_leftover_outline = vax_strategy(vax_delivery_group             = 'general_public',
                                                vax_dose_strategy              = vax_doses_general,       
                                                vax_strategy_roll_out_speed    = vax_strategy_toggles$vax_strategy_roll_out_speed,
                                                vax_roll_out_speed_modifier    = limiter,
                                                vax_strategy_num_doses         = leftover_doses ,
                                                vax_strategy_start_date        = vax_strategy_toggles$vax_strategy_start_date,
                                                vax_age_strategy               = vax_strategy_toggles$vax_age_strategy,            
                                                vax_strategy_vaccine_type      = vax_strategy_toggles$vax_strategy_vaccine_type,            
                                                vax_strategy_vaccine_interval  = vax_strategy_toggles$vax_strategy_vaccine_interval,            
                                                vax_strategy_max_expected_cov  = vax_strategy_toggles$vax_strategy_max_expected_cov,
                                                vax_end_hist   = vaccine_coverage_end_history_UPDATED
  )
  generalPublic_leftover_outline = generalPublic_leftover_outline %>% mutate(risk_group = 'general_public')
  
  ggplot(generalPublic_leftover_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  sum(generalPublic_leftover_outline$doses_delivered_this_date)
  total = generalPublic_leftover_outline %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep")
  ggplot(total) + geom_point(aes(x=date,y=doses_delivered_this_date))  
  
  vaccination_history_MODF = bind_rows(vaccination_history_TRUE,at_risk_delivery_outline,generalPublic_leftover_outline)

  ###CHECKS
  #CHECK 1: total doses delivered <= total doses available
  if (risk_group_accessibility == FALSE){
    check_df = bind_rows(at_risk_delivery_outline,generalPublic_leftover_outline)
    if (round(sum(check_df$doses_delivered_this_date)) > vax_strategy_toggles$vax_strategy_num_doses){stop('Total doses delivered > total doses avaliable!')}
    
    check_df_daily = check_df %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label ='total')
    if (nrow(check_df_daily[round(check_df_daily$doses_delivered_this_date)>vax_strategy_toggles$vax_strategy_roll_out_speed,])>0){
      ggplot(check_df_daily) + geom_point(aes(x=date,y=doses_delivered_this_date))
      
      df1 = at_risk_delivery_outline  %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label = 'at risk delivery')
      df2 = generalPublic_leftover_outline  %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label = 'general public leftover')
      check_df_daily_comp = bind_rows(df1,df2,check_df_daily)
      ggplot(check_df_daily_comp) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(label)))
      ggplot(check_df_daily_comp[check_df_daily_comp$label != 'total',]) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(label)))
      
      nrow_greater = nrow(check_df_daily[round(check_df_daily$doses_delivered_this_date)>vax_strategy_toggles$vax_strategy_roll_out_speed,])
      
      warning('Caution, more doses delivered per day than capacity for',paste(nrow_greater),' days')
      
      if (nrow_greater>10){
        stop('Error maximum daily doses exceeded for more than 10 days')
      }
      if ((max(check_df_daily$doses_delivered_this_date)-vax_strategy_toggles$vax_strategy_roll_out_speed )/vax_strategy_toggles$vax_strategy_roll_out_speed > 0.1){
        stop('Error maximum daily doses exceeded by 10%') #NOTE - 10% leeway, unless want more complicated limiter on booster dose
      }
      
    }
    ggplot(check_df[check_df$risk_group == 'general_public',]) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
    ggplot(check_df[check_df$risk_group != 'general_public',]) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  }
  #CHECKS: for non-prioritised dose delivery with risk groups
  if (vax_risk_strategy == "N"){
    hypoth_doses = vaccination_history_MODF %>% 
      filter(! age_group %in% c('0 to 4','5 to 9','10 to 17')) %>%
      mutate(dose = case_when(
        dose == 8 ~ 3,
        dose == 2 & vaccine_type == 'Johnson & Johnson' ~ 3,
        TRUE ~dose
      )) %>%
      group_by(risk_group,age_group,dose) %>%
      summarise(doses = sum(doses_delivered_this_date),.groups = "keep") %>%
      left_join(pop_risk_group_dn, by = c("risk_group", "age_group")) %>%
      mutate(cov=doses/pop) %>%
      arrange(dose,age_group)
    
    expected_cov = unique(c(vax_strategy_toggles$vax_strategy_max_expected_cov,risk_group_acceptability))
    if (length(unique(na.omit(round(hypoth_doses$cov[hypoth_doses$dose == 1],digits=2)))) != length(expected_cov) |
        unique(na.omit(round(hypoth_doses$cov[hypoth_doses$dose == 1],digits=2)))[1] != expected_cov[1]){
      warning('not all who are willing have recieved the first dose')
    }
    if (length(unique(na.omit(round(hypoth_doses$cov[hypoth_doses$dose == 1],digits=2))))>1){
      if (!'vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
        stop('hypoth dose one not equal across risk groups')
      }
    } 
    if (vax_doses_risk==vax_doses_general & !'vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
      if (max(vaccination_history_TRUE$date[vaccination_history_TRUE$risk_group == 'general_public']) !=
          max(vaccination_history_TRUE$date[vaccination_history_TRUE$risk_group == risk_group_name])){
        stop('Existing rollout doesnt align between risk groups')
      }
      if (max(vaccination_history_MODF$date[vaccination_history_MODF$risk_group == 'general_public']) !=
            max(vaccination_history_MODF$date[vaccination_history_MODF$risk_group == risk_group_name])){
        if (vax_doses_risk==1){stop('max delivery dates dont align between risk groups (d=1)')
        } else if (abs(max(vaccination_history_MODF$date[vaccination_history_MODF$risk_group == 'general_public']) -
                        max(vaccination_history_MODF$date[vaccination_history_MODF$risk_group == risk_group_name])) > 60){
          warning('max delivery dates dont align between risk groups (d>1)')
        }
      }
    }
    if (vax_doses_general >= 2){
      if (unique(na.omit(round(hypoth_doses$cov[hypoth_doses$dose == 3 & hypoth_doses$risk_group == 'general_public'],digits=2))) != vax_strategy_toggles$vax_strategy_max_expected_cov){
        warning('not all who are willing in the general public have recieved the booster dose')
      }
    }
    if (vax_doses_risk >= 2){
      if (unique(na.omit(round(hypoth_doses$cov[hypoth_doses$dose == 3 & hypoth_doses$risk_group == risk_group_name],digits=2))) != vax_strategy_toggles$vax_strategy_max_expected_cov){
        warning('not all who are willing in the risk group have recieved the booster dose')
      }
    }
  }
  
  return(vaccination_history_MODF)
}
