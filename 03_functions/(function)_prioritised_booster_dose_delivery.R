### This function allocates future doses by risk_group using repeated calls of the booster_strategy() function

booster_strategy_prioritised <- function(
    booster_risk_strategy,             # options: 'Y','N'
    booster_risk_proportion,           # value between 0-1 (equivalent to %) of doses prioritised to the at risk group
    risk_group_accessibility = FALSE, #TRUE means opportunistic vaccination possible for this risk group,
    booster_toggles_import = booster_toggles
){
  
  ### WARNINGS 
  if (!booster_risk_strategy %in% c('Y','N')){stop('Is the vax strategy on or off? Y/N')}
  if (booster_risk_proportion<0 | booster_risk_proportion>1){stop('booster_risk_proportion must be between 0-1 (0%-100%)')}
  if (booster_risk_proportion == 0 & booster_risk_strategy == 'Y'){
    warning("Giving 0% priority is not priority! I have overwritten booster_risk_strategy to EQ 'N'")
    booster_risk_strategy = "N"
  }
  if (length(booster_toggles_import$delivery_risk_group)==1){
    stop("Can't prioritise within a single risk group!")
  }
  #_______________________________________________________________________________
  
  
  
  ### BRANCH ONE: Are we prioritising the at risk group at all?
  if (booster_risk_strategy == "N"){
    #if not prioritised, set booster_risk_proportion to population proportion of risk group
    risk_split = age_risk_split %>% group_by(risk_group) %>% summarise(sum = sum(split),.groups='keep')
    booster_risk_proportion = risk_split$sum[risk_split$risk_group == risk_group_name]
  }
  #___________________________________________________________
  
  
  
  ### BRANCH TWO: What % priority are the at risk group receiving? 
  if (risk_group_accessibility == FALSE){
    speed_risk_group_rollout =  booster_toggles_import$rollout_speed*booster_risk_proportion
  } else if (risk_group_accessibility == TRUE){
    speed_risk_group_rollout = round(sum(pop_high_risk$pop)/(365/12*4.2)) #median first antenatal care visit
  }

  at_risk_delivery_outline = booster_strategy(booster_strategy_start_date = booster_toggles_import$start_date,       # start of hypothetical vaccination program
                                              booster_dose_allocation     = booster_toggles_import$dose_allocation,  # num of doses avaliable
                                              
                                              booster_rollout_speed       = speed_risk_group_rollout,         # doses delivered per day
                                              booster_delivery_risk_group = risk_group_name,
                                              
                                              booster_delivery_includes_previously_boosted = booster_toggles_import$delivery_includes_previously_boosted,
                                              booster_age_strategy        = booster_toggles_import$age_strategy,     # options: "oldest", "youngest","50_down","uniform"
                                              booster_strategy_vaccine_type = booster_toggles_import$vaccine_type,   # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"  
                                              booster_strategy_vaccine_interval = booster_toggles_import$vaccine_interval,
                                              
                                              booster_prioritised = 'Y'
  )
  at_risk_delivery_outline = at_risk_delivery_outline %>% mutate(risk_group = risk_group_name)
  
  ggplot(at_risk_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  sum(at_risk_delivery_outline$doses_delivered_this_date)
  total = at_risk_delivery_outline %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep")
  ggplot(total) + geom_point(aes(x=date,y=doses_delivered_this_date))
  
  vaccination_history_FINAL = rbind(vaccination_history_FINAL,at_risk_delivery_outline)
  #___________________________________________________________
  
  
  
  ### BRANCH THREE: Distribute remaining doses to general population
  generalPublic_leftover_outline = booster_strategy(booster_strategy_start_date = booster_toggles_import$start_date,       # start of hypothetical vaccination program
                   booster_dose_allocation     = booster_toggles_import$dose_allocation - sum(at_risk_delivery_outline$doses_delivered_this_date),  # num of doses avaliable
                   booster_rollout_speed       = booster_toggles_import$rollout_speed,         # doses delivered per day
                   booster_delivery_risk_group = "general_public",
                   booster_delivery_includes_previously_boosted = booster_toggles_import$delivery_includes_previously_boosted,
                   booster_age_strategy        = booster_toggles_import$age_strategy,     # options: "oldest", "youngest","50_down","uniform"
                   booster_strategy_vaccine_type = booster_toggles_import$vaccine_type,   # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"  
                   booster_strategy_vaccine_interval = booster_toggles_import$vaccine_interval,
                   booster_prioritised = 'Y',
                   vaccination_history_FINAL_local = vaccination_history_FINAL
  )
  generalPublic_leftover_outline = generalPublic_leftover_outline %>% mutate(risk_group = "general_public")
  
  
  ggplot(generalPublic_leftover_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  sum(generalPublic_leftover_outline$doses_delivered_this_date)
  total = generalPublic_leftover_outline %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date),.groups = "keep")
  ggplot(total) + geom_point(aes(x=date,y=doses_delivered_this_date))  
  
  vaccination_history_FINAL = rbind(vaccination_history_FINAL,generalPublic_leftover_outline)
  
  ###CHECKS
  #CHECK 1: total doses delivered <= total doses available
  if (risk_group_accessibility == FALSE){
    check_df = bind_rows(at_risk_delivery_outline,generalPublic_leftover_outline)
    if (round(sum(check_df$doses_delivered_this_date)) > booster_toggles_import$dose_allocation){stop('Total doses delivered > total doses avaliable!')}

    check_df_daily = check_df %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label ='total')
    if (nrow(check_df_daily[round(check_df_daily$doses_delivered_this_date)>booster_toggles_import$rollout_speed,])>0){
      ggplot(check_df_daily) + geom_point(aes(x=date,y=doses_delivered_this_date))

      df1 = at_risk_delivery_outline  %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label = 'at risk delivery')
      df2 = generalPublic_leftover_outline  %>% group_by(date) %>% summarise(doses_delivered_this_date = sum(doses_delivered_this_date)) %>% mutate(label = 'general public leftover')
      check_df_daily_comp = bind_rows(df1,df2,check_df_daily)
      ggplot(check_df_daily_comp) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(label)))
      ggplot(check_df_daily_comp[check_df_daily_comp$label != 'total',]) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(label)))

      nrow_greater = nrow(check_df_daily[round(check_df_daily$doses_delivered_this_date)>booster_toggles_import$rollout_speed,])

      warning('Caution, more doses delivered per day than capacity for',paste(nrow_greater),' days')

      if (nrow_greater>10){
        stop('Error maximum daily doses exceeded for more than 10 days')
      }
      if ((max(check_df_daily$doses_delivered_this_date)-booster_toggles_import$rollout_speed )/booster_toggles_import$rollout_speed > 0.1){
        stop('Error maximum daily doses exceeded by 10%') #NOTE - 10% leeway, unless want more complicated limiter on booster dose
      }

    }
    ggplot(at_risk_delivery_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
    ggplot(generalPublic_leftover_outline) + geom_point(aes(x=date,y=doses_delivered_this_date,color=as.factor(age_group),shape=as.factor(dose)))
  }
  #CHECK: for non-prioritised dose delivery with risk groups
  if (booster_risk_strategy == "N"){
    hypoth_doses = vaccination_history_FINAL %>% 
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
    
    if (max(vaccination_history_FINAL$date[vaccination_history_FINAL$risk_group == 'general_public']) !=
        max(vaccination_history_FINAL$date[vaccination_history_FINAL$risk_group == risk_group_name])) {
        warning('max delivery dates dont align between risk groups for additional booster dose')
    }
  }
  
  #CHECK:gaps in dates in vaccination_history_FINAL
  filler_dates = seq(min(vaccination_history_FINAL$date),max(vaccination_history_FINAL$date), by = 'days') #expected dates!
  filler_dates = filler_dates[! filler_dates %in% vaccination_history_FINAL$date]
  copy = vaccination_history_FINAL %>% ungroup %>% filter(date == max(vaccination_history_FINAL$date))
  #COMEBACK - this is computationally inefficient
  for (i in 1:length(filler_dates)){
    copy = copy %>% mutate(date = filler_dates[i], doses_delivered_this_date = 0)
    vaccination_history_FINAL = rbind(vaccination_history_FINAL,copy)
  }
  
  return(vaccination_history_FINAL)
}
