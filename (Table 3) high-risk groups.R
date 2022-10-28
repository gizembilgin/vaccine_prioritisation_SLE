# This program runs the 'CommandDeck' multiple times with varying prioritisation strategies involving high-risk adults.
# It compiles results from these runs in Table 3 of the main paper.


### (1) Overarching trackers #####################################################################################################
warehouse_table = data.frame() 
warehouse_plot  = data.frame()

#sensitivity analysis
if (risk_group_name == "pregnant_women"){
  RR_estimate  = RR_default =  2.4
} else if (risk_group_name == "adults_with_comorbidities"){
  RR_estimate  = RR_default = 1.95
}

#resets
vax_strategy_toggles = vax_strategy_toggles_CURRENT_TARGET
outbreak_timing = "off"
vax_strategy_toggle = "on"
vax_risk_strategy_toggle = "on"
risk_group_toggle = "on" 
risk_group_prioritisation_to_date = NA
default_prioritisation_proportion = 0.5




### (2) Queue strategies to run ##################################################################################################
queue = list()

### Baseline - no prioritisation
this_run = list(
  vax_risk_strategy = 'N',           
  vax_risk_proportion = 0,         
  vax_doses_general = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy,               
  vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy,
  risk_group_acceptability = vax_strategy_toggles$vax_strategy_max_expected_cov,
  risk_group_accessibility = FALSE,
  risk_group_age_broaden = FALSE
)
if ('vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
  this_run$risk_group_acceptability = sensitivity_analysis_toggles$vax_hesistancy_risk_group
}

queue[[1]] = list(vax_strategy_description = "no prioritisation",
                  apply_risk_strategy_toggles = this_run)


#(A/C) Prioritisation of primary doses
this_run$vax_risk_strategy = 'Y'

### single dose
this_run$vax_risk_proportion = 0.25
queue[[2]] = list(vax_strategy_description = '25% prioritisation',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak

this_run$vax_risk_proportion = 0.5
queue[[3]] = list(vax_strategy_description = '50% prioritisation',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak

this_run$vax_risk_proportion = 0.75
queue[[4]] = list(vax_strategy_description = '75% prioritisation',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak

this_run$vax_risk_proportion = default_prioritisation_proportion
###______________

### booster doses
this_run$vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy + 1
queue[[5]] = list(vax_strategy_description = 'booster at three months (at risk only)',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak
this_run$vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy
#___________________________________


#(B/C) Providing (opportunistic) additional doses
this_run$risk_group_accessibility = TRUE
### single dose
this_run$vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy
queue[[6]] = list(vax_strategy_description = 'additional primary doses',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak
###______________

### booster doses
this_run$vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy + 1
queue[[7]] = list(vax_strategy_description = 'additional booster doses',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak
#___________________________________


#(C/C) Broaden to include high-risk children
this_run$vax_doses_risk = vax_strategy_toggles_CURRENT_TARGET$vax_dose_strategy
this_run$risk_group_age_broaden = TRUE
queue[[8]] = list(vax_strategy_description = 'broaden to <18 pregnant individuals',
                  apply_risk_strategy_toggles = this_run)  #roll out vaccine DURING outbreak
#___________________________________



### (3) Run  ##################################################################################################
if ('additional_doses' %in% names(sensitivity_analysis_toggles)){
  queue = queue[c(1:4)]
}
for (ticket in 1:length(queue)){
  
  commands = queue[[ticket]]
  
  VE_loop = 0
  vax_strategy_description    = commands$vax_strategy_description
  apply_risk_strategy_toggles = commands$apply_risk_strategy_toggles
  vax_strategy_toggles = vax_strategy_toggles_CURRENT_TARGET
  
  #make additional booster doses equal priority as previous schedule
  if (length(booster_prioritised_strategies)>1){
    booster_prioritised_strategies$risk_proportion = apply_risk_strategy_toggles$vax_risk_proportion
  }
  if('additional_doses' %in% names(sensitivity_analysis_toggles)){
    if (sensitivity_analysis_toggles$additional_doses == 'booster_doses_2023'){
      apply_risk_strategy_toggles$vax_risk_strategy = "N"
    }
  }

  source(paste(getwd(),"/CommandDeck.R",sep=""))
  
  severe_outcome_projections = severe_outcome_log %>% 
    mutate(label = vax_strategy_description, day = as.numeric(date - date_start ))
  warehouse_plot = rbind(warehouse_plot,severe_outcome_projections)
  row = row %>% 
    mutate(scenario                    = vax_strategy_description,
           date_complete_at_risk_group = date_complete_at_risk_group) %>% 
    relocate(scenario, .before = colnames(row)[[1]])
  warehouse_table = rbind(warehouse_table,row)
  
  
  ### code for RR sensitivity analysis
  if ('RR_risk_group' %in% names(sensitivity_analysis_toggles)){
    RR_to_test_list = sensitivity_analysis_toggles$RR_risk_group 
    if (ticket == 1){SA_RR_warehouse_table = data.frame()} 
    
    for (RR_loop in 1:length(RR_to_test_list)){
      RR_estimate = RR_to_test_list[[RR_loop]]
      
      source(paste(getwd(),"/(5)_severe_outcomes_calc.R",sep="")) 
      source(paste(getwd(),"/(function)_severe_outcome_proj.R",sep=""))
      
      row = row %>% mutate(scenario = vax_strategy_description,
                           RR_estimate = RR_estimate) %>% 
        relocate(scenario, .before = colnames(row)[[1]])
      SA_RR_warehouse_table = rbind(SA_RR_warehouse_table,row)
    }
    save(SA_RR_warehouse_table,file =  paste(rootpath,"x_results/sensitivity_analysis_RR_",Sys.Date(),".Rdata",sep=''))
    RR_estimate = RR_default
  }
  
  ### code for VE sensitivity analysis
  if ('VE_older_adults' %in% names(sensitivity_analysis_toggles)){
    save_toggles = sensitivity_analysis_toggles
    
    if ('VE_adults_comorb' %in% names(sensitivity_analysis_toggles)){VE_loop_length = 2
    } else{VE_loop_length = 1}
      
    if (ticket == 1){SA_VE_warehouse_table = warehouse_table %>% 
      select(-date_complete_at_risk_group) %>%
      mutate(VE_mod = 'none')
    } else{
      workshop = row %>% 
        select(-date_complete_at_risk_group) %>%
        mutate(VE_mod = 'none')
      SA_VE_warehouse_table = rbind(SA_VE_warehouse_table,workshop)
    }
    
    for (VE_loop in 1:VE_loop_length){
      if ('VE_older_adults' %in% names(sensitivity_analysis_toggles) & 'VE_adults_comorb' %in% names(sensitivity_analysis_toggles)){this_sensitivity_analysis = 'VE_comorb'
      } else{this_sensitivity_analysis = 'VE_older_adults'} #COMEBACK - there MUST be a simpler way to do this in R, but a quick 30min search yielded no results
      
      source(paste(getwd(),"/(5)_severe_outcomes_calc.R",sep="")) 
      source(paste(getwd(),"/(function)_severe_outcome_proj.R",sep=""))
      
      row = row %>% mutate(scenario = vax_strategy_description,
                           VE_mod = this_sensitivity_analysis) %>% 
        relocate(scenario, .before = colnames(row)[[1]])
      SA_VE_warehouse_table = rbind(SA_VE_warehouse_table,row)
      
      sensitivity_analysis_toggles = sensitivity_analysis_toggles[!names(sensitivity_analysis_toggles) %in% c('VE_adults_comorb')] #remove for second loop
    }
    VE_loop = 0
    
    save(SA_VE_warehouse_table,file =  paste(rootpath,"x_results/sensitivity_analysis_VE_",Sys.Date(),".Rdata",sep=''))
    sensitivity_analysis_toggles = save_toggles
  }
  
  rm(severe_outcome_this_run, severe_outcome_log, severe_outcome_log_tidy)
}
#____________________________________________________________________________________________________________________________________



### (4) Save outputs  ##################################################################################################
results_warehouse_entry = list()
results_warehouse_entry[[1]] = warehouse_table
results_warehouse_entry[[2]] = warehouse_plot

### (A/B) Quick plot for sense checking
warehouse_plot = warehouse_plot %>% 
  mutate(time = day) %>%
  filter(time>=0)

if (risk_group_name == 'adults_with_comorbidities'){warehouse_plot = warehouse_plot[! warehouse_plot$label %in% c('broaden to <18 pregnant individuals'),]}

section_1 = c("no prioritisation" ,"25% prioritisation" ,"50% prioritisation" ,"75% prioritisation", "booster at three months (at risk only)" )
section_2 = c("no prioritisation","additional primary doses","additional booster doses","broaden to <18 pregnant individuals")
section_list = list(section_1,section_2)

for (section in 1:length(section_list)){
  list_plot_commands = section_list[[section]]
  workshop = warehouse_plot[warehouse_plot$label %in% list_plot_commands, ]
  
  abs_plot_list = list()
  for (i in 1:length(unique(workshop$outcome))){
    outcome = unique(workshop$outcome)[i]
    abs_plot_list [[i]] <- ggplot(data=workshop[workshop$outcome==outcome,]) + 
      geom_line(aes(x=time,y=proj,color=as.factor(label))) +
      labs(title=paste(outcome)) +
      labs(colour = "") +
      theme_bw() + 
      xlab("") + 
      ylab("")}
  
  cum_plot_list = list()
  for (i in 1:length(unique(workshop$outcome))){
    outcome = unique(workshop$outcome)[i]
    cum_plot_list [[i]] <- ggplot(data=workshop[workshop$outcome==outcome,]) + 
      geom_line(aes(x=time,y=proj_cum,color=as.factor(label))) +
      labs(title=paste(outcome)) +
      labs(colour = "") +
      theme_bw() + 
      xlab("") + 
      ylab("")
  }
  
  # 1 = death, 2 = hosp, 3 = severe_disease, 4 = YLL, 5 = cases
  ggarrange(abs_plot_list[[5]],cum_plot_list[[5]],
                   abs_plot_list[[2]],cum_plot_list[[2]],
                   abs_plot_list[[3]], cum_plot_list[[3]],
                   abs_plot_list[[1]],  cum_plot_list[[1]],
                   common.legend = TRUE,
                   legend="bottom",
                   ncol = 2,
                   nrow = 4)
}


### (B/B) Create table
baseline_to_compare = "no prioritisation"

iteration_num = 1
if ('RR_risk_group' %in% names(sensitivity_analysis_toggles)){iteration_num = length(RR_to_test_list)}
if ('VE_adults_comorb' %in% names(sensitivity_analysis_toggles)){iteration_num = 3
} else if ('VE_older_adults' %in% names(sensitivity_analysis_toggles)){iteration_num = 2 #include default
}
  
for (i in 1:iteration_num){
  
  if ('RR_risk_group' %in% names(sensitivity_analysis_toggles)){
    workshop = SA_RR_warehouse_table %>% 
      filter(RR_estimate == RR_to_test_list[[i]] ) %>% 
      select(-RR_estimate)
  } else if ('VE_older_adults' %in% names(sensitivity_analysis_toggles)){
    workshop = SA_VE_warehouse_table %>% 
      filter(VE_mod == unique(SA_VE_warehouse_table$VE_mod)[i] ) %>% 
      select(-VE_mod)
  } else {
    workshop = warehouse_table %>% 
    select(-date_complete_at_risk_group)
  }
  
  table3 = workshop %>% 
    pivot_longer(
      cols = 2:ncol(workshop) ,
      names_to = 'outcome',
      values_to = 'num'
    ) 
  
  #compare to baseline
  baseline = table3 %>% 
    filter(scenario == baseline_to_compare) %>%
    rename(baseline_num=num) %>%
    select(-scenario)
  
  table3 = table3 %>%
    left_join(baseline,by=c('outcome'))  %>%
    mutate(abs_reduction = num - baseline_num,
           rel_reduction = 100*(num - baseline_num)/baseline_num)
  
  if (risk_group_name == "pregnant_women"){
    table3$outcome = factor(table3$outcome,levels=c('cases','severe_disease','hosp','death','YLL','neonatal_deaths'))
  } else{
    table3$outcome = factor(table3$outcome,levels=c('cases','severe_disease','hosp','death','YLL'))
  }
  
  table3$scenario = factor(table3$scenario, levels = 
                             c("no prioritisation",                      "25% prioritisation",                     "50% prioritisation",                    
                               "75% prioritisation",                     "booster at three months (at risk only)", "additional primary doses",              
                               "additional booster doses",               "broaden to <18 pregnant individuals" ))
  table3 = table3 %>% arrange(scenario,outcome)
  
  options(scipen = 1000)
  print = table3 %>% 
    filter(! scenario %in% c(baseline_to_compare))  %>%
    mutate(abs_reduction = round(abs_reduction),
           rel_reduction = round(rel_reduction,digits=1),
           together_value = paste(format(abs_reduction, format="f", big.mark=",", digits=1),
                                  ' (',rel_reduction,'%)',sep=''),
           together_outcome = paste(outcome,sep='_')) %>%
    ungroup() %>%
    select(-num,-baseline_num,-abs_reduction,-rel_reduction,-outcome) %>%
    pivot_wider(
      id_cols = scenario,
      names_from = together_outcome,
      values_from = together_value)
  

  #Save results!  
  time = Sys.time()
  time = gsub(':','-',time)
  if ('RR_risk_group' %in% names(sensitivity_analysis_toggles)){
    this_RR = RR_to_test_list[[i]]
    write.csv(print,file=paste(rootpath,'x_results/table3',vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,risk_group_name,'RR',this_RR,time,'.csv',sep=''))
    
  } else if ('VE_older_adults' %in% names(sensitivity_analysis_toggles)){
    this_VE_mod = unique(SA_VE_warehouse_table$VE_mod)[i] 
    write.csv(print,file=paste(rootpath,'x_results/table3',vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,this_VE_mod,time,'.csv',sep=''))
    
  } else if ('vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
    write.csv(print,file=paste(rootpath,'x_results/table3',vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,risk_group_name,' vax_hest ',time,'.csv',sep=''))
    results_warehouse_entry[[4]]= print
    
  } else if('additional_doses' %in% names(sensitivity_analysis_toggles)){
    write.csv(print,file=paste(rootpath,'x_results/table3_SA_additional_doses',sensitivity_analysis_toggles$additional_doses,risk_group_name,time,'.csv',sep=''))
    
  } else{
    write.csv(print,file=paste(rootpath,'x_results/table3',vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,risk_group_name,gov_target,time,'.csv',sep=''))
    results_warehouse_entry[[4]]= print
    results_warehouse[[receipt]] = results_warehouse_entry
  }
}



