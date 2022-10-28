### The 'CommandDeck' runs all sub-scripts of the COVID-19 transmission model to
### complete one standard 'run' of the disease model.

rm(list=ls())
options(warn = 0)

#### SETUP ##################################################################################
#load libraries
library(readr)
library(deSolve)
library(rvest)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyverse)

load(file = '1_inputs/last_fit_date.Rdata')
rootpath = str_replace(getwd(), "GitHub_vaxAllocation","") #Note: x_results not stored within GitHub repository

# Load functions
source(paste(getwd(),"/(function)_COVID_ODE.R",sep=""))
source(paste(getwd(),"/(function)_VE_time_step.R",sep=""))
source(paste(getwd(),"/(function)_rho_time_step.R",sep=""))
source(paste(getwd(),"/(function)_vax_strategies.R",sep=""))
source(paste(getwd(),"/(function)_vax_strategies_risk.R",sep=""))
if (exists("VE_estimates_imputed") == FALSE){load(file='1_inputs/VE_estimates_imputed.Rdata')}
#_____________________________________________________________________________________________


#### FUNCTION BODY ###########################################################################
CommandDeck <- function(
    
  ticket = 1,
  
  date_start = fitted_max_date,  
  model_weeks = 5,
  
  outbreak_timing = "off",
  strain_inital = 'omicron',     
    
  vax_strategy_toggle = "off",
  vax_risk_strategy_toggle = "off",
  risk_group_toggle = "on", 
  risk_group_name = "adults_with_comorbidities", #options: pregnant_women, adults_with_comorbidities
  RR_estimate  = 2,
  risk_group_prioritisation_to_date = NA,
  default_prioritisation_proportion = 0.5,
  risk_group_lower_cov_ratio = NA,
  sensitivity_analysis_toggles = list(),
  
  vax_strategy_toggles =
    list(vax_strategy_start_date        = date_start+30,
         vax_strategy_num_doses         = as.integer(1642011),
         vax_strategy_roll_out_speed    = 11075 ,                           # doses delivered per day
         vax_delivery_group             = 'universal',
         vax_age_strategy               = "uniform_no_children",            # options: "oldest", "youngest","50_down","uniform"
         vax_dose_strategy              = 1,                                # options: 1,2
         vax_strategy_vaccine_type      = "Johnson & Johnson" ,             # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"
         vax_strategy_vaccine_interval  = c(30*3) ,                         # (days) interval between doses, you must specify multiple intervals if multiple doses e.g. c(21,90)
         vax_strategy_max_expected_cov  = 0.88                              # value between 0-1 of age group willing to be vaccinated
    ),
  
  apply_risk_strategy_toggles = list(
    vax_risk_strategy = 'Y',             # options: 'Y','N'
    vax_risk_proportion = 0.8,           # value between 0-1 (equivalent to %) of doses prioritised to the at risk group
    vax_doses_general = 1,               # number of doses delivered to general pop
    vax_doses_risk = 2                   # number of doses delivered to risk group
  ),
  
  waning_toggle_acqusition = TRUE,
  waning_toggle_severe_outcome = TRUE, 
  waning_toggle_rho_acqusition = TRUE,
    
  setting = "SLE",
  
  age_split_results = "N",
  fitting = "off"
 
){

  #       (1/3) Initalise model state                 
  ####################################################################
  strain_now = 'omicron'   
  complete_model_runs = 1 # when >1 samples randomly from distribution of parameters (where available)
  
  if (risk_group_toggle == "on"){
    num_risk_groups = 2
  } else{ num_risk_groups = 1; vax_risk_strategy_toggle = "off"}

  # load refit, or refit if not recent enough
  if (fitting == "on"){
    warning('Fitting is on')
  } else if ( ! 'vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
    if (as.numeric(abs(fitted_max_date - Sys.Date()))>30){ 
      warning('refitting model as fitted_max_date over one month since today!')
      source(paste(getwd(),"/(0)_fitting_model.R",sep=""))
    } else{
      load(file = '1_inputs/fitted_results.Rdata')
      
      if('additional_doses' %in% names(sensitivity_analysis_toggles)){
        if (sensitivity_analysis_toggles$additional_doses == 'start_2022'){
          load(file = '1_inputs/fitted_results_SA_2022.Rdata')
        }
      }
      
      if (risk_group_toggle == "off"){
        loaded_fit = fitted_results[[1]]
      } else if (risk_group_name == 'pregnant_women'){
        loaded_fit = fitted_results[[2]]
      } else if (risk_group_name == 'adults_with_comorbidities'){
        loaded_fit = fitted_results[[3]]
      }
      rm(fitted_results)
      
      if (risk_group_toggle == "on"){if(!loaded_fit[[5]] == risk_group_name){stop('risk group name != fitted risk group name')}}
      
      parameters = loaded_fit[[1]]
      fitted_next_state = loaded_fit[[2]]
      fitted_incidence_log_tidy = loaded_fit[[3]]
      fitted_incidence_log = loaded_fit[[4]]
      rm(loaded_fit)
      
      fitted_incidence_log_tidy = fitted_incidence_log_tidy %>% filter(date <= date_start) 
      fitted_incidence_log = fitted_incidence_log %>% filter(date <= date_start)
      
      if (risk_group_toggle == "on"){
        if ((is.na(risk_group_prioritisation_to_date) == FALSE) ){
          stop('no fitted result avaliable for this risk group characteristic')
        }
      }
    }
  } else if('vax_hesistancy_risk_group' %in% names(sensitivity_analysis_toggles)){
    
    if (! risk_group_name == 'pregnant_women'){stop('havent configured vax hesistance sensitivity analysis for another risk group')}
    
    load(file = '1_inputs/SA_vaxHest_fitted_results.Rdata')
    loaded_fit = SA_vaxHest_fitted_results
    
    parameters = loaded_fit[[1]]
    fitted_next_state = loaded_fit[[2]]
    fitted_incidence_log_tidy = loaded_fit[[3]]
    fitted_incidence_log = loaded_fit[[4]]
    rm(loaded_fit)
    
    fitted_incidence_log_tidy = fitted_incidence_log_tidy %>% filter(date <= date_start) # CHECKED last of fitted log = first of new log
    fitted_incidence_log = fitted_incidence_log %>% filter(date <= date_start)
  } 
  
  
  if ( fitting == "on"){
    Reff_tracker = data.frame()
    rho_tracker_dataframe = data.frame()
    VE_tracker_dataframe = data.frame()
  }
  #__________________________________________________________________
  
  

  #       (2/3) Run model            
  #####################################################################
  source(paste(getwd(), "/(1)_simulate_setting.R", sep = "")) #load setting stats if new setting
  
  #making some interim variables to assist with configuring states
  num_disease_classes = 4
  disease_class_list = c('S', 'E', 'I', 'R')
  num_vax_doses = D = length(unique(vaccination_history_TRUE$dose))
  vax_type_list = sort(unique(vaccination_history_TRUE$vaccine_type))
  num_vax_types = T = length(unique(vaccination_history_TRUE$vaccine_type))
  num_vax_classes = num_vax_doses * num_vax_types + 1 # + 1 for unvaccinated
  
  source(paste(getwd(), "/(3)_disease_characteristics.R", sep = ""))
  source(paste(getwd(), "/(2)_inital_state.R", sep = ""))
  source(paste(getwd(), "/(4)_time_step.R", sep = ""))
  source(paste(getwd(), "/(5)_severe_outcomes_calc.R", sep = "")) # COMEBACK - should this just save its results somewhere?
  source(paste(getwd(), "/(function)_severe_outcome_proj.R", sep = ""))
  #__________________________________________________________________
  

    
  #       (3/3) Basic plots            
  #####################################################################
  plot_standard = theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'))

  if (fitting == "on"){
    plot1 <- 
      ggplot() + 
      geom_line(data=incidence_log_plot,aes(x=date,y=rolling_average_percentage),na.rm=TRUE) +
      xlab("") +
      scale_x_date(date_breaks="1 month", date_labels="%b") +
      ylab("daily cases % whole pop") +
      plot_standard
    
    plot2 <- ggplot() + 
      geom_line(data=Reff_tracker,aes(x=date,y=Reff),na.rm=TRUE) +
      xlab("") + 
      scale_x_date(date_breaks="1 month", date_labels="%b") +
      ylab("Reff") +
      plot_standard
    
    plot3<- ggplot() + 
      geom_line(data=incidence_log_plot,aes(x=date,y=cumulative_incidence_percentage),na.rm=TRUE) +
      xlab("") + 
      scale_x_date(date_breaks="1 month", date_labels="%b") +
      ylab("cumulative cases % whole pop") +
      plot_standard

    plot4 = ggplot(rho_tracker_dataframe) + geom_line(aes(x=date,y=rho))
    
    plot5 = ggplot(VE_tracker_dataframe) + geom_line(aes(x=date,y=VE,color=as.factor(dose)))
  } 
  #__________________________________________________________________ 

  
  #copy list of variables to the global environment
  vars = list(severe_outcome_log = severe_outcome_log,
              row = row)
  
  #additional variables used in (5) and (6) only for S.A.
  if ('RR_risk_group' %in% names(sensitivity_analysis_toggles) | 'VE_older_adults' %in% names(sensitivity_analysis_toggles)){
    vars_extra <- list(
      save_VE_waning_distribution = save_VE_waning_distribution,
      strain_now = strain_now,
      date_start = date_start,
      risk_group_labels = risk_group_labels,
      D = D,
      vax_type_list = vax_type_list,
      age_group_labels = age_group_labels,
      severe_outcome_this_run = severe_outcome_this_run,
      severe_outcome_log_tidy = severe_outcome_log_tidy,
      seed_date = seed_date,
      exposed_log = exposed_log)
    vars = append(vars,vars_extra)
    #COMEBACK - include a rm() in Table 2 or 3 for these variables when implemented
  }
  
  for (i in 1:length(vars)){
    assign(names(vars)[[i]], vars[[i]], envir = .GlobalEnv)
  }
  
  return()
}

CommandDeck()