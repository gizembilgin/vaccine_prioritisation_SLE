### This program will generate sensitivity analysis for the prioritisation of high-risk groups with additional doses under two cirumstances:
# (1) running the model from start of 2022
# (2) a booster dose eligibilty to all adults is introducted


#### Running the model from start of 2022 #####################################################################
rm(list=ls())

sensitivity_analysis_toggles = list(additional_doses = 'start_2022')

date_start = as.Date('2022-01-01')
strain_inital = strain_now = 'omicron' 
outbreak_timing = "off" #roll-out during steady state
model_weeks = 52

#turn on waning of all immunity
waning_toggle_acqusition = TRUE
waning_toggle_severe_outcome = TRUE
waning_toggle_rho_acqusition = TRUE

risk_group_toggle = 'on'
risk_group_lower_cov_ratio = NA
risk_group_prioritisation_to_date = NA

#set up setting
setting = "SLE"
source(paste(getwd(),"/(1)_simulate_setting.R",sep=""))
if (setting == "SLE"){
  gov_target = 0.516
  workshop_doses = gov_target - sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose == 1])/100
  workshop_doses = round(workshop_doses * sum(pop))
  
  vax_strategy_toggles_CURRENT_TARGET =
    list(vax_strategy_start_date        = date_start,
         vax_strategy_num_doses         = as.integer(workshop_doses),
         vax_strategy_roll_out_speed    = 11075 ,                           # doses delivered per day
         vax_delivery_group             = 'universal',
         vax_age_strategy               = "uniform_no_children",            # options: "oldest", "youngest","50_down","uniform"
         vax_dose_strategy              = 1,                                # options: 1,2
         vax_strategy_vaccine_type      = "Johnson & Johnson" ,             # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"
         vax_strategy_vaccine_interval  = c(90) ,                           #  (days) interval between doses, you must specify multiple intervals if multiple doses e.g. c(21,90)
         vax_strategy_max_expected_cov  = 0.88                              # value between 0-1 of age group willing to be vaccinated
    )
} else { stop ('pick a valid setting!')}
#________________________________________________________________________________________________________________


###(Table 3) prioritisation of high-risk groups
receipt = 2
risk_group_name = "pregnant_women"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))
#________________________________________________________________________________________________________________

receipt = 3
risk_group_name = "adults_with_comorbidities"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 
#________________________________________________________________________________________________________________
#################################################################################################################



#### Expand to booster dose to everyone in 2023###################################################################
rm(list=ls())

sensitivity_analysis_toggles = list(additional_doses = 'booster_doses_2023')  #to skip queue 5 onwards in Table 3 call

load(file = '1_inputs/last_fit_date.Rdata')
date_start = fitted_max_date 

strain_inital = strain_now = 'omicron' 
outbreak_timing = "off" #roll-out during steady state
model_weeks = round(as.numeric(as.Date('2023-01-01') - date_start)/7)+52

#turn on waning of all immunity
waning_toggle_acqusition = TRUE
waning_toggle_severe_outcome = TRUE
waning_toggle_rho_acqusition = TRUE

risk_group_toggle = 'on'
risk_group_lower_cov_ratio = NA
risk_group_prioritisation_to_date = NA

#set up setting
setting = "SLE"
source(paste(getwd(),"/(1)_simulate_setting.R",sep=""))
if (setting == "SLE"){
  gov_target = 0.516
  workshop_doses = gov_target - sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose == 1])/100
  workshop_doses = round(workshop_doses * sum(pop))
  
  vax_strategy_toggles_CURRENT_TARGET =
    list(vax_strategy_start_date        = date_start,
         vax_strategy_num_doses         = as.integer(workshop_doses),
         vax_strategy_roll_out_speed    = 11075 ,                           # doses delivered per day
         vax_delivery_group             = 'universal',
         vax_age_strategy               = "uniform_no_children",            # options: "oldest", "youngest","50_down","uniform"
         vax_dose_strategy              = 1,                                # options: 1,2
         vax_strategy_vaccine_type      = "Johnson & Johnson" ,             # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"
         vax_strategy_vaccine_interval  = c(90) ,                           #  (days) interval between doses, you must specify multiple intervals if multiple doses e.g. c(21,90)
         vax_strategy_max_expected_cov  = 0.88                              # value between 0-1 of age group willing to be vaccinated
    )
} else { stop ('pick a valid setting!')}


#________________________________________________________________________________________________________________


###(Table 3) prioritisation of high-risk groups
receipt = 2
risk_group_name = "pregnant_women"

booster_toggles = list(start_date = as.Date('2023-01-01'),
                       dose_allocation = 9999999999,
                       rollout_speed = vax_strategy_toggles_CURRENT_TARGET$vax_strategy_roll_out_speed,
                       delivery_risk_group = c('general_public',risk_group_name),
                       delivery_includes_previously_boosted = 'N',
                       age_strategy = vax_strategy_toggles_CURRENT_TARGET$vax_age_strategy,
                       vaccine_type = vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,
                       vaccine_interval = 90)
booster_prioritised_strategies = list(strategy = 'Y',
                                      risk_proportion = 99)

source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))
#________________________________________________________________________________________________________________

receipt = 3
risk_group_name = "adults_with_comorbidities"

booster_toggles = list(start_date = as.Date('2023-01-01'),
                       dose_allocation = 9999999999,
                       rollout_speed = vax_strategy_toggles_CURRENT_TARGET$vax_strategy_roll_out_speed,
                       delivery_risk_group = c('general_public',risk_group_name),
                       delivery_includes_previously_boosted = 'N',
                       age_strategy = vax_strategy_toggles_CURRENT_TARGET$vax_age_strategy,
                       vaccine_type = vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,
                       vaccine_interval = 90)

source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 
#________________________________________________________________________________________________________________
#################################################################################################################






