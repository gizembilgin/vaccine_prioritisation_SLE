### This program conducts sensitivity analysis for lower and higher levels of seroprevalence
### in the inital state of the model. It's aim to gauge the effect of the model's fit
### on the paper's results.

### DEPENDENCIES: nil!
rm(list=ls())


### Setup ____________________________________________________________________________________________________________
#start timing
time.start.FleetAdmiral=proc.time()[[3]]

#load latest fit
load(file = '1_inputs/last_fit_date.Rdata')
date_start = fitted_max_date ##latest fit date

#initialise length of model run and circulating strain
strain_inital = strain_now = 'omicron' 
outbreak_timing = "off" #roll-out during steady state
model_weeks = 52
age_split_results = "N"

#turn on waning of all immunity
waning_toggle_acqusition = TRUE
waning_toggle_severe_outcome = TRUE
waning_toggle_rho_acqusition = TRUE

#turn off risk groups to start with
risk_group_toggle = "off"
vax_risk_strategy_toggle = "off"
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

#initialise data frame
results_warehouse = list()
sensitivity_analysis_toggles = list()
queue = list(
  list(modification_factor_on_preexisting_immunity = 0.5),
  list(modification_factor_on_preexisting_immunity = 1.5)
)
receipt = 0
#_______________________________________________________________________________

for (place_in_queue in 1:length(queue)){
  
  sensitivity_analysis_toggles = queue[[place_in_queue]]
  
  ###(Table 2) prioritisation strategies including children 5 to 17
  receipt = receipt + 1
  source(paste(getwd(),"/(Table 2)_varying_eligb_age.R",sep=""))
  #_____________________________________________________________________________
  
  
  
  ###(Table 3) prioritisation of high-risk groups
  receipt = receipt + 1
  risk_group_name = "pregnant_women"
  source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))
  
  receipt = receipt + 1
  risk_group_name = "adults_with_comorbidities"
  source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 
  
  risk_group_toggle = "off"
  vax_risk_strategy_toggle = "off"
  #_____________________________________________________________________________
}




save.image(file = paste(rootpath,"x_results/sensitivity_analysis_effect_of_fit_",Sys.Date(),".Rdata",sep=''))

time.end.FleetAdmiral=proc.time()[[3]]
time.end.FleetAdmiral-time.start.FleetAdmiral 