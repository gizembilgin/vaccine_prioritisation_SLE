### This program is an alternative version of 'FleetAdmiral' whereby vaccines donated to our study setting are
### of the type Pfizer (double-dose) instead of Johnson & Johnson (single-dose). 
### This program re-runs the main results of our paper with Pfizer.

### DEPENDENCIES: nil!
rm(list=ls())


### Setup ____________________________________________________________________________________________________________
#start timing
time.start.FleetAdmiral=proc.time()[[3]]

#load latest fit
load(file = '01_inputs/last_fit_date.Rdata')
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
  workshop_doses = round(workshop_doses * sum(pop)) * 2 #ATTENTION: times two for double-dose vaccine!!!
  
  vax_strategy_toggles_CURRENT_TARGET =
    list(vax_strategy_start_date        = date_start,
         vax_strategy_num_doses         = as.integer(workshop_doses),
         vax_strategy_roll_out_speed    = 11075 ,               # doses delivered per day
         vax_delivery_group             = 'universal',
         vax_age_strategy               = "uniform_no_children",            # options: "oldest", "youngest","50_down","uniform"
         vax_dose_strategy              = 2,                                # options: 1,2
         vax_strategy_vaccine_type      = "Pfizer" ,                        # options: "Moderna","Pfizer","AstraZeneca","Johnson & Johnson","Sinopharm","Sinovac"
         vax_strategy_vaccine_interval  = c(21,90) ,                        # (days) intervals between doses
         vax_strategy_max_expected_cov  = 0.88                              # value between 0-1 of age group willing to be vaccinated (vaccine hesitancy est in discussion)
    )
} else { stop ('pick a valid setting!')}

#initialise data frame
results_warehouse = list()
sensitivity_analysis_toggles = list()
#________________________________________________________________________________________________________________



###(Table 2) prioritisation strategies including children 5 to 17
receipt = 1
source(paste(getwd(),"/(Table 2)_varying_eligb_age.R",sep=""))
#________________________________________________________________________________________________________________



###(Table 3) prioritisation of high-risk groups
receipt = 2
risk_group_name = "pregnant_women"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))

receipt = 3
risk_group_name = "adults_with_comorbidities"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 

risk_group_toggle = "off"
vax_risk_strategy_toggle = "off"
#________________________________________________________________________________________________________________


current_coverage = c(sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose ==1]),
                     sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose ==2])) #COMEBACK - if J&J in use!
if ("Johnson & Johnson" %in% unique(vaccination_history_POP$vaccine_type)){warning('True vaccine coverage MUST consider J&J dose 1')}

save.image(file = paste(rootpath,"x_results/sensitivity_analysis_Pfizer_",Sys.Date(),".Rdata",sep=''))

time.end.FleetAdmiral=proc.time()[[3]]
time.end.FleetAdmiral-time.start.FleetAdmiral 