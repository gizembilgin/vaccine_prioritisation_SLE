### The 'FleetAdmiral' is the highest overarching level of this R Project.
### This program runs all (Table *) and (sensitivity)* scripts required to generate the results
### of the vaccine allocation paper and Supplementary Material.

### DEPENDENCIES: nil!
rm(list=ls())


### Setup ____________________________________________________________________________________________________________
#start timing
time.start.FleetAdmiral=proc.time()[[3]]

#load latest fit
load(file = '1_inputs/last_fit_date.Rdata')
date_start = fitted_max_date 

#initialise length of model run and circulating strain
strain_inital = strain_now = 'omicron' 
outbreak_timing = "off" #roll-out during steady state
model_weeks = 52

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



###(Figure S4.2) outbreak of a new immune-escape variant
receipt = 4
source(paste(getwd(),"/(sensitivity)_new_variant_timing_of_outbreak.R",sep=""))
receipt = 5
source(paste(getwd(),"/(sensitivity)_new_variant_rollout_speed.R",sep=""))

results_warehouse[[4]][[3]]
results_warehouse[[5]][[3]]
#________________________________________________________________________________________________________________



###(sensitivity analysis) varying level of risk in high-risk groups
receipt = 6
sensitivity_analysis_toggles = list(RR_risk_group = list(1,1.5,2.4,3))
risk_group_name = "pregnant_women"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))

receipt = 7
sensitivity_analysis_toggles = list(RR_risk_group = list(1,1.5,1.95,3))
risk_group_name = "adults_with_comorbidities"
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 
#________________________________________________________________________________________________________________



###(sensitivity analysis) reduced VE in older adults and/or adults with comorbidities
#(Table S4.4) allocating vaccines to children with reduced vaccine effectiveness in older adults
receipt = 8
sensitivity_analysis_toggles = list(VE_older_adults = "reduced")
source(paste(getwd(),"/(Table 2)_varying_eligb_age.R",sep=""))

#(Table S4.5) allocating vaccines to pregnant women with reduced vaccine effectiveness in older adults
receipt = 9
risk_group_name = "pregnant_women"
sensitivity_analysis_toggles = list(VE_older_adults = "reduced")
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep=""))

#(Table S4.5 & S4.6) allocating vaccines to high-risk adults with reduced vaccine effectiveness in older adults and/or adults with comorbidities
receipt = 10
risk_group_name = "adults_with_comorbidities"
sensitivity_analysis_toggles = list(VE_older_adults = "reduced",VE_adults_comorb = 0.9)
source(paste(getwd(),"/(Table 3) high-risk groups.R",sep="")) 

sensitivity_analysis_toggles = list()
#________________________________________________________________________________________________________________



current_coverage = c(sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose ==1]),
                     sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose ==2])) #COMEBACK - if J&J in use!
if ("Johnson & Johnson" %in% unique(vaccination_history_POP$vaccine_type)){warning('True vaccine coverage MUST consider J&J dose 1')}

save.image(file = paste(rootpath,"x_results/complete_model_run_",Sys.Date(),".Rdata",sep=''))

time.end.FleetAdmiral=proc.time()[[3]]
time.end.FleetAdmiral-time.start.FleetAdmiral 


# time = Sys.time()
# time = gsub(':','-',time)
# file_name = paste(rootpath,"x_results/Vaccine allocation project results",time)
# file_name = gsub(' ','_',file_name)
#
# library(rmarkdown); library(tinytex)
# render('FleetAdmiral_compiler.Rmd',output_file = file_name)
# render('FleetAdmiral_compiler.Rmd',output_file = file_name, output_format = "pdf_document")
# render('FleetAdmiral_compiler.Rmd',output_file = file_name, output_format = "word_document")
