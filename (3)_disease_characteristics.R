### This program sets the disease parameters to the currently circulating strain of COVID-19

strain_now = strain_inital

### (A/D) Transmission
load(file = "1_inputs/param_age.Rdata")
suscept = param_age$value[param_age$param == 'susceptibility']    # (i) age-specific susceptibility to infection
gamma   = param_age$value[param_age$param == 'prop_sympt']        # (ii) proportion of cases symptomatic
lota    = 0.5                                                     # (iii) modification factor on infectiousness of asymptomatic cases


### (B/D) Latent period 
if (run_type == "point"){  
  if (strain_inital == 'delta' | strain_inital == 'WT'){
    AverageLatentPeriod = 3.71
  } else if (strain_inital == 'omicron'){
    AverageLatentPeriod = 2.22
  }
}else if (run_type == "rand"){AverageLatentPeriod = rlnorm(1,meanlog = 1.3, sd=0.2)} 

lambda = 1/AverageLatentPeriod


### (C/D) Symptomatic period
if (run_type == "point"){
  if (strain_inital == 'delta' | strain_inital == 'WT'){
    AverageSymptomaticPeriod = 10.9
  } else if (strain_inital == 'omicron'){
    AverageSymptomaticPeriod = 9.87
  }
}else if (run_type == "rand"){AverageSymptomaticPeriod = runif(1,min=7,max=14)} #taken from Zachreson et al., 2021, COMEBACK use ACT bounds

delta = 1/AverageSymptomaticPeriod


### (D/D) Waning of infection-derived immunity
lengthInfectionDerivedImmunity = 180 #days
omega = 1/lengthInfectionDerivedImmunity

