### This program runs the model over each time_step, updating vaccination rates daily

options(warn = 0) #options = 0 to turn off, 2 to stop at first warning

sol_log = data.frame()
sol_log_unedited = data.frame()
incidence_log = data.frame()

time_step = 1
num_time_steps = model_weeks *7
if (outbreak_timing == "after"){
  num_time_steps = model_weeks *7 + as.numeric(seed_date-date_start) -7
}
  

for (increments_number in 1:num_time_steps){

  if (fitting == "on" & increments_number == 1){
    
    parameters = c(
      suscept = suscept,
      beta=beta,
      NPI=NPI_inital,
      contact_matrix=contact_matrix,
      lota=lota,
      gamma=gamma,
      lambda=lambda,
      delta=delta,
      omega=omega,
      rho=rho_inital,
      age_group_labels=age_group_labels,
      risk_group_labels = risk_group_labels,
      VE=VE_inital,
      # VE_onwards=VE_onwards,
      num_age_groups=num_age_groups,
      num_risk_groups = num_risk_groups,
      num_disease_classes = num_disease_classes,
      num_vax_types=num_vax_types,
      num_vax_doses=num_vax_doses)
    rm(rho_inital,NPI_inital,VE_inital)
    
    sol = as.data.frame(ode(y=state,times=(seq(0,time_step,by=1)),func=covidODE,parms=parameters))
    
    sol_log <- sol
    sol_log_unedited <- sol
   
    Reff <- NA
    Reff_tracker = rbind(Reff_tracker,Reff)
    colnames(Reff_tracker) <- c('Reff')
    
  } else{
    
    date_now = date_start + increments_number*time_step
    
    if (increments_number > 1){ 
    
      if (date_now <= max(NPI_estimates$date)){
        NPI_this_step <- NPI_estimates$NPI[NPI_estimates$date == date_now]/100
        parameters$NPI = NPI_this_step
      } #i.e. assume after end date that NPI constant
        
      if ((date_now - min(vaxCovDelay$delay))>= min(vaccination_history_FINAL$date)){
        parameters$VE = VE_time_step(strain_now,date_now,'any_infection')
      }
      
      if (waning_toggle_rho_acqusition == TRUE ){
        parameters$rho = rho_time_step(date_now)
        rho = parameters$rho
      }
      
        state_working=tail.matrix(sol,1)
        state_working=select(state_working,-time) #remove column with time
        state_working=as.vector(state_working)
        
        # lets reconstruct our tidy matrix (easier to work with)
        A=RISK*J*(T*D+1) # +1 is unvax
        
        S = as.matrix(state_working[1:A])
        E = as.matrix(state_working[(A+1):(2*A)])
        I = as.matrix(state_working[(2*A+1):(3*A)])
        R = as.matrix(state_working[(3*A+1):(4*A)])
        Incid = as.matrix(state_working[(4*A+1):(5*A)]) 
        
        prev_state = data.frame()
        class_list = list(S,E,I,R)
        class_name_list = c('S','E','I','R')

        for (i in 1:num_disease_classes){
          workshop = data.frame(pop = class_list[[i]])
          row.names(workshop) = NULL
          workshop = workshop %>% mutate(class =  class_name_list[i])
          workshop$temp = rep(seq(1,(num_age_groups*num_vax_classes)),RISK)
          workshop$age_group = rep(age_group_labels,num_vax_classes*RISK)
          workshop$dose = 0
          workshop$vaccine_type = "unvaccinated"
          for (d in 1:num_vax_doses){
            workshop$dose[workshop$temp %in% c((T*(d-1)+1)*J+1):((T*d+1)*J)] = d
            for (t in 1:num_vax_types){
              workshop$vaccine_type[workshop$temp %in% c((((t-1)+(d-1)*T+1)*J+1):(((t-1)+(d-1)*T+2)*J))] = vax_type_list[t]
            }
          }
          workshop$risk_group = 'general_public'
          if (RISK>1){
            workshop$risk_group[(num_age_groups*num_vax_classes+1):(num_age_groups*num_vax_classes*2)] = risk_group_name
          }
          prev_state = rbind(prev_state,workshop)
        }
        prev_state$pop  = as.numeric(prev_state$pop)
        if (round(sum(prev_state$pop))!= sum(pop)){stop('prev state not equal to pop size! (~line 105 in time step)')}
        
        next_state=prev_state # initialise next state
        
        
      ### Include today's vaccinations
      for (r in 1:RISK){
        this_risk_group = risk_group_labels[r]   
         for (t in 1:num_vax_types){ #iterating over vaccine types
           this_vax = vax_type_list[t]
           
           this_vax_history = vaccination_history_FINAL[vaccination_history_FINAL$vaccine_type == this_vax & vaccination_history_FINAL$risk_group == this_risk_group,]
           
           # (1/3) recorded vax
           VR_this_step = crossing(dose = seq(1:D),
                                   age_group = age_group_labels,
                                   doses = 0)
           for (d in 1:D){
             for (i in 2:J){ #ASSUMPTION that don't vaccinate 0-4  
               if (nrow(this_vax_history[this_vax_history$dose == d & this_vax_history$date == as.Date(date_now) - vaxCovDelay$delay[vaxCovDelay$dose == d],]) >0){
                 VR_this_step$doses[VR_this_step$dose == d & VR_this_step$age_group == age_group_labels[i]] =
                   this_vax_history$doses_delivered_this_date[this_vax_history$date ==  as.Date(date_now) - vaxCovDelay$delay[vaxCovDelay$dose == d] & 
                                                                        this_vax_history$dose==d &
                                                                        this_vax_history$age_group == age_group_labels[i]]
      
               }
             }
           }

            for (i in 1:num_age_groups){ # across age groups
              
              increase = rep(0,num_vax_doses)
              for (d in 1:D){
                increase[d] = VR_this_step$doses[VR_this_step$dose == d & VR_this_step$age_group == age_group_labels[i]] 
              }
              
             for (j in 1:4){ #ASSUMPTION: all SEIR vaccinated
               class=class_name_list[j]
               
               prop = rep(0,num_vax_doses) #prop in S,E,I or R in vaccine groups
               for (d in 0:D){
                 if (d==0){
                   prop[d+1] = prev_state$pop[prev_state$class == class & prev_state$risk_group == this_risk_group & prev_state$vaccine_type == "unvaccinated" & prev_state$age_group == age_group_labels[i]]/
                     sum(prev_state$pop[prev_state$vaccine_type == "unvaccinated" & prev_state$risk_group == this_risk_group & prev_state$age_group == age_group_labels[i]])
                 } else{
                   prop[d+1] = prev_state$pop[prev_state$class == class & prev_state$risk_group == this_risk_group & prev_state$vaccine_type == this_vax &  prev_state$dose == d & prev_state$age_group == age_group_labels[i]]/
                     sum(prev_state$pop[prev_state$risk_group == this_risk_group & prev_state$vaccine_type == this_vax & prev_state$dose == d & prev_state$age_group == age_group_labels[i]])
                 }
                 if (is.nan(prop[d+1]) == TRUE){prop[d+1]=0}
               }

               next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == "unvaccinated" & next_state$age_group == age_group_labels[i]] =
                 next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == "unvaccinated" & next_state$age_group == age_group_labels[i]] - increase[1]* prop[1]
               
               for (d in 1:(D-1)){
                 next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == d & next_state$age_group == age_group_labels[i]] =
                   next_state$pop[ next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == d & next_state$age_group == age_group_labels[i]] + increase[d]*prop[d]-increase[d+1]*prop[d+1]
               }
               for (d in D){
                 next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == D & next_state$age_group == age_group_labels[i]] =
                   next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == D & next_state$age_group == age_group_labels[i]] + increase[D] * prop[D]
               }
             }
           }
         }
      }
      
        
        
      #### BOOSTER TO PREVIOUS PRIMARY
      #NB: using '8' as a flag for a booster to a previously primary delivered individaul
      if (nrow(vaccination_history_FINAL[vaccination_history_FINAL$dose == 8,])>0){
        for (r in 1:RISK){
          this_risk_group = risk_group_labels[r]   
          for (t in 1:num_vax_types){ #iterating over vaccine types FROM
            for (d in 1:D){
              this_vax = vax_type_list[t]
              
              this_vax_history = vaccination_history_FINAL %>%
                filter(FROM_vaccine_type == this_vax & FROM_dose == d & risk_group == this_risk_group)
              
              if (nrow(this_vax_history)>0){
                if (nrow(this_vax_history[this_vax_history$dose != '8',])>0){stop('incorrect dose flagged as a booster')}
                
                # (1/3) recorded vax
                VR_this_step = crossing(age_group = age_group_labels,
                                        doses = 0)
                for (i in 2:J){ #COMEBACK - could be faster with less for loop, assumption that don't vaccinate 0-4  
                  if (nrow(this_vax_history[this_vax_history$date == as.Date(date_now) - vaxCovDelay$delay[vaxCovDelay$dose == booster_dose_number],]) >0){
                    VR_this_step$doses[VR_this_step$age_group == age_group_labels[i]] =
                      this_vax_history$doses_delivered_this_date[this_vax_history$date ==  as.Date(date_now) - vaxCovDelay$delay[vaxCovDelay$dose == booster_dose_number] &
                                                                   this_vax_history$age_group == age_group_labels[i]]
                  }
                }
                
                for (i in 1:num_age_groups){ # across age groups
                  increase = VR_this_step$doses[VR_this_step$age_group == age_group_labels[i]] 
                  
                  for (j in 1:4){ #let's assume all SEIR vaccinated
                    class=class_name_list[j]
                    
                    prop = prev_state$pop[prev_state$class == class & prev_state$risk_group == this_risk_group & prev_state$vaccine_type == this_vax &  prev_state$dose == d & prev_state$age_group == age_group_labels[i]]/
                      sum(prev_state$pop[prev_state$risk_group == this_risk_group & prev_state$vaccine_type == this_vax & prev_state$dose == d & prev_state$age_group == age_group_labels[i]])
                    if (is.nan(prop) == TRUE){prop=0}
                    
                    next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == booster_type & next_state$dose == booster_dose_number & next_state$age_group == age_group_labels[i]] =
                      next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == booster_type & next_state$dose == booster_dose_number & next_state$age_group == age_group_labels[i]] + increase * prop
                    
                    next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == d & next_state$age_group == age_group_labels[i]] =
                      next_state$pop[next_state$class == class & next_state$risk_group == this_risk_group & next_state$vaccine_type == this_vax & next_state$dose == d & next_state$age_group == age_group_labels[i]] - increase * prop
                    
                  }
                }
              }
            }
          }
        }
      }
        rm(this_vax_history, VR_this_step, this_risk_group, this_vax, increase, class, prop)
        
        if (fitting == "on" & date_now == as.Date('2021-11-14')){next_state_FIT = next_state} #savings to compare against known point of seroprevalence
        if (! date_start %in% seed_date & date_now %in% seed_date){
          if (fitting == "on"){
            if (date_now == seed_date[1]){
              strain_now = 'delta'
            } else if (date_now == seed_date[2]){
              strain_now = 'omicron'
              parameters$lambda = 1/2.22 #COMEBACK - hard coded :(
              parameters$delta = 1/9.87
            }
            parameters$beta = rep(beta_fitted_values$beta_optimised[beta_fitted_values$strain == strain_now],num_age_groups)
          } 
          
          seed.Infected = seed*AverageSymptomaticPeriod/(AverageSymptomaticPeriod+AverageLatentPeriod)
          seed.Exposed  = seed*AverageLatentPeriod/(AverageSymptomaticPeriod+AverageLatentPeriod)
    
          #ASSUMPTION: uniform across age groups
          seed.Infected = round(seed.Infected * pop/sum(pop))
          seed.Exposed = round(seed.Exposed * pop/sum(pop))
    
          for (i in 1:num_age_groups){ # across age groups
            infect = seed.Infected[i]
            expose = seed.Exposed[i]
            
            infect = round(infect * next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]] /sum(next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]]))
            expose = round(expose * next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]] /sum(next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]]))
            
            next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]] = next_state$pop[next_state$class == 'S'  & next_state$age_group == age_group_labels[i]] - infect - expose
            next_state$pop[next_state$class == 'I'  & next_state$age_group == age_group_labels[i]] = next_state$pop[next_state$class == 'I'  & next_state$age_group == age_group_labels[i]] + infect
            next_state$pop[next_state$class == 'E'  & next_state$age_group == age_group_labels[i]] = next_state$pop[next_state$class == 'E'  & next_state$age_group == age_group_labels[i]] + expose
            
          }
        }
        if (fitting == "off"){
          if (date_now>=seed_date & (! as.Date('1900-01-01') %in% seed_date)){
            parameters$VE$VE = parameters$VE$VE * 0.9
            parameters$rho = parameters$rho * 0.9
          }
        }
    
        if (round(sum(next_state$pop))!= round(sum(prev_state$pop))){stop('pop not retained between next_state and prev_state!')}
        if (round(sum(next_state$pop))!= sum(pop)) {stop('pop in next_state not equal to setting population')}
        if(nrow(next_state[round(next_state$pop)<0,])>0){ stop('(4)_time_step line ~290')}
        
        
        # convert back into long vector form for ODE solver
        workshop = next_state
        workshop$class = factor(workshop$class, levels = disease_class_list)
        workshop$risk_group = factor(workshop$risk_group, levels = risk_group_labels)
        workshop$vaccine_type = factor(workshop$vaccine_type, levels = vax_type_list)
        workshop$age_group = factor(workshop$age_group, levels = age_group_labels)
        
        next_state = workshop %>% arrange(class,risk_group,dose,vaccine_type,age_group)
        
        S_next=next_state$pop[next_state$class == 'S']
        E_next=next_state$pop[next_state$class == 'E']
        I_next=next_state$pop[next_state$class == 'I']
        R_next=next_state$pop[next_state$class == 'R']
        
        next_state_FINAL=as.numeric(c(S_next,E_next,I_next,R_next,
                                      Incidence_inital,Exposed_incidence_inital)) #setting Incid to repeated 0s
        rm(S_next,E_next,I_next,R_next)
    }     
    
    # next time_step!
    sol <- as.data.frame(ode(y=next_state_FINAL,times=(seq(0,time_step,by=1)),func=covidODE,parms=parameters))
    
    sol[,1]=sol[,1]+time_step*(increments_number-1) #make times correct
    
    if (increments_number>1){ sol_log=head(sol_log,-1) }#remove last entry from sol_log (overlap of runs)
    sol_log <- rbind(sol_log,sol)
    sol_log_unedited <- rbind(sol_log_unedited,sol)


  ### INCIDENCE CALCULATIONS 
  J=num_age_groups
  T=num_vax_types
  D=num_vax_doses
  RISK=num_risk_groups
  A=RISK*J*(T*D+1) # +1 is unvax
  
  incidence_log_unedited <- sol_log_unedited[, c(1,(A*num_disease_classes+2):(A*(num_disease_classes+1)+1))]
  incidence_log_unedited <- incidence_log_unedited %>% filter (time %% time_step == 0, rowSums(incidence_log_unedited) != time)
  incidence_log_unedited <- distinct(round(incidence_log_unedited,digits=2))
  
  incidence_log_unedited$date <- date_start + incidence_log_unedited$time
  incidence_log_unedited$daily_cases  <- rowSums(incidence_log_unedited[,2:(A+1)])
  
  incidence_log <- incidence_log_unedited %>% select(date,daily_cases) 
  if (! fitting == "on"){incidence_log = rbind(fitted_incidence_log,incidence_log)}
  
  incidence_log = incidence_log %>%
    mutate(rolling_average = (daily_cases + lag(daily_cases) + lag(daily_cases,n=2)+lag(daily_cases,n=3)
                              +lag(daily_cases,n=4)+lag(daily_cases,n=5)+lag(daily_cases,n=6))/7,
           rolling_average_percentage = 100*rolling_average/sum(pop),
           cumulative_incidence = cumsum(daily_cases),
           cumulative_incidence_percentage = 100*cumsum(daily_cases)/sum(pop))
  
  if (fitting == "on"){
    if (fitting == "off" & increments_number == 1){
    } else{
      Reff <- Reff_time_step(parameters,next_state)
      Reff_tracker = rbind(Reff_tracker,Reff)
    }
    
    rho_tracker_dataframe = rbind(rho_tracker_dataframe,parameters$rho) 
    
    workshop = parameters$VE
    workshop = workshop[workshop$VE>0,] %>% mutate(immunity = paste(vaccine_type,dose))
    if (nrow(workshop)>0){
      workshop = aggregate(workshop$VE, by=list(category=workshop$immunity), FUN=mean)
      colnames(workshop) = c('dose','VE')
      workshop$date = date_now
      VE_tracker_dataframe = rbind(VE_tracker_dataframe,workshop)
    }
  }
  }
} ### END INCREMENT (#incidence log moved within loop to allow rho_time_step to access)
if (fitting == "off"){
  rm(fitted_incidence_log, sol_log, covidODE, rho_time_step, Reff_time_step, NPI, NGM_R0)
}

check <- sol_log_unedited
check$Incid = rowSums(sol_log_unedited[,(A*4+2):(A*5+1)])
check$R = rowSums(sol_log_unedited[,(A*3+2):(A*4+1)])
check$I = rowSums(sol_log_unedited[,(A*2+2):(A*3+1)])
check$E = rowSums(sol_log_unedited[,(A+2):(A*2+1)])
check$S = rowSums(sol_log_unedited[,2:(A+1)])

check = check %>% select(S,E,I,R,Incid) %>%
  mutate(pop = S + E + I + R)
if (round(check$pop[1]) !=sum(pop)){stop('population does not stay constant!')}


### INCIDENCE LOG TIDY 
workshop = subset(incidence_log_unedited, select=-c(time,daily_cases))
workshop = pivot_longer(
  workshop,
  cols = paste((num_disease_classes)*(num_age_groups*num_vax_classes)*RISK+1):paste((num_disease_classes+1)*(num_age_groups*num_vax_classes)*RISK),
  names_to = 'temp',
  values_to = 'incidence'
)
workshop$temp = as.numeric(workshop$temp) - (num_disease_classes)*(num_age_groups*num_vax_classes)*RISK

workshop2=as.data.frame(unique(workshop$temp)); colnames(workshop2)=c('temp')
if (RISK == 1){
  workshop2 = workshop2 %>%   mutate(temp_risk = temp, risk_group = risk_group_labels[[1]])
} else {
  workshop2 = workshop2 %>%
    mutate(temp_risk = case_when(
      temp <= max(workshop$temp)/2 ~ temp,
      temp > max(workshop$temp)/2  ~ temp - max(workshop$temp)/2 
    ),
    risk_group = case_when(
      temp <= max(workshop$temp)/2 ~ risk_group_labels[[1]],
      temp > max(workshop$temp)/2  ~ risk_group_labels[[2]]
    ))
}
workshop2$age_group = rep(age_group_labels,num_vax_classes) #smallest subdivision is age
workshop2$dose = 0                                          #then dose
workshop2$vaccine_type = "unvaccinated"                     #then vaccine type
for (d in 1:num_vax_doses){
  workshop2$dose[workshop2$temp_risk %in% c((T*(d-1)+1)*J+1):((T*d+1)*J)] = d
  for (t in 1:num_vax_types){
    workshop2$vaccine_type[workshop2$temp_risk %in% c((((t-1)+(d-1)*T+1)*J+1):(((t-1)+(d-1)*T+2)*J))] = vax_type_list[t]
  }
}
#View(workshop2) #CHECKED: yes aligns as expected

incidence_log_tidy = workshop %>% 
  left_join(workshop2, by = "temp")
incidence_log_tidy = subset(incidence_log_tidy,select=-c(temp))

if (! fitting == "on"){
  incidence_log_tidy = rbind(fitted_incidence_log_tidy,incidence_log_tidy)
  rm(fitted_incidence_log_tidy)
}



### EXPOSED LOG TIDY -> reinfection ratio (reduced severe outcomes)
skip = (num_disease_classes+1)*(num_age_groups*num_vax_classes)*RISK
exposed_log = sol_log_unedited %>% 
  select(1, (skip + 2):(skip + 2*J + 1))
rm(sol_log_unedited)

exposed_log = exposed_log %>%
  filter (time %% time_step == 0, rowSums(exposed_log) != time)
  
exposed_log$date <- date_start + exposed_log$time

workshop = subset(exposed_log, select = -c(time))
workshop = pivot_longer(
  workshop,
  cols = colnames(workshop)[1]:colnames(workshop)[2*J],
  names_to = 'temp',
  values_to = 'exposed'
)
workshop$temp = as.numeric(workshop$temp) - skip

workshop2=as.data.frame(unique(workshop$temp)); colnames(workshop2)=c('temp')
workshop2$age_group = rep(age_group_labels,2) #smallest subdivision is age
workshop2 = workshop2 %>% mutate(infection_type = case_when(
  temp <= num_age_groups ~ "new_infection",
  temp > num_age_groups ~ "reinfection"))

exposed_log_tidy = workshop %>% 
  left_join(workshop2, by = "temp") %>%
  select(-temp)

exposed_log = exposed_log_tidy %>% ungroup() %>% pivot_wider(
  id_cols = c(date,age_group),
  names_from = infection_type,
  values_from = exposed) %>% 
  mutate(reinfection_ratio = reinfection/(new_infection+reinfection))
#ggplot(exposed_log) + geom_line(aes(x=date,y=reinfection_ratio,color=as.factor(age_group)))


if ( fitting == "on"){
  colnames(rho_tracker_dataframe) = c('rho')
  rho_tracker_dataframe = cbind(rho = rho_tracker_dataframe, date = seq(date_start+1,date_start+nrow(rho_tracker_dataframe),by="days"))
  if (fitting == "on"){Reff_tracker <- cbind(Reff = Reff_tracker, date = seq(date_start+1,date_start+nrow(rho_tracker_dataframe)+1,by="days"))}
  colnames(Reff_tracker) = c('Reff','date')
}


