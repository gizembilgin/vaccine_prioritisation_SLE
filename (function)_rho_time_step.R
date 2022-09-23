### This script loads the data and function needed for dynamic infection-derived protection against reinfection.
### Dynamic infection-derived protection (rho) requires:
# (1) the shape of waning rho
# (2) tracking of days since infection of people in Recovered class
# (3) tracking of the introduction of new immune escape variants
# LIMITATION: our model does not track multiple infections, and hence assumes no additional immunity from repeat infections



### (1) Shape of waning ###############################################################################
#Exponential fit separately, see Supplementary Material for data source
rho_dn = crossing(type = c('new_to_new','prev_to_new'),
                  outcome = c('symptomatic_disease','severe_outcome'),
                  days = seq(0,365*2)
)
rho_dn = rho_dn %>% mutate(
  protection = case_when(
    outcome == 'symptomatic_disease' & type == 'prev_to_new' ~ 0.9932*exp(-0.004*days),
    outcome == 'symptomatic_disease' & type == 'new_to_new' ~ 0.9654*exp(-0.0002*days),
    outcome == 'severe_outcome' ~ 0.88 
  )
)

rho_est_plot = ggplot() + 
  geom_line(data = rho_dn[rho_dn$outcome == 'symptomatic_disease',],aes(days,protection,color=as.factor(type))) + 
  ylim(0,1) +labs(title = 'symptomatic_disease') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())   

rho_dn_wide = rho_dn %>% 
  filter(outcome == 'symptomatic_disease') %>%
  pivot_wider(
    id_cols = days,
    names_from = type,
    values_from = protection
  )
#######################################################################################################



### (2) Shape of introduction ###############################################################################
intro_raw <- read.csv("1_inputs/GISAID_omicron_intro_Africa.csv", header=TRUE)
intro_raw$date = as.Date(intro_raw$date, '%d/%m/%Y')

intro_raw = intro_raw %>% 
  mutate(days = as.numeric(date - min(intro_raw$date)),
         percentage = percentage/100)

ggplot(intro_raw) + geom_point(aes(x=days,y=percentage))

days_to_predict = seq(0,365)
smoothed_spline <- smooth.spline(x = intro_raw$days, y = intro_raw$percentage, df = 9)
fitted.results = predict(smoothed_spline,days_to_predict,deriv = 0)
fit = data.frame(fitted.results)
colnames(fit) = c('days','percentage')

ggplot() + geom_point(data = intro_raw, aes(x=days,y=percentage)) +
  geom_line(data = fit, aes(x=days,y=percentage)) 

synthetic_strain_shift = fit %>% filter(percentage >=0 & percentage <=1) 
synthetic_strain_shift$days = synthetic_strain_shift$days - min(synthetic_strain_shift$days)

ggplot() + geom_point(data = intro_raw, aes(x=days,y=percentage)) +
  geom_line(data = synthetic_strain_shift, aes(x=days,y=percentage)) +
  ylab('percentage of circulating strains') +
  xlab('days since introduction')+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = 'black'))
#######################################################################################################



### (3) Tracking of days since infection (proxy) for people in Recovered class ########################
rho_acq = rho_dn %>% filter(outcome == 'symptomatic_disease')

# This function calculates the population-level infection-derived protection against reinfection at any given time by age_group
rho_time_step <- function(date_now){
  
  #create data set of historical cases by day in recovery class
  if (exists("incidence_log") == FALSE){ incidence_log = data.frame() }
  if (nrow(incidence_log)>0){
    workshop = incidence_log %>% select(date,daily_cases)
    workshop = rbind(workshop,hist_cases)
  } else{
    workshop = hist_cases
  }
  #______________________
  
  
  #calculate proportion of individuals in recovery class by day since recovery
  workshop = workshop %>%
    mutate(date = date + round(AverageSymptomaticPeriod)) %>%
    filter(date <= date_now & date > (date_now - 1/omega)) %>%
    mutate(days = round(as.numeric(date_now - date)))
  workshop$date = as.Date(workshop$date, '%Y-%m-%d')
  workshop = workshop %>% mutate(prop_window = daily_cases/sum(workshop$daily_cases))   
  if (round(sum(workshop$prop_window),digits=5) != 1){stop('error in rho_time_step')}
  
  ggplot(workshop) + geom_line(aes(date,prop_window)) 
  #______________________
  
  
  #calculate protection to latest versus previous strains
  if (length(variant_change_date[variant_change_date<date_now])>0){ # do we have any new variants?
    
    num_variants_introduced = length(variant_change_date[variant_change_date<date_now])
    
    #rho = how protected are people today from the circulating strains today?
    #    = %new circulating today * protection to new today + %old circulating today * protection to old today
    
    for (variant_num in 1:num_variants_introduced){
      
      this_variant_introduction = variant_change_date[variant_num]
      this_variant_shift = synthetic_strain_shift %>% 
        select(days,percentage) %>%
        mutate(date = days + this_variant_introduction) %>%
        select(-days)
      this_variant_shift$date = as.Date(this_variant_shift$date, format =  '%Y-%m-%d')
      
      if (variant_num == 1){
        protection_to_old = workshop %>% 
          left_join(rho_acq[rho_acq$type == 'new_to_new',], by = "days") %>% 
          mutate(interim = protection * prop_window)
        protection_to_old = sum(protection_to_old$interim,na.rm=TRUE)
      } else{
        protection_to_old = protection
      }
      
      protection_to_new = workshop %>%
        left_join(rho_dn_wide, by = "days") %>% 
        left_join(this_variant_shift,by="date") %>%
        mutate(interim = case_when(
          date < this_variant_introduction ~ prev_to_new * prop_window,
          date >= this_variant_introduction & date <= max(this_variant_shift$date)~(prev_to_new*(1-percentage)+new_to_new*percentage)*prop_window,
          date > max(this_variant_shift$date) ~ new_to_new*prop_window
        ))
      protection_to_new = sum(protection_to_new$interim,na.rm=TRUE)
      
      if (date_now %in% this_variant_shift$date){
        percentage_new = this_variant_shift$percentage[this_variant_shift$date == date_now]
      } else{
        percentage_new = this_variant_shift$percentage[this_variant_shift$date == max(this_variant_shift$date)]
      }
      
      protection = percentage_new * protection_to_new + (1-percentage_new) * protection_to_old
    }
    
    function_result = protection
    
  } else{
    workshop = workshop %>% 
      left_join(rho_acq[rho_acq$type == 'new_to_new',], by = "days") %>%
      mutate(interim = protection * prop_window)
    
    function_result = sum(workshop$interim,na.rm=TRUE) 
  }
  
  return(function_result)
}