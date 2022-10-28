### This function calculates the incidence of severe outcomes in our Study setting

#LIMITATIONS:
# - age distribution from Qatar (where care seeking patterns may differ, underlying prevalence of risk factors may vary)
# - age distribution from WT, may have changed with Omicron!
# - age distribution adapted from 10 age year distributions, hence don't include detail on 0-4

##### setup
if (exists("age_split_results") == FALSE){ age_split_results = "N"}
plot_standard =     theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))
if (as.Date('1900-01-01') %in% seed_date){
  use_date = date_start +1
} else{
  use_date = seed_date[1]
}


##### infection-derived immunity against severe outcomes
rho_SO_est = unique(rho_dn$protection[rho_dn$outcome == 'severe_outcome'])
if (length(rho_SO_est)>1){stop('rho against severe outcome not constant as currently assumed, see (function)_severe_outcome_proj')}

reinfection_protection = exposed_log %>%
  mutate(protection = reinfection_ratio * rho_SO_est) %>%
  select(date,age_group,protection)
#ggplot(reinfection_protection) + geom_point(aes(x=date,y=protection,color=as.factor(age_group)))


##### projection from incidence to severe outcomes
#(A/D) Join incidence_log_tidy with severe outcome incidence by vax status
workshop = severe_outcome_this_run %>%
  left_join(incidence_log_tidy, by = c("date", "age_group", "risk_group", "vaccine_type", "dose")) %>%
  mutate(proj = incidence*percentage) %>%
  left_join(reinfection_protection, by = c("date", "age_group")) %>%
  mutate(proj = proj*(1-protection))
if(!nrow(severe_outcome_this_run[severe_outcome_this_run$date <= max(incidence_log_tidy$date),]) == nrow(workshop)){stop('have lost rows in severe outcome proj')} 

severe_outcome_log_tidy = workshop %>% 
  select(date,risk_group,age_group,dose,vaccine_type,outcome,proj) %>%
  filter(date >= date_start+1)


#(B/D) Sum across age groups, doses and vaccination status to get overall severe incidence per day
if (age_split_results == "N"){
  workshop_2 = workshop %>%
    group_by(date,outcome) %>%
    summarise(proj=sum(proj),.groups = "keep")
 
   
  #(C/D) Add cases as an outcome
  workshop_incid =  incidence_log_unedited %>%
    mutate(outcome ='cases',proj = daily_cases) %>%
    select(date,outcome,proj)
  rm(incidence_log_unedited)
  workshop_2 = rbind(workshop_2,workshop_incid)
  
  #(D/D) Calculate cumulative severe outcomes by outcome type
  severe_outcome_log = workshop_2 %>% 
    filter(date > use_date) %>%
    group_by(outcome) %>%
    mutate(proj_cum = cumsum(proj))
  
  severe_outcome_log_plot = workshop_2 %>% 
    filter(date >= date_start) %>%
    group_by(outcome) %>%
    mutate(proj = case_when(
      is.na(proj) ~ 0,
      TRUE ~ proj),
      proj_cum = cumsum(proj)) %>%
    filter(! outcome == 'cases' & !date %in% seed_date)
  
  plot1 <- 
    ggplot() + 
    geom_line(data=severe_outcome_log_plot[severe_outcome_log_plot$outcome != 'cases',],
              aes(x=date,y=proj,color=as.factor(outcome)),na.rm=TRUE) +
    xlab("") + 
    scale_x_date(date_breaks="1 month", date_labels="%b") +
    #ylim(0,40) +
    labs(color="")+
    ylab("incidence") + 
    plot_standard

  plot2 <- 
    ggplot() + 
    geom_line(data=severe_outcome_log_plot[severe_outcome_log_plot$outcome != 'cases',],aes(x=date,y=proj_cum,color=as.factor(outcome)),na.rm=TRUE) +
    xlab("") + 
    scale_x_date(date_breaks="1 month", date_labels="%b") +
    ylab("cumulative incidence") +
    plot_standard + 
    labs(color="")
  
  grid.arrange(plot1, plot2)
  

} else if (age_split_results == "Y"){
  workshop_2 = workshop %>%
    mutate(macro_age_group = case_when(
      age_group %in% c('0 to 4') ~ 'children <5',
      age_group %in% c('5 to 9','10 to 17') ~ 'children 5-17',
      TRUE ~ 'adults'
    )) %>%
    group_by(date,macro_age_group,outcome) %>%
    summarise(proj=sum(proj),.groups = "keep")
  
  #(C/D) Add cases as an outcome
  workshop_incid =  incidence_log_tidy[,c('date','incidence','age_group')] %>%
    mutate(outcome ='cases') %>%
    mutate(macro_age_group = case_when(
      age_group %in% c('0 to 4') ~ 'children <5',
      age_group %in% c('5 to 9','10 to 17') ~ 'children 5-17',
      TRUE ~ 'adults'
    )) %>%
    group_by(date,macro_age_group,outcome) %>%
    summarise(proj=sum(incidence),.groups = "keep") %>%
    select(date,macro_age_group,outcome,proj)
  
  workshop_2 = rbind(workshop_2,workshop_incid)
  
  #(D/D) Calculate cumulative severe outcomes by outcome type
  severe_outcome_log = workshop_2 %>% 
    filter(date > use_date)%>%
    group_by(outcome,macro_age_group) %>%
    mutate(proj_cum = cumsum(proj))
  
  severe_outcome_log_plot = workshop_2 %>% 
    filter(date >= date_start) %>%
    group_by(outcome,macro_age_group) %>%
    mutate(proj_cum = cumsum(proj))
    
  plot_children_0_4 <- 
    ggplot() + 
    geom_line(data=severe_outcome_log_plot[severe_outcome_log_plot$outcome != 'cases' & severe_outcome_log_plot$macro_age_group == 'children <5',]
              ,aes(x=date,y=proj,color=as.factor(outcome)),na.rm=TRUE) +
    labs(title='children 0-4')+
    plot_standard +
    xlab("") +
    ylab("incidence") +
    labs(color="") +
    scale_x_date(date_breaks="1 month", date_labels="%b")
  
  plot_children_5_17 <- 
    ggplot() + 
    geom_line(data=severe_outcome_log_plot[severe_outcome_log_plot$outcome != 'cases' & severe_outcome_log_plot$macro_age_group == 'children 5-17',]
              ,aes(x=date,y=proj,color=as.factor(outcome)),na.rm=TRUE) +
    labs(title='children 5-17')+
    plot_standard +
      xlab("") +
      ylab("incidence") +
      labs(color="") +
      scale_x_date(date_breaks="1 month", date_labels="%b")
  
  plot_adults <- 
    ggplot() + 
    geom_line(data=severe_outcome_log_plot[severe_outcome_log_plot$outcome != 'cases' & severe_outcome_log_plot$macro_age_group == 'adults',]
              ,aes(x=date,y=proj,color=as.factor(outcome)),na.rm=TRUE) +
    labs(title='adults')+
    plot_standard +
    xlab("") +
    ylab("incidence") +
    labs(color="") +
    scale_x_date(date_breaks="1 month", date_labels="%b")
 
  grid.arrange(plot_children_0_4, plot_children_5_17, plot_adults)
  
}

#create row for table comparing vaccine strategies
row = severe_outcome_log %>% 
  filter(date == max(severe_outcome_log$date)) %>%
  select(-c(proj,date)) %>%
  pivot_wider(names_from=outcome,
              values_from=proj_cum) 