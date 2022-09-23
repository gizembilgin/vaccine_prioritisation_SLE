### This (mech shop) applies age distribution to severe outcomes estimates.
### Creates: age_dn_severe_outcomes

age_dn_severe_outcomes = data.frame()
for (setting in c('SLE')){
  
  ### (1) Importing raw age distributions ########################################################
  #(A) importing raw values from Seedat et al. (Ayoub et al underlying data)
  workshop = read.csv('1_inputs/severe_outcome_age_distribution_RAW.csv')
  
  #(B) join overall value back on to age-specific values
  overall = workshop %>%
    filter(agegroup == 'overall') %>% 
    mutate(overall=rate) %>% 
    select(outcome,overall)
  workshop = workshop %>%
    filter(agegroup != 'overall')  %>%
    left_join(overall) %>% 
    select(agegroup,outcome,rate,overall) 
  
  #(C) Check zero values, and manually overwrite them
  workshop[workshop$rate == 0,] #n=3
  workshop$rate[workshop$agegroup == '0 to 9' & workshop$outcome=='critical_disease'] = 
    workshop$rate[workshop$agegroup == '10 to 19' & workshop$outcom=='critical_disease'] * 1.5
  workshop$rate[workshop$agegroup == '10 to 19' & workshop$outcom=='death'] = 
    workshop$rate[workshop$agegroup == '20 to 29' & workshop$outcom=='death']
  workshop$rate[workshop$agegroup == '0 to 9' & workshop$outcom=='death'] = 
    workshop$rate[workshop$agegroup == '10 to 19' & workshop$outcom=='death'] * 1.5
  workshop[workshop$rate == 0,] #n=0
  
  #(D) calculate RR, clean age groups
  workshop = workshop %>%
    mutate(RR = rate/overall) %>%
    filter(agegroup != '80+') #remove dodgy age group!
  workshop$agegroup[workshop$agegroup == '70 to 79'] = '70 to 100'
  
  #(E) plot to see!
  plot_list = list()
  for (i in 1:length(unique(workshop$outcome))){
    outcome = unique(workshop$outcome)[i]
    plot_list [[i]] <- ggplot() + 
      geom_point(data=workshop[workshop$outcome==outcome,],aes(x=agegroup,y=RR)) +
      labs(title=paste("distribution of all ",outcome,sep=""))
  }
  gridExtra::grid.arrange(grobs=plot_list)
  
  plot_list = list()
  for (i in 1:length(unique(workshop$outcome))){
    outcome = unique(workshop$outcome)[i]
    plot_list [[i]] <- ggplot() + 
      geom_point(data=workshop[workshop$outcome==outcome,],aes(x=agegroup,y=RR)) +
      labs(title=paste("distribution of all ",outcome,"(log-scale)",sep=""))+ scale_y_log10()
  }
  gridExtra::grid.arrange(grobs=plot_list)
  ########################################################################################################
  
  
  
  ### (2) Adjust age-distributions to setting ###########################################################
  #(A) Multiply RR by pop-level value
  pop_level <- read.csv('1_inputs/severe_outcome_country_level.csv')
  pop_level = pop_level %>%
    filter(country == setting) %>%
    mutate(pop_est = percentage) %>%
    select(outcome,pop_est)
  workshop = workshop %>% left_join(pop_level) %>%
    mutate(age_est = pop_est * RR)
  
  #(B) Calculate % of pop by 10 year age bands
  age_groups_10 = c(0,9,19,29,39,49,59,69,100)
  age_group_labels_10 = c('0 to 9','10 to 19','20 to 29','30 to 39','40 to 49','50 to 59','60 to 69','70 to 100')
  
  pop_10_bands <- read.csv(paste(rootpath,"inputs/pop_estimates.csv",sep=''), header=TRUE)
  pop_10_bands <- pop_10_bands %>%
    filter(country == setting) %>%
    mutate(agegroup = cut(age,breaks = age_groups_10, include.lowest = T,labels = age_group_labels_10)) %>%
    group_by(agegroup) %>%
    summarise(pop = sum(population)) 
  pop_10_bands = pop_10_bands %>%  mutate(pop_percentage = pop/sum(pop_10_bands$pop))  
  
  
  #(C) Hence, rederive pop-level est by age_est * % pop
  workshop = workshop %>% left_join(pop_10_bands) %>%
    mutate(interim = pop_percentage * age_est) 
  workshop_sum = workshop %>%
    group_by(outcome) %>%
    summarise(interim_est = sum(interim))

  #(D) Calculate correction per outcome
  workshop_sum = workshop_sum %>% left_join(pop_level) %>%
    mutate(correction = pop_est/interim_est)
  
  #(E) Apply correction, and recalculate age- and pop- estimates
  workshop = workshop %>% left_join(workshop_sum) %>%
    mutate(RR2 = RR * correction,
           age_est2 = pop_est * RR2,
           interim2 = pop_percentage * age_est2)
  
  #(F) Check!
  workshop %>%
    group_by(outcome) %>%
    summarise(interim_est = sum(interim2)) %>% 
    left_join(pop_level) %>%
    mutate(correction = pop_est/interim_est)
  ########################################################################################################
  
  
  
  ### (3) Convert age-distributions to model age groups ###################################################
  #(A) remove redundant columns from workshop
  workshop = workshop %>% 
    select(agegroup,outcome,RR2) %>%
    rename(agegroup_10 = agegroup,
           RR = RR2)
  
  #(B) align 10-year age bands with model age bands 
  age_groups_10 = c(0,9,19,29,39,49,59,69,100)
  age_group_labels_10 = c('0 to 9','10 to 19','20 to 29','30 to 39','40 to 49','50 to 59','60 to 69','70 to 100')
  
  pop_w <- read.csv(paste(rootpath,"inputs/pop_estimates.csv",sep=''), header=TRUE)
  pop_w <- pop_w[pop_w$country == setting,] %>%
    mutate(agegroup_10 = cut(age,breaks = age_groups_10, include.lowest = T,labels = age_group_labels_10),
           agegroup_model = cut(age,breaks = age_groups_num, include.lowest = T,labels = age_group_labels)) %>%
    ungroup() %>%
    group_by(agegroup_model) %>%
    mutate(group_percent = population/sum(population)) %>%
    select(age,agegroup_10,agegroup_model,group_percent)
  
  #(C) calculate model age group specific RR by weighted within group %
  workshop = workshop %>% left_join(pop_w) %>%
    mutate(interim = RR*group_percent)
  workshop_sum =  workshop %>% 
    group_by(outcome,agegroup_model) %>%
    summarise(RR = sum(interim)) %>%
    rename(age_group = agegroup_model)
  
  #(D) CHECK
  #find pop % within model agroups
  workshop_sum = workshop_sum %>% left_join(pop_setting) %>%
    mutate(pop_percent = pop/sum(pop_setting$pop),
           interim =RR*pop_percent)
  if(sum(workshop_sum$pop_percent) != length(unique(workshop_sum$outcome))){stop('ah!')}
  
  aggregate(workshop_sum$interim, by=list(category= workshop_sum$outcome), FUN=sum)
  #1! perfect!
  #######################################################################################################
  
  
  ###(4) Finalise
  age_distribution_save = workshop_sum %>%
    select(outcome,age_group,RR) %>%
    mutate(setting=setting)
  age_dn_severe_outcomes = rbind(age_dn_severe_outcomes,age_distribution_save)
}


#plot settings
include = c('severe_disease','hosp','death')
workshop = age_dn_severe_outcomes [age_dn_severe_outcomes$outcome %in% include,]
plot_list = list()
for (i in 1:length(unique(workshop$outcome))){
  outcome = unique(workshop$outcome)[i]
  plot_list [[i]] <- ggplot() + 
    #geom_point(data=workshop[workshop$outcome==outcome,],aes(x=age_group,y=RR,color=as.factor(setting))) +
    geom_point(data=workshop[workshop$outcome==outcome,],aes(x=age_group,y=RR)) +
    labs(title=paste(outcome,sep=""))
}
gridExtra::grid.arrange(grobs=plot_list)

plot_list = list()
for (i in 1:length(unique(workshop$outcome))){
  outcome = unique(workshop$outcome)[i]
  plot_list [[i]] <- ggplot() + 
    #geom_point(data=workshop[workshop$outcome==outcome,],aes(x=age_group,y=RR,color=as.factor(setting))) +
    geom_point(data=workshop[workshop$outcome==outcome,],aes(x=age_group,y=RR)) +
    labs(title=paste(outcome,"(log-scale)",sep=""))+ scale_y_log10()
}
gridExtra::grid.arrange(grobs=plot_list)


save(age_dn_severe_outcomes, file = "1_inputs/severe_outcome_age_distribution.Rdata")
rm(workshop_sum,workshop,age_group_labels_10,age_groups_10,pop_w,overall,plot_list)
