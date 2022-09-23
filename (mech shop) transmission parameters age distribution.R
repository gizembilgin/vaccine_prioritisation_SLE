### This (mech shop) translates transmission parameter estimates to model age groups.

raw = read.csv("1_inputs/model_param_raw.csv",header=TRUE)
raw = raw %>% rename(agegroup_RAW = age_group)
underlying_age_grouping <- c(0,9,19,29,39,49,59,69,110)

setting_list = c('SLE')
param_age = data.frame()

for (i in 1:length(setting_list)){
  setting = setting_list[i]
  
  pop_RAW =  pop_orig %>%
    filter(country == setting) %>%
    mutate(agegroup_RAW = cut(age,breaks = underlying_age_grouping, include.lowest = T, labels = unique(raw$agegroup_RAW)),
           agegroup_MODEL = cut(age,breaks = age_groups_num, include.lowest = T, labels = age_group_labels)) %>%
    ungroup() %>%
    group_by(agegroup_MODEL) %>%
    mutate(model_group_percent = population/sum(population))
  
  workshop = pop_RAW %>% left_join(raw) %>% 
    mutate(interim = model_group_percent * value) %>%
    group_by(param,agegroup_MODEL) %>%
    summarise(value = sum(interim)) %>%
    rename(agegroup = agegroup_MODEL)

  param_age = rbind(param_age,workshop)
}


save(param_age, file = "1_inputs/param_age.Rdata")

