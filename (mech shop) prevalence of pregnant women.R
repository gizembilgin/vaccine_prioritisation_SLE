### This (mech shop) calculates the % per age group of pregnant women.
### Creates: prevalence_pregnancy

#LIMITATION: the model takes a one year snapshot, if a longer snapshot is required women will have to cycle in and out of this group, then:
# *3/4 to get women currently pregnant?  or *1/4 for third trimester? or * ~ 2 for all lactating women?


### read in data
ASFR = read.csv("1_inputs/DHS_ASFR.csv",header=TRUE)
women_pop = read.csv(paste(rootpath,"inputs/pop_estimates_female.csv",sep=''),header=TRUE)


### add 10-14 pregnancy as reported in DHS 2019 with retrospective data
row_10_14 = data.frame(' 10-14 ',4/1000,NA,NA)
colnames(row_10_14) = colnames(ASFR)
ASFR = rbind(row_10_14,ASFR)


### plot ASFR
#View(ASFR)
ggplot(data=ASFR) + 
  geom_pointrange(aes(x=ASFR*100,y=AGE,xmin=LCI*100,xmax=UCI*100)) +
 # xlim(0,1) +
  xlab("Age-specific fertility ratio (%)") + 
  ylab("") + 
  labs(title="") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = 'black'))


### calculate and plot female ratio estimates
#pop_setting_orig - colnames: age, country, population
#women_pop - colnames: age, country, population, population_thousands
women_pop = women_pop %>% 
  rename(pop_women = population) %>% 
  select(-population_thousands)
pop_together = pop_setting_orig %>% 
  left_join(women_pop) %>%
  mutate(female_prop = pop_women/population) %>%
  select(-pop_women)
ggplot(data=pop_together) + 
  geom_point(aes(x=female_prop*100,y=age)) +
  xlim(0,100) +
  ylim(15,49)


### convert ASFR to whole-population values (apply female ratio estimates)
ASFR_labels = ASFR$AGE
ASFR_breaks = c(10,14,19,24,29,34,39,44,49)
pop_together = pop_together %>%
  mutate(agegroup_ASFR = cut(age,breaks = ASFR_breaks, include.lowest = T, labels = ASFR_labels)) %>%
  ungroup() %>%
  group_by(agegroup_ASFR) %>%
  mutate(ASFR_group_percent = population/sum(population),
         interim = ASFR_group_percent * female_prop) 

ASFR_group_ratios = pop_together %>%
  group_by(agegroup_ASFR) %>%
  summarise(female_prop = sum(interim))

Pop_ASFR = ASFR %>%
  rename(agegroup_ASFR = AGE)%>% 
  left_join(ASFR_group_ratios) %>%
  mutate(ASFR = ASFR * female_prop) %>%
  select(agegroup_ASFR,ASFR) 
         

### adapt ASFR to model age groups         
pop_conversion = pop_setting_orig %>%
  mutate(agegroup_ASFR = cut(age,breaks = ASFR_breaks, include.lowest = T, labels = ASFR_labels),
         agegroup_MODEL = cut(age,breaks = age_groups_num, include.lowest = T, labels = age_group_labels)) %>%
  left_join(Pop_ASFR) %>%
  select(-agegroup_ASFR) %>%
  ungroup() %>% group_by(agegroup_MODEL) %>%
  mutate(agegroup_percent = population/sum(population),
         interim = agegroup_percent * ASFR) 
pop_conversion$interim[is.na(pop_conversion$interim)]=0
model_pregnancy_agegroups = aggregate(pop_conversion$interim, 
                              by=list(category= pop_conversion$agegroup_MODEL), FUN=sum)
colnames(model_pregnancy_agegroups) = c('age_group','prop')   


### save
#colnames: risk_group, age_group, prop, source
prevalence_pregnancy = model_pregnancy_agegroups %>% 
  mutate(risk_group = 'pregnant_women',
         source = 'DHS analysis + UN Pop prospects female ratio')

ggplot(data=prevalence_pregnancy) + 
  geom_point(aes(x=prop*100,y=age_group)) +
  # xlim(0,1) +
  xlab("Age-specific fertility ratio (%)") + 
  ylab("model age groups") + 
  labs(title="") 

save(prevalence_pregnancy, file = "1_inputs/prevalence_pregnancy.Rdata")