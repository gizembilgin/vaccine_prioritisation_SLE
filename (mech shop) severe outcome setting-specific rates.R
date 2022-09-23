### This (mech shop) creates setting-specific estimates of severe outcomes. 
# COMEBACK - could include uncertainty in both country-specific estimate and variant multiplier.
# NOTE - this includes the incidence of hospital admission NOT occupancy, the WHO defines severe disease as - 
# "a patient with severe acute respiratory illness (fever and at least one sign or symptom of respiratory disease), AND requiring hospitalization"
# Hence, hosp - severe = unmet need!

VOC = strain_inital

##### (1/7) Load population-level wild-type estimate of severe outcomes
severe_outcome_0 <- read.csv('1_inputs/severe_outcome_country_level.csv')
severe_outcome_0$percentage = severe_outcome_0$percentage/100 #make it between 0-1
severe_outcome_0 <- severe_outcome_0[severe_outcome_0$outcome %in% c('death','severe_disease','hosp') &
                                       severe_outcome_0$country == setting
                                     ,-c(1,5)] #dropping ICU and ICR as we won't use them, removing source and country column
#_______________________________________________________________________________


#####(2/7) Load variant-specific multipliers
workshop <- read.csv('1_inputs/severe_outcome_variant_multiplier.csv')

if (VOC != 'WT'){
  #<interlude for omicron>
  workshop2 <- read.csv('1_inputs/severe_outcome_variant_multiplier_complex.csv') #omicron vs delta
  omicron_basis = workshop[workshop$variant == 'delta',]
  omicron_basis$variant = 'omicron'
  omicron_basis$source = paste(omicron_basis$source,'/',workshop2$source)
  omicron_basis <- omicron_basis %>%
    mutate(multiplier = case_when(
      outcome == 'hosp' ~ multiplier*workshop2$multiplier[workshop2$outcome == 'hosp'],
      outcome %in% c('ICU','death') ~ multiplier*workshop2$multiplier[workshop2$outcome == 'hosp_long']))
  #ASSUMPTION: hosp_long proportional to ICU and death
  workshop = rbind(workshop,omicron_basis)
  #<fin>
  
  workshop = workshop[workshop$variant == VOC,c('outcome','multiplier')]
  #_______________________________________________________________________________
  
  #####(3/7) Calculating population-level variant-specific estimate of severe outcomes
  #could be made faster, but the assumptions were are making would be less obvious
  severe_outcome_1 <- severe_outcome_0 %>%
    mutate(percentage = case_when(
      outcome == 'death' ~ percentage * workshop$multiplier[workshop$outcome == 'death'],
      outcome == 'severe_disease' ~ percentage * workshop$multiplier[workshop$outcome == 'ICU'], #ASSUMPTION
      outcome == 'hosp' ~ percentage * workshop$multiplier[workshop$outcome == 'hosp']
    ),variant=VOC)
  rm (omicron_basis,workshop2)
  
} else if (VOC == 'WT'){
  severe_outcome_1 = severe_outcome_0 %>% mutate(variant = VOC)
}
#_______________________________________________________________________________



#####(4/7) Calculating age-specific estimates of severe outcomes
load(file = '1_inputs/severe_outcome_age_distribution.Rdata') #adjusted values from Qatar
workshop = age_dn_severe_outcomes
workshop = workshop[workshop$setting == setting,]

severe_outcome_2 <- severe_outcome_1 %>%  
  left_join(workshop) %>% mutate(percentage=percentage*RR)
severe_outcome_FINAL <- severe_outcome_2 %>%
  select(outcome,outcome_long,age_group,percentage) 

rm(severe_outcome_0,severe_outcome_2)
#_______________________________________________________________________________



#####(5/7) Calculating YLL from death
#requires average age in age-group and life_expectancy of this study setting
workshop <- pop_setting_orig %>%
  mutate(agegroup = cut(age,breaks = age_groups_num, include.lowest = T,labels = age_group_labels)) 

workshop_sum <- aggregate(workshop$population, by=list(category=workshop$agegroup), FUN=sum)
colnames(workshop_sum) <-c('agegroup','pop')

workshop <- workshop %>% 
  left_join(workshop_sum) %>% 
  mutate(weight=population/pop,
         age_weight = (age+0.5)*weight) %>%
  group_by(agegroup) %>%
  summarise(average_age = sum(age_weight))

lifeExpect <- read.csv('1_inputs/UN_life_expectancy_est.csv')
lifeExpect = lifeExpect[lifeExpect$setting == setting,]

closest_age = data.frame()
for (i in 1:nrow(workshop)){
  closest_age = rbind(closest_age,lifeExpect$age[which.min(abs(lifeExpect$age-workshop$average_age[i]))])
}
colnames(closest_age) = c('age')

workshop = cbind(workshop,closest_age)

workshop <- workshop %>%
  left_join(lifeExpect) %>%
  select(agegroup,average_age,age,life_expectancy)
#"The average number of remaining years of life expected by a hypothetical cohort of individuals alive at age x 
# who would be subject during the remaining of their lives to the mortality rates of a given period."
# https://population.un.org/wpp/Download/Standard/Mortality/

YLL_FINAL = workshop %>%
  select(agegroup,life_expectancy) %>%
  rename(age_group = agegroup, YLL = life_expectancy)

#apply discounting using continuous approach, as per larson et al.
if (discounting_rate >0){YLL_FINAL$life_expectancy = (1/discounting_rate)*(1-exp(-discounting_rate*YLL_FINAL$life_expectancy ))}


YLL_row = severe_outcome_FINAL %>%
  filter(outcome == 'death') %>%
  mutate(outcome = 'YLL',
         outcome_long = 'YLL per death in this age_group multiplied by death rate') %>%
  left_join(YLL_FINAL) %>%
  mutate(percentage = percentage*YLL)
YLL_row = YLL_row[,c(1:4)]

severe_outcome_FINAL = rbind(severe_outcome_FINAL,YLL_row)

severe_outcome_FINAL = severe_outcome_FINAL %>%
  mutate(outcome_VE = case_when(
    outcome %in% c('death','YLL') ~ 'death',
    outcome %in% c('hosp','severe_disease') ~ 'severe_disease'
  ))

ggplot() + 
  geom_point(data=severe_outcome_FINAL[severe_outcome_FINAL$outcome != 'YLL',],
             aes(x=factor(age_group,level=age_group_labels),
                                           y=percentage,color=as.factor(outcome)),na.rm=TRUE) +
  xlab('age group') +
  labs(color='outcome') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'))
#_______________________________________________________________________________

save(severe_outcome_FINAL, file = "1_inputs/severe_outcome_FINAL.Rdata")