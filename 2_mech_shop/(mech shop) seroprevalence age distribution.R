### This (mech shop) takes existing seroprevalence studies and estimates the seroprevalence for our model age groups.
### This work needs to be re-run every time the age groups in the model change.
### Seroprevalence is used ONLY when fitting the model.
### Creates: seroprev.
#_________________________________________________________________________________________________________________________________________________________________________

seroprev_raw = read.csv("1_inputs/seroprevalence_RAW.csv",header=TRUE)
setting = "SLE"


#(A) select setting specific seroprevalence and split by overall/age-specific
SLE = seroprev_raw[seroprev_raw$setting == setting,]
SLE_overall = SLE[SLE$age_group == "overall",]
SLE_age = SLE[!(SLE$age_group == "overall"),]
names(SLE_age)[names(SLE_age) == 'age_group'] <- 'agegroup_SERO' #change column name


#(B) examine underlying population size of these age groupings
SLE_pop <- pop_orig[pop_orig$country == "SLE",]
underlying_age_grouping <- c(0,9,19,39,59,110)
SLE_pop <- SLE_pop %>%
  mutate(agegroup_SERO = cut(age,breaks = underlying_age_grouping, include.lowest = T, labels = SLE_age$agegroup_SERO),
         agegroup_MODEL = cut(age,breaks = age_groups_num, include.lowest = T, labels = age_group_labels)) %>%
  ungroup() %>%
  group_by(agegroup_MODEL) %>%
  mutate(model_group_percent = population/sum(population))
pop_SERO = SLE_pop %>%
  group_by(agegroup_SERO) %>%
  summarise(pop = sum(population)) 
pop_SERO = pop_SERO%>% 
  mutate(percentage = pop/sum(pop_SERO$pop))


#(C) check that age-specific aligns with pop-level
check = SLE_age %>% left_join(pop_SERO) %>% select(agegroup_SERO,seroprev,percentage) %>%
  mutate(interim = seroprev*percentage)
sum(check$interim) # = 2.22 NOT reported 2.6 -> adjust!


#(D) adjust age-specific to align with population estimate
factor = SLE_overall$seroprev/sum(check$interim)
factor # = 1.169
SLE_age$seroprev = SLE_age$seroprev * factor #CHECKED: re-runing (C) gives 2.6 :)


#(E) translate to model age groups
workshop = SLE_pop %>% left_join(SLE_age) %>% 
  mutate(interim = seroprev*model_group_percent) %>%
  group_by(agegroup_MODEL) %>%
  summarise(seroprev = sum(interim))

#CHECK
pop_MODEL = SLE_pop %>%
  group_by(agegroup_MODEL) %>%
  summarise(pop = sum(population))
check = workshop %>% left_join(pop_MODEL) %>% mutate(interim = seroprev * pop /sum(pop_MODEL$pop))
sum(check$interim) #=2.6

#SAVE
seroprev = workshop %>% 
  rename(age_group = agegroup_MODEL) %>%
  mutate(year = 2021, setting = "SLE") %>%
  select(year,setting,age_group,seroprev)

save(seroprev, file = "1_inputs/seroprev.Rdata")

rm(check,workshop,seroprev_raw,factor,
   pop_MODEL,pop_SERO)




