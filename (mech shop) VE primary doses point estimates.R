###This (mech shop) creates point estimates of VE by dose, vaccine_type, outcome and strain
# strains  = "delta"   "omicron"
# outcomes = "any_infection"       "death"               "severe_disease"      "symptomatic_disease"
# vaccine type = "AstraZeneca"       "Johnson & Johnson" "Moderna"           "Pfizer"            "Sinopharm"         "Sinovac"    

### Creates: VE_estimates_imputed

require(ggpubr); require(readr);require(ggplot2); require(tidyverse)

##### (1/2) Inital estimates from IVAC living systematic review ##########4#########################################################################
VE_estimates <- read.csv("1_inputs/VE_WHO_forest_plot.csv",header=TRUE)
colnames(VE_estimates)

VE_estimates = VE_estimates  %>% select(strain, vaccine_type, dose, outcome,VE,lower_est,upper_est)

this_strain = 'delta'
this_strain = 'omicron'
to_plot = VE_estimates[VE_estimates$strain == this_strain,]
plot_list = list()
for (i in 1:length(unique(to_plot$outcome))){
  outcome = unique(to_plot$outcome)[i]
  plot_list [[i]] <- ggplot(data=to_plot[to_plot$outcome==outcome,]) + 
    #geom_point(aes(x=VE,y=vaccine_type, color=as.factor(dose),shape=strain)) +
    geom_pointrange(aes(x=VE,y=vaccine_type,color=as.factor(dose),shape=strain,xmin=lower_est,xmax=upper_est)) +
    xlim(0,100) +
    xlab("") + 
    ylab("") + 
    labs(title=paste("VE against ",outcome,sep=""))
}
plot_group = ggarrange(plot_list[[1]],plot_list[[4]],plot_list[[2]],plot_list[[3]],
          common.legend = TRUE,
          legend="bottom"
)
annotate_figure(plot_group,top=text_grob(paste(this_strain)))

#NOTES ON ODD BEHAVIOUR: Sinopharm -> J&J 3rd is higher than Sinopharm 1 & 2 but lower than J&J 2nd
#repeat Sinopharm 3rd with J&J as J&J booster
workshop = VE_estimates %>% 
  filter(vaccine_type == 'Johnson & Johnson' & dose == 3 & outcome %in% c('severe_disease','death'))
workshop$dose = 2
VE_estimates = rbind(VE_estimates,workshop)


###Calculate ratios:
#(A/D) compare VE against death (where avaliable) to VE against severe disease
death = VE_estimates[VE_estimates$outcome == 'death',] %>% 
  select(strain,vaccine_type,dose,VE) %>%
  rename(death = VE)
severe_disease = VE_estimates[VE_estimates$outcome == 'severe_disease',] %>% 
  select(strain,vaccine_type,dose,VE) %>%
  rename(severe_disease = VE)

workshop = death %>% 
  left_join(severe_disease, by = c("strain", "vaccine_type", "dose")) %>%
  mutate(ratio = severe_disease/death)
workshop %>% summarise(mean = mean(ratio,na.rm=TRUE),
                       sd = sd(ratio,na.rm=TRUE))
# mean        sd
# 1 0.9667011 0.1131279
workshop %>% group_by(dose) %>%
  summarise(mean = mean(ratio,na.rm=TRUE),
                       sd = sd(ratio,na.rm=TRUE))
# dose  mean     sd
#   1     1 0.929 0.145 
#   2     2 1.02  0.0171
#close enough to one!
#____________________

#(B/D) compare VE against any infection (where avaliable) to symptomatic disease
symptomatic_disease = VE_estimates[VE_estimates$outcome == 'symptomatic_disease',] %>% 
  select(strain,vaccine_type,dose,VE) %>%
  rename(symptomatic_disease = VE)
any_infection = VE_estimates[VE_estimates$outcome == 'any_infection',] %>% 
  select(strain,vaccine_type,dose,VE) %>%
  rename(any_infection = VE)

workshop = symptomatic_disease %>% 
  left_join(any_infection, by = c("strain", "vaccine_type", "dose")) %>%
  mutate(ratio = any_infection/symptomatic_disease)
workshop %>% summarise(mean = mean(ratio,na.rm=TRUE),
                       sd = sd(ratio,na.rm=TRUE))
# mean        sd
# 1 0.9017756 0.2476893
workshop %>% group_by(dose) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# dose count    mean      sd
# 1     9   1.11   0.0834
# 2     5   0.779  0.231 
workshop %>% group_by(strain) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# strain  count  mean     sd
# 1 delta       5 1.03  0.109 
# 2 omicron    13 0.530 0.0450
infection_sympt_ratio = workshop %>% group_by(strain) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
#_________________________


#(C/D) compare VE against omicron (where avaliable) to delta
delta = VE_estimates[VE_estimates$strain == 'delta',] %>% 
  select(outcome,vaccine_type,dose,VE) %>%
  rename(delta = VE)
omicron = VE_estimates[VE_estimates$strain == 'omicron',] %>% 
  select(outcome,vaccine_type,dose,VE) %>%
  rename(omicron = VE)

workshop = delta %>% 
  left_join(omicron, by = c("outcome", "vaccine_type", "dose")) %>%
  mutate(ratio = omicron/delta)
workshop %>% summarise(mean = mean(ratio,na.rm=TRUE),
                       sd = sd(ratio,na.rm=TRUE))
# mean        sd
# 1 0.6210314 0.1543049
workshop %>% group_by(dose) %>%
  summarise(mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# dose  mean    sd
# 1     1 0.540 0.171
# 2     2 0.651 0.148
workshop %>% group_by(outcome) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# outcome             count    mean      sd
# 1 any_infection           9   0.447  0.0225
# 2 death                  11 NaN     NA     
# 3 severe_disease          6   0.620  0.164 
# 4 symptomatic_disease     7   0.710  0.112 
delta_omicron_ratio = workshop %>% group_by(outcome) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
delta_omicron_ratio$mean[delta_omicron_ratio$outcome == 'death'] = delta_omicron_ratio$mean[delta_omicron_ratio$outcome == 'severe_disease']
#since death is missing
#_________________________


#(D/D) dose one to two
dose_one = VE_estimates[VE_estimates$dose == 1,] %>% 
  select(strain,outcome,vaccine_type,VE) %>%
  rename(dose_one = VE)
dose_two = VE_estimates[VE_estimates$dose == 2,] %>% 
  select(strain,outcome,vaccine_type,VE) %>%
  rename(dose_two = VE)

workshop = dose_one %>% 
  left_join(dose_two, by = c("strain", "outcome", "vaccine_type")) %>%
  mutate(ratio = dose_one/dose_two)
# mean        sd
# 1 0.7445608 0.1569164
workshop %>% summarise(mean = mean(ratio,na.rm=TRUE),
                       sd = sd(ratio,na.rm=TRUE))

workshop %>% group_by(strain) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# strain  count  mean    sd
# 1 delta       6 0.755 0.139
# 2 omicron    16 0.709 0.232

workshop %>% group_by(vaccine_type) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# vaccine_type count  mean     sd
# 1 AstraZeneca      4 0.783 0.122 
# 2 Johnson & Johnson4 0.755 0.258 
# 2 Moderna          5 0.798 0.0949
# 3 Pfizer           2 0.735 0.203 
# 4 Sinopharm        5 0.642 0.164 
# 5 Sinovac          6 0.769 0.220 

workshop %>% group_by(outcome) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
# outcome             count  mean     sd
# 1 any_infection           6 0.698 0.146 
# 2 death                   7 0.885 0.0500
# 3 severe_disease          4 0.830 0.129 
# 4 symptomatic_disease     5 0.596 0.106 

dose_ratio = workshop %>% group_by(outcome) %>%
  summarise(count = sum(is.na(ratio)),
            mean = mean(ratio,na.rm=TRUE),
            sd = sd(ratio,na.rm=TRUE))
#_________________________
#_________________________________________________________________________________________________________________________________________




##### (2/2) Impute missing values based on previous analysis ###################################################################################
workshop = VE_estimates %>%
  mutate(source = case_when(
    is.na(VE) == FALSE ~ 'literature',
    TRUE ~ 'imputed'
  ))

#Step One: estimate dose one from dose two
for (s in unique(workshop$strain)){
  for (t in unique(workshop$vaccine_type)){
    for (o in unique(workshop$outcome)){
      workshop_rows = workshop[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == o, ]
      if (workshop_rows$source[workshop_rows$dose == 1] == "imputed" & 
          workshop_rows$source[workshop_rows$dose == 2] == "literature"){
        estimate = workshop_rows$VE[workshop_rows$dose == 2] * dose_ratio$mean[dose_ratio$outcome == o]
        workshop$VE[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == o & workshop$dose == 1] = estimate
        workshop$source_extend[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == o & workshop$dose == 1] = "dose one estimated from dose two"
      }
    }
  }
}
sum(workshop$source_extend == "dose one estimated from dose two",na.rm=TRUE)
#n=8
  
#Step Two: estimate severe_disease <-> death, and acquisition <-> symptomatic
#(A/B) severe_disease <-> death
#close enough to one, thus set equal
for (s in unique(workshop$strain)){
  for (t in unique(workshop$vaccine_type)){  
    for (d in c(1,2)){
        workshop_rows = workshop[workshop$strain == s & workshop$vaccine_type == t & workshop$dose == d & workshop$outcome %in% c('death','severe_disease'), ]
        #severe_outcome -> death
        if (workshop_rows$source[workshop_rows$outcome == 'death'] == "imputed" &
            is.na(workshop_rows$source_extend[workshop_rows$outcome == 'death']) &
            is.na(workshop_rows$VE[workshop_rows$outcome == 'severe_disease']) == FALSE){ #NOTE ASSUMPTION - may be imputing from already imputed value!!!
          estimate = workshop_rows$VE[workshop_rows$outcome == 'severe_disease'] 
          workshop$VE[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'death' & workshop$dose == d] = estimate
          workshop$source_extend[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'death' & workshop$dose == d] = "death inferred from severe_disease"
        }
        #death -> severe_outcome
        if (workshop_rows$source[workshop_rows$outcome == 'severe_disease'] == "imputed" &
            is.na(workshop_rows$source_extend[workshop_rows$outcome == 'severe_disease']) &
            is.na(workshop_rows$VE[workshop_rows$outcome == 'death']) == FALSE){ #NOTE ASSUMPTION - may be imputing from already imputed value!!!
          estimate = workshop_rows$VE[workshop_rows$outcome == 'death'] 
          workshop$VE[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'severe_disease' & workshop$dose == d] = estimate
          workshop$source_extend[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'severe_disease' & workshop$dose == d] = "severe_disease inferred from death"
        }
    }
  }
}
sum(workshop$source_extend == "death inferred from severe_disease",na.rm=TRUE) # n = 12
sum(workshop$source_extend == "severe_disease inferred from death",na.rm=TRUE) # n = 0


#(B/B) acquisition <-> symptomatic
for (s in unique(workshop$strain)){
    for (t in unique(workshop$vaccine_type)){  
      for (d in c(1,2)){
        workshop_rows = workshop[workshop$strain == s & workshop$vaccine_type == t & workshop$dose == d & workshop$outcome %in% c('any_infection','symptomatic_disease'), ]
        #any infection -> symptomatic infection
        if (workshop_rows$source[workshop_rows$outcome == 'any_infection'] == "imputed" &
            is.na(workshop_rows$source_extend[workshop_rows$outcome == 'any_infection']) &
            is.na(workshop_rows$VE[workshop_rows$outcome == 'symptomatic_disease']) == FALSE){ #NOTE ASSUMPTION - may be imputing from already imputed value!!!
          estimate = workshop_rows$VE[workshop_rows$outcome == 'symptomatic_disease'] * infection_sympt_ratio$mean[infection_sympt_ratio$strain == s]
          workshop$VE[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'any_infection' & workshop$dose == d] = estimate
          workshop$source_extend[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'any_infection' & workshop$dose == d] = "any_infection inferred from symptomatic_disease"
        }
        #symptomatic infection -> any infection
        if (workshop_rows$source[workshop_rows$outcome == 'symptomatic_disease'] == "imputed" &
            is.na(workshop_rows$source_extend[workshop_rows$outcome == 'symptomatic_disease']) &
            is.na(workshop_rows$VE[workshop_rows$outcome == 'any_infection']) == FALSE){ #NOTE ASSUMPTION - may be imputing from already imputed value!!!
          estimate = workshop_rows$VE[workshop_rows$outcome == 'any_infection'] * 1/infection_sympt_ratio$mean[infection_sympt_ratio$strain == s]
          workshop$VE[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'symptomatic_disease' & workshop$dose == d] = estimate
          workshop$source_extend[workshop$strain == s & workshop$vaccine_type == t & workshop$outcome == 'symptomatic_disease' & workshop$dose == d] = "symptomatic_disease inferred from any_infection"
        }  
    }
  }
}
sum(workshop$source_extend == "any_infection inferred from symptomatic_disease",na.rm=TRUE) # n = 4
sum(workshop$source_extend == "symptomatic_disease inferred from any_infection",na.rm=TRUE) # n = 8


#Step Three: estimate omicron from delta
for (o in unique(workshop$outcome)){
  for (t in unique(workshop$vaccine_type)){
      for (d in c(1,2)){
        workshop_rows = workshop[workshop$dose == d & workshop$vaccine_type == t & workshop$outcome == o, ]
        if (workshop_rows$source[workshop_rows$strain == "omicron"] == "imputed" & 
            is.na(workshop_rows$source_extend[workshop_rows$strain == 'omicron']) &
            is.na(workshop_rows$VE[workshop_rows$strain == 'delta']) == FALSE){
          estimate = workshop_rows$VE[workshop_rows$strain == 'delta'] * delta_omicron_ratio$mean[delta_omicron_ratio$outcome == o]
          workshop$VE[workshop$dose == d & workshop$vaccine_type == t & workshop$outcome == o & workshop$strain == "omicron"] = estimate
          workshop$source_extend[workshop$dose == d & workshop$vaccine_type == t & workshop$outcome == o & workshop$strain == "omicron"] = "omicron estimated from delta"
        }
      }
  }
}
sum(workshop$source_extend =="omicron estimated from delta",na.rm=TRUE) # n=12
if (nrow(workshop[is.na(workshop$VE),])>0){stop('Some values to go!')}


VE_estimates_imputed = workshop %>%
  mutate(vaccine_type_long = case_when(
    vaccine_type == "AstraZeneca" ~ "AstraZeneca Vaxzevria (ChAdOx1)",
    vaccine_type == "Johnson & Johnson" ~ "Johnson & Johnson Janssen (Ad26.COV2.S)",
    vaccine_type == "Moderna" ~ "Moderna Spikevax (mRNA-1273)",
    vaccine_type == "Pfizer" ~ "Pfizer-BioNTech Comirnaty (BNT162b2)",
    vaccine_type == "Sinopharm" ~ "Sinopharm BIBP vaccine",
    vaccine_type == "Sinovac" ~ "Sinovac Biotech CoronaVac"         
  ),
  vaccine_mode = case_when(
      vaccine_type == 'Pfizer' ~ 'mRNA',
      vaccine_type == 'Moderna' ~ 'mRNA',
      vaccine_type == 'AstraZeneca' ~ 'viral',
      vaccine_type == 'Sinopharm' ~ 'viral',
      vaccine_type == 'Sinovac' ~ 'viral',
      vaccine_type == 'Johnson & Johnson' ~ 'viral'
    ),
  outcome_family = case_when(
    outcome %in% c('any_infection','symptomatic_disease') ~ 'acquisition',
    outcome %in% c('severe_disease','death') ~ 'severe_outcome'
    
  ))


to_plot = VE_estimates_imputed %>%
  filter(strain == 'omicron' & dose !=3)
#to_plot = VE_estimates_imputed %>% filter(vaccine_type %in% c('Moderna','Pfizer','AstraZeneca'), strain == 'omicron')

plot_list = list()
for (i in 1:length(unique(to_plot$outcome))){
  outcome = unique(to_plot$outcome)[i]
  plot_list [[i]] <- ggplot(data=to_plot[to_plot$outcome==outcome,]) + 
    geom_pointrange(aes(x=VE,y=vaccine_type_long,color=as.factor(dose),shape=source,xmin=lower_est,xmax=upper_est)) +
    xlim(0,100) +
    xlab("") +
    theme_bw() + 
    scale_shape_manual(values=c(1,19)) +
    ylab("") + 
    labs(title=paste("VE against ",outcome,sep="")) + 
    theme(text=element_text(size=10), 
          plot.title=element_text(size=12))
}
plot_VE_point_estimates = ggarrange(plot_list[[1]],plot_list[[4]],plot_list[[3]],plot_list[[2]],
                                    common.legend = TRUE,
                                    legend="bottom")
plot_VE_point_estimates

save(VE_estimates_imputed,file = "1_inputs/VE_estimates_imputed.Rdata")


