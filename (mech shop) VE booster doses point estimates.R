### This (mech shop) creates point estimates of d = 3 (booster) heterologous combinations with Pfizer
# strains  = "omicron" (no booster doses till omicron in SLE)
# outcomes = "any_infection"       "death"               "severe_disease"      "symptomatic_disease"
# vaccine type =  "Pfizer"  (with previous primary doses in AZ, Sinovac or Pfizer)

### Creates: VE_booster_estimates


require(ggpubr); require(readr);require(ggplot2); require(tidyverse)

#####  Visualise estimates from IVAC living systematic review ##########4#########################################################################
VE_raw <- read.csv("1_inputs/VE_WHO_forest_plot.csv",header=TRUE)

VE_raw = VE_raw %>% 
  filter(dose == 3 & vaccine_type == 'Pfizer') %>% 
  select(strain, vaccine_type, primary_if_booster, outcome,VE,lower_est,upper_est) %>%
  mutate(primary_if_booster_long = case_when(
    primary_if_booster == "AstraZeneca" ~ "AstraZeneca Vaxzevria (ChAdOx1)",
    primary_if_booster == "Johnson & Johnson" ~ "Johnson & Johnson Janssen (Ad26.COV2.S)",
    primary_if_booster == "Moderna" ~ "Moderna Spikevax (mRNA-1273)",
    primary_if_booster == "Pfizer" ~ "Pfizer-BioNTech Comirnaty (BNT162b2)",
    primary_if_booster == "Sinopharm" ~ "Sinopharm BIBP vaccine",
    primary_if_booster == "Sinovac" ~ "Sinovac Biotech CoronaVac"         
  ))

this_strain = 'omicron'
to_plot = VE_raw[VE_raw$strain == this_strain,]
plot_list = list()
for (i in 1:length(unique(to_plot$outcome))){
  outcome = unique(to_plot$outcome)[i]
  plot_list [[i]] <- ggplot(data=to_plot[to_plot$outcome==outcome,]) + 
    geom_pointrange(aes(x=VE,y=primary_if_booster_long,xmin=lower_est,xmax=upper_est)) +
    xlim(0,100) +
    theme_bw() + 
    xlab("") + 
    ylab("") + 
    theme(text=element_text(size=10), 
          plot.title=element_text(size=12))+
    labs(title=paste("VE against ",outcome,sep=""))
}
ggarrange(plot_list[[1]],plot_list[[4]],plot_list[[2]],plot_list[[3]],
                       common.legend = TRUE,
                       legend="bottom"
)




#####  Average across estimates from IVAC living systematic review ##########4#########################################################################
VE_estimates = VE_raw %>% 
  select(strain, vaccine_type, primary_if_booster, outcome,VE) %>% 
  group_by(strain, vaccine_type, primary_if_booster, outcome) %>%
  summarise(VE = sum(VE)/n())

plot_list = list()
for (i in 1:length(unique(to_plot$outcome))){
  outcome = unique(to_plot$outcome)[i]
  plot_list [[i]] <- ggplot() + 
    geom_point(data = VE_raw[VE_raw$outcome == outcome,], aes(x=VE,y=primary_if_booster),shape=1) +
    geom_point(data = VE_estimates[VE_estimates$outcome == outcome,],aes(x=VE,y=primary_if_booster)) +
    xlim(0,100) +
    theme_bw() + 
    xlab("") + 
    ylab("") + 
    theme(text=element_text(size=10), 
          plot.title=element_text(size=12)) +
    labs(title=paste("VE against ",outcome,sep=""))
}
plot_group = ggarrange(plot_list[[1]],plot_list[[4]],plot_list[[2]],plot_list[[3]],
                       common.legend = TRUE,
                       legend="bottom"
)
annotate_figure(plot_group,top=text_grob('booster dose'))



#####  Impute missing values based on previous analysis ########################################################################################
##only 'any_infection' needs imputing
imputed_rows = VE_estimates %>% 
  filter(outcome == 'symptomatic_disease' & (!primary_if_booster %in% unique(VE_estimates$primary_if_booster[VE_estimates$outcome == 'any_infection']))) %>%
  mutate(VE = VE * infection_sympt_ratio$mean[infection_sympt_ratio$strain == 'omicron'],
         outcome = 'any_infection')
VE_estimates = rbind(VE_estimates,imputed_rows)
#_________________________________________________________________________________________________________________________________________



#####  Compare to dose 1 and 2 ########################################################################################
to_plot = VE_estimates_imputed[VE_estimates_imputed$strain == 'omicron' & VE_estimates_imputed$dose !=3,]
plot_list = list()
for (i in 1:length(unique(to_plot$outcome))){
  outcome = unique(to_plot$outcome)[i]
  plot_list [[i]] <- ggplot(data=to_plot[to_plot$outcome==outcome,]) + 
    geom_pointrange(aes(x=VE,y=vaccine_type_long,color=as.factor(dose),shape=source,xmin=lower_est,xmax=upper_est)) +
    geom_point(data = VE_estimates[VE_estimates$outcome == outcome,],aes(x=VE,y=primary_if_booster)) +
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
#CHECKED: all d3>d2
#_________________________________________________________________________________________________________________________________________



#####  Save point estimate for booster doses ################################################################################
VE_booster_estimates = VE_estimates %>%
  mutate(dose = 3,
  vaccine_type_long = case_when(
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

save(VE_booster_estimates,file = "1_inputs/VE_booster_estimates.Rdata")



