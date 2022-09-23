# This program runs the 'CommandDeck' multiple times with varying prioritisation strategies including children 5 to 17.
# It compiles the results from these runs in Table 2 of the main paper.


### (1) Overarching trackers #####################################################################################################
warehouse_table = data.frame() 
warehouse_plot = data.frame()
vax_strategy_toggle = "on"
age_split_results = "Y"


### (2) Planning #################################################################################################################
#calculating poss level of coverage
sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose == 1]) #0.318 current coverage
sum(pop[3:num_age_groups])/sum(pop) #0.54 max for SLE with 18+ eligibility
sum(pop[2:num_age_groups])/sum(pop)  #0.86 max for SLE with 5+ eligibility


### (3) Queue strategies to run ##################################################################################################
queue = list()

#Section (1/3) run baseline - current vaccination targets
queue[[1]] = list(vax_strategy_description = 'current vaccination targets (51.6%)',
                  vax_strategy_toggles = vax_strategy_toggles_CURRENT_TARGET) 


#Section (2/3) - current eligibility (~50%) -> then expand to kids 60,70,80%
target_list = list(0.6,0.7,0.8)

for (i in 1:length(target_list)){
  target = target_list[[i]]
  target_percentage = target * 100
  workshop_doses = target - sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose == 1])/100
  workshop_doses = round(workshop_doses * sum(pop))
  this_vax_strategy = vax_strategy_toggles_CURRENT_TARGET
  this_vax_strategy$vax_strategy_num_doses = as.integer(workshop_doses)
  this_vax_strategy$vax_age_strategy = "adults_then_children"
  
  queue[[(1+i)]] = list(vax_strategy_description = paste('current then expand to children ',target_percentage,'%',sep=''),
                        vax_strategy_toggles = this_vax_strategy) 
}


#Section (3/3) - expand now to include kids 
target_list = list(0.516,0.6,0.7,0.8)

for (i in 1:length(target_list)){
  target = target_list[[i]]
  target_percentage = target * 100
  workshop_doses = target - sum(vaccination_history_POP$coverage_this_date[vaccination_history_POP$date == max(vaccination_history_POP$date) & vaccination_history_POP$dose == 1])/100
  workshop_doses = round(workshop_doses * sum(pop))
  this_vax_strategy = vax_strategy_toggles_CURRENT_TARGET
  this_vax_strategy$vax_strategy_num_doses = as.integer(workshop_doses)
  this_vax_strategy$vax_age_strategy = "uniform"
  
  queue[[(4+i)]] = list(vax_strategy_description = paste('expand to children now ',target_percentage,'%',sep=''),
                        vax_strategy_toggles = this_vax_strategy) 
}



### (4) Run  ##################################################################################################
for (ticket in 1:length(queue)){
  
  commands = queue[[ticket]]
  
  vax_strategy_description = commands$vax_strategy_description
  vax_strategy_toggles = commands$vax_strategy_toggles
  
  source(paste(getwd(),"/CommandDeck.R",sep=""))
  
  severe_outcome_projections = severe_outcome_log %>% 
    mutate(label = vax_strategy_description, day = as.numeric(date - date_start ))
  warehouse_plot = rbind(warehouse_plot,severe_outcome_projections)
  
  row = row %>% mutate(scenario = vax_strategy_description) %>% relocate(scenario, .before = colnames(row)[[1]])
  warehouse_table = rbind(warehouse_table,row)
}

age_split_results = "N"
#____________________________________________________________________________________________________________________________________



### (5) Save outputs  ##################################################################################################
results_warehouse_entry = list()
results_warehouse_entry[[1]] = warehouse_table
results_warehouse_entry[[2]] = warehouse_plot

scenario_table_list = list()
scenario_plot_list = list()

#(A/B) absolute outcome plot #################
if (outbreak_timing == "after"){
  warehouse_plot = warehouse_plot %>% mutate(time = day)
} else if (outbreak_timing %in% c("during","off")){
  warehouse_plot = warehouse_plot %>% mutate(time = date)
}

section_1 = c("current vaccination targets (51.6%)", "current then expand to children 60%", "current then expand to children 70%", "current then expand to children 80%")
section_2 = c("current vaccination targets (51.6%)","expand to children now 51.6%","expand to children now 60%","expand to children now 70%","expand to children now 80%")
section_list = list(section_1,section_2)

for (section in 1:length(section_list)){
  list_plot_commands = section_list[[section]]
  workshop = warehouse_plot[warehouse_plot$label %in% list_plot_commands, ]
  
  #let's still plot with net outcomes (not differentiated by children and adults)
  workshop = workshop %>% 
    group_by(date,outcome,label,day,time) %>%
    summarise(proj = sum(proj),
              proj_cum = sum(proj_cum),.groups = "keep")
  
  plot_list = list()
  for (i in 1:length(unique(workshop$outcome))){
    outcome = unique(workshop$outcome)[i]
    plot_list [[i]] <- ggplot(data=workshop[workshop$outcome==outcome,]) + 
      geom_point(aes(x=time,y=proj,color=as.factor(label))) +
      labs(title=paste(outcome)) +
      theme_bw() + 
      xlab("") + 
      ylab("")}
  # 1 = death, 2 = hosp, 3 = severe_disease, 4 = YLL, 5 = cases
  plot = ggarrange(plot_list[[5]], plot_list[[1]], 
                   common.legend = TRUE,
                   legend="bottom")
  annotate_figure(plot, top = text_grob('absolute outcome by scenario', face = 'bold', size = 16))
  scenario_plot_list[[section]] = plot
}

results_warehouse_entry[[3]]= scenario_plot_list
#________________________________________________


#(B/B) cumulative outcome table ############
baseline_to_compare = "current vaccination targets (51.6%)"

table2 = warehouse_table
table2 = table2 %>% 
  pivot_longer(
    cols = 3:ncol(table2) ,
    names_to = 'outcome',
    values_to = 'num'
  ) 

#add overall
workshop = table2 %>%
  group_by(scenario,outcome) %>%
  summarise(num=sum(num),.groups = "keep") %>%
  mutate(macro_age_group = 'overall')
table2 = rbind(table2,workshop)

#compare to baseline
baseline = table2 %>% 
  filter(scenario == baseline_to_compare) %>%
  rename(baseline_num=num) %>%
  select(-scenario)

table2 = table2 %>%
  left_join(baseline,by=c('macro_age_group','outcome'))  %>%
  mutate(abs_reduction = num - baseline_num,
         rel_reduction = 100*(num - baseline_num)/baseline_num)

#reorder
table2$outcome = factor(table2$outcome,levels=c('cases','severe_disease','hosp','death','YLL'))
table2$macro_age_group = factor(table2$macro_age_group,levels=c('children <5','children 5-17','adults','overall'))
table2$scenario = factor(table2$scenario, levels = 
                           c(baseline_to_compare,
                             "expand to children now 51.6%","expand to children now 60%",
                             "expand to children now 70%","expand to children now 80%",
                             "current then expand to children 60%","current then expand to children 70%", 
                             "current then expand to children 80%")   )
table2 = table2 %>% arrange(scenario,macro_age_group,outcome)

#print
options(scipen = 1000)
print = table2 %>%
  filter(! scenario %in% c(baseline_to_compare))  %>%
  mutate(abs_reduction = round(abs_reduction),
         rel_reduction = round(rel_reduction,digits=1),
         together_outcome = paste(outcome,macro_age_group,sep='_')) %>%
  ungroup() %>%
  select(-num,-baseline_num,-macro_age_group,-outcome) 

print_num = print %>% mutate(value = abs_reduction,type = 'absolute') %>%
  pivot_wider(
    id_cols = scenario,
    names_from = together_outcome,
    values_from = value)
print_abs = print %>% mutate(value = paste(rel_reduction,"%",sep=''),type = 'relative') %>%
  pivot_wider(
    id_cols = scenario,
    names_from = together_outcome,
    values_from = value)

print = rbind(print_num,print_abs) %>% arrange(scenario)

time = Sys.time()
time = gsub(':','-',time)
write.csv(print,file=paste(rootpath,'x_results/table2',vax_strategy_toggles_CURRENT_TARGET$vax_strategy_vaccine_type,time,'.csv',sep=''))
results_warehouse_entry[[4]]= print
results_warehouse[[receipt]] = results_warehouse_entry
#____________________________________________________________________________________________________________________________________



