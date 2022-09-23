### This program runs 'CommandDeck' multiple times to understand the impact of current program targets:
### (1) during an outbreak, and (2) if rolled out prior to an outbreak.
### It compiles the results from these runs into panel A of Figure S4.2 in the Supplementary Material.


### (1) Overarching trackers #####################################################################################################
warehouse_table = data.frame() 
warehouse_plot = data.frame()
baseline_to_compare = "no further vaccine rollout"



### (2) Queue strategies to run ##################################################################################################
queue = list()

queue[[1]] = list(vax_strategy_description = "no further vaccine rollout",
                  vax_strategy_toggle = "off",
                  outbreak_timing = "during")

queue[[2]] = list(vax_strategy_description = 'current roll-out DURING outbreak',
                  vax_strategy_toggle = "on",
                  outbreak_timing = "during",
                  vax_strategy_toggles = vax_strategy_toggles_CURRENT_TARGET)  #roll out vaccine DURING outbreak

queue[[3]] = list(vax_strategy_description = 'current roll-out PRIOR to outbreak',
                  vax_strategy_toggle = "on",
                  outbreak_timing = "after",
                  vax_strategy_toggles = vax_strategy_toggles_CURRENT_TARGET)  #roll out vaccine PRIOR TO outbreak



### (3) Run  ##################################################################################################
for (ticket in 1:length(queue)){
  commands = queue[[ticket]]
  vax_strategy_description = commands$vax_strategy_description
  vax_strategy_toggle      = commands$vax_strategy_toggle
  outbreak_timing          = commands$outbreak_timing
  if ('vax_strategy_toggles' %in% names(commands)){
    vax_strategy_toggles = commands$vax_strategy_toggles
  }
  
  source(paste(getwd(),"/CommandDeck.R",sep=""))
  
  severe_outcome_projections = severe_outcome_log %>% 
    mutate(label = vax_strategy_description, 
           day = as.numeric(date - max(seed_date) )) # day = start of outbreak
  warehouse_plot = rbind(warehouse_plot,severe_outcome_projections)
  
  row = row %>% mutate(scenario = vax_strategy_description) %>% relocate(scenario, .before = colnames(row)[[1]])
  warehouse_table = rbind(warehouse_table,row)
}
#____________________________________________________________________________________________________________________________________



### (4) Save outputs  ##################################################################################################
results_warehouse_entry = list()
results_warehouse_entry[[1]] = warehouse_table
results_warehouse_entry[[2]] = warehouse_plot

#(A/B) Plot
warehouse_plot = warehouse_plot %>% 
  mutate(time = day) %>%
  filter(time>=0)

abs_plot_list = list()
for (i in 1:length(unique(warehouse_plot$outcome))){
  outcome = unique(warehouse_plot$outcome)[i]
  abs_plot_list [[i]] <- ggplot(data=warehouse_plot[warehouse_plot$outcome==outcome,]) + 
    geom_line(aes(x=time,y=proj,color=as.factor(label))) +
    labs(title=paste(outcome)) +
    labs(colour = "") +
    theme_bw() + 
    xlab("") + 
    ylab("")}

cum_plot_list = list()
for (i in 1:length(unique(warehouse_plot$outcome))){
  outcome = unique(warehouse_plot$outcome)[i]
  cum_plot_list [[i]] <- ggplot(data=warehouse_plot[warehouse_plot$outcome==outcome,]) + 
    geom_line(aes(x=time,y=proj_cum,color=as.factor(label))) +
    labs(title=paste(outcome)) +
    labs(colour = "") +
    theme_bw() + 
    xlab("") + 
    ylab("")
}


# 1 = death, 2 = hosp, 3 = severe_disease, 4 = YLL, 5 = cases
plot = ggarrange(abs_plot_list[[5]],cum_plot_list[[5]],
                 abs_plot_list[[2]],cum_plot_list[[2]],
                 abs_plot_list[[3]], cum_plot_list[[3]],
                 abs_plot_list[[1]],  cum_plot_list[[1]],
                 common.legend = TRUE,
                 legend="bottom",
                 ncol = 2,
                 nrow = 4)
results_warehouse_entry[[3]]= plot

#(B/B) Table
averted_table = warehouse_table[warehouse_table$scenario != baseline_to_compare,]
averted_table_rel = averted_table
for (i in 1:(length(queue)-1)){
  end = (length(unique(warehouse_plot$outcome))+1)
  averted_table[i,c(2:end)] = 
    warehouse_table[warehouse_table$scenario == baseline_to_compare,c(2:end)] -
    averted_table[i,c(2:end)] 
  
  averted_table_rel[i,c(2:end)] = 100 * averted_table[i,c(2:end)]/
    warehouse_table[warehouse_table$scenario == baseline_to_compare,c(2:end)]
}
table_list = list(absolute = averted_table, 
                  relative = averted_table_rel) #COMEBACK could be merged

results_warehouse_entry[[4]]= table_list

#____________________________________________________________________________________________________________________________________

results_warehouse[[receipt]] = results_warehouse_entry

