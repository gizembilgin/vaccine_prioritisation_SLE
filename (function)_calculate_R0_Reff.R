### This program contains a function to fit beta to the basic reproduction number (R0), and
### a function to calculate the effective reproduction number (Reff).
### These two functions are used to adjust transmission to circulating variants, and estimate
### the transmission advantage of immune-escape variants.

#### R0 #########################################################################
# R0 = p(NGM) - spectral radius, i.e., absolute maximum eigenvalue
# NGM = Next Generation Matrix = contact_matrix_adjusted * diag{beta*suscept*(gamma+lota(1-gamme)*1/delta)}

#(A) setup
beta_fitted_values = data.frame()
strain_list = c('WT','delta','omicron')
NGM_R0_list = list()

for (loop in 1:length(strain_list)){
  strain = strain_list[loop]
  if (strain == 'WT'){R0_to_fit = 2.79
  } else if (strain %in% c("delta","omicron")){R0_to_fit = 5.08} #ASSUMPTION: transmission advantage of Omicron entirely immune-escape!
  
  contact_matrix_adjust = matrix(data = 0, nrow = num_age_groups, ncol = num_age_groups)
  for (i in 1:num_age_groups){
    for (j in 1:num_age_groups){
      contact_matrix_adjust[i,j] = contact_matrix[i,j] * pop[i]/pop[j]
    }
  }
  
  #(B) ballpark beta
  contact_sum = rowSums(contact_matrix[,1:ncol(contact_matrix)])
  contact_ave=sum(contact_sum*(pop/sum(pop))) 
  beta_ball_park = R0_to_fit/(contact_ave*AverageSymptomaticPeriod)
  
  #(C) fit!
  minimise_this <- function(beta) {
    diag_matrix = beta*AverageSymptomaticPeriod*(gamma+(1-gamma)*lota)
    diag_matrix = diag(diag_matrix,num_age_groups)
    diag_matrix = suscept*diag_matrix
    
    NGM_R0 <- contact_matrix_adjust %*% diag_matrix
    R0_beta <- abs(eigen(NGM_R0)$values[1])
    
    fit = abs(R0_to_fit-R0_beta)
    
    return(fit)
  }
  
  #(D) check fit!
  beta_optimised = optimize(minimise_this,c(beta_ball_park*1/4,beta_ball_park*4))$minimum
  beta_optimised = rep(beta_optimised,num_age_groups)
  
  beta_check = beta_optimised[1]
  diag_matrix = beta_check*AverageSymptomaticPeriod*(gamma+(1-gamma)*lota)
  diag_matrix = diag(diag_matrix,num_age_groups)
  diag_matrix = suscept*diag_matrix
  NGM_R0 <- contact_matrix_adjust %*% diag_matrix
  R0_beta <- abs(eigen(NGM_R0)$values[1])
  
  if (! round(R0_beta,digits = 2) == round(R0_to_fit,digits = 2)){stop('beta fitting is not working!')}

  row = data.frame(beta_optimised[1],strain)
  colnames(row) = c('beta_optimised','strain')
  beta_fitted_values = rbind(beta_fitted_values,row)
  
  NGM_R0_list[[loop]] = NGM_R0
}
#_______________________________________________________________________________


#### Reff ######################################################################
Reff_time_step <- function(parameters,next_state){
  
  #NGM_modified = NGM_R0 * NPI * (1-pre-existing immunity)
  
  num = which(strain_list == strain_now)
  NGM_R0 = NGM_R0_list[[num]]
  
  #(A) NPI  
  NGM_modified = NGM_R0 * parameters$NPI  
    
  #(B) pre-existing immunity, i.e., vax and recovery
  if (round(sum(next_state$pop) - sum(pop),digits =0 ) == 0) {
    next_state_immunity = next_state %>%
      mutate(pop = 
               case_when(
                 class == 'R' ~ pop * (1-parameters$rho), #immunity from infection
                 TRUE ~ pop
                 ))
    
    next_state_immunity = next_state_immunity %>% 
      left_join(parameters$VE, by = c("age_group", "dose", "vaccine_type", "risk_group"))
    next_state_immunity$VE[is.na(next_state_immunity$VE) == TRUE] = 0
    next_state_immunity = next_state_immunity %>% mutate(pop = pop * (1-VE)) #immunity from vaccination
    next_state_immunity = sum(next_state_immunity$pop)/sum(pop)
    
    NGM_modified = NGM_modified * next_state_immunity
    
  } else{
    stop("issue with time_step population fluctating")
  }
  
    Reff  <- abs(eigen(NGM_modified)$values[1]) 
    return(Reff)
}