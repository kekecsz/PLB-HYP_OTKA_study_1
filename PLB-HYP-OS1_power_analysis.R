
# This is the R script for the power analysis of the PLB-HYP OTKA Study 1 project.

library(BayesFactor) # for ttestBF
library(MASS) # for mvrnorm
library(tidyverse) # for tidy coding


simulate_data <- function(mean_sham_as_true_embedded,
                          mean_sham_as_true_whitenoise,
                          mean_true_as_true_relaxation,
                          mean_true_as_true_confusion,
                          
                          mean_sham_as_control_embedded,
                          mean_sham_as_control_whitenoise,
                          mean_true_as_control_relaxation,
                          mean_true_as_control_confusion,
                          
                          sd,
                          
                          r_as_trues,
                          r_as_controls,
                          r_trues_with_controls,
                          
                          n_per_subgroup){
  
  Sigma = matrix(c(1, r_trues_with_controls, r_as_trues, r_trues_with_controls,
                   r_trues_with_controls, 1, r_trues_with_controls, r_as_controls,
                   r_as_trues, r_trues_with_controls, 1, r_trues_with_controls,
                   r_trues_with_controls, r_as_controls, r_trues_with_controls, 1), nrow = 4)
  
  my_data_pre = as.data.frame(
    mvrnorm(n = n_per_subgroup*4,
            mu = c(0, 0, 0, 0),
            Sigma = Sigma)
  )
  
  my_data_pre_realistic_sd = my_data_pre * sd
  
  sham_as_true = rep(c("embedded", "whitenoise"), each = n_per_subgroup*2)
  sham_as_control = rep(c("whitenoise", "embedded"), each = n_per_subgroup*2)
  true_as_true = rep(c("relaxation", "confusion"), each = n_per_subgroup*2)
  true_as_control = rep(c("confusion", "relaxation"), each = n_per_subgroup*2)
  
  
  my_data_pre2 = cbind(my_data_pre_realistic_sd, sham_as_true, sham_as_control, true_as_true, true_as_control)
  
  my_col_names = c("expectancy_sham_as_true", "expectancy_sham_as_control", "expectancy_true_as_true", "expectancy_true_as_control", "sham_as_true", "sham_as_control", "true_as_true", "true_as_control")
  
  names(my_data_pre2) = my_col_names
  
  my_data_pre2[my_data_pre2[,"sham_as_true"] == "embedded", "expectancy_sham_as_true"] = my_data_pre2[my_data_pre2[,"sham_as_true"] == "embedded", "expectancy_sham_as_true"] + mean_sham_as_true_embedded
  my_data_pre2[my_data_pre2[,"sham_as_true"] == "whitenoise", "expectancy_sham_as_true"] = my_data_pre2[my_data_pre2[,"sham_as_true"] == "whitenoise", "expectancy_sham_as_true"] + mean_sham_as_true_whitenoise
  
  my_data_pre2[my_data_pre2[,"sham_as_control"] == "embedded", "expectancy_sham_as_control"] = my_data_pre2[my_data_pre2[,"sham_as_control"] == "embedded", "expectancy_sham_as_control"] + mean_sham_as_control_embedded
  my_data_pre2[my_data_pre2[,"sham_as_control"] == "whitenoise", "expectancy_sham_as_control"] = my_data_pre2[my_data_pre2[,"sham_as_control"] == "whitenoise", "expectancy_sham_as_control"] + mean_sham_as_control_whitenoise
  
  my_data_pre2[my_data_pre2[,"true_as_true"] == "relaxation", "expectancy_true_as_true"] = my_data_pre2[my_data_pre2[,"true_as_true"] == "relaxation", "expectancy_true_as_true"] + mean_true_as_true_relaxation
  my_data_pre2[my_data_pre2[,"true_as_true"] == "confusion", "expectancy_true_as_true"] = my_data_pre2[my_data_pre2[,"true_as_true"] == "confusion", "expectancy_true_as_true"] + mean_true_as_true_confusion
  
  my_data_pre2[my_data_pre2[,"true_as_control"] == "relaxation", "expectancy_true_as_control"] = my_data_pre2[my_data_pre2[,"true_as_control"] == "relaxation", "expectancy_true_as_control"] + mean_true_as_control_relaxation
  my_data_pre2[my_data_pre2[,"true_as_control"] == "confusion", "expectancy_true_as_control"] = my_data_pre2[my_data_pre2[,"true_as_control"] == "confusion", "expectancy_true_as_control"] + mean_true_as_control_confusion
  
  return(my_data_pre2)
}


analysis_code_main_lmBF_factors_none = function(data, rscale){
  data_with_ID = data %>% 
    mutate(ID = factor(paste0("ID_", 1:nrow(data))))
  
  data_long = data_with_ID %>% 
    gather(key = "expectancy_type", value = "expectancy", expectancy_sham_as_true:expectancy_true_as_control) %>% 
    filter(expectancy_type != "expectancy_sham_as_control") %>% 
    filter(expectancy_type != "expectancy_true_as_control") %>% 
    mutate(expectancy_type = factor(expectancy_type),
           sham_as_true  = factor(sham_as_true ),
           true_as_true = factor(true_as_true))
  
  bf_full = lmBF(expectancy ~ expectancy_type
                 , whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf = 1/matrix(bf_full)
  
  return(bf)
}

analysis_code_main_lmBF_factors_sham_as_true_true_as_true = function(data, rscale){
  data_with_ID = data %>% 
    mutate(ID = factor(paste0("ID_", 1:nrow(data))))
  
  data_long = data_with_ID %>% 
    gather(key = "expectancy_type", value = "expectancy", expectancy_sham_as_true:expectancy_true_as_control) %>% 
    filter(expectancy_type != "expectancy_sham_as_control") %>% 
    filter(expectancy_type != "expectancy_true_as_control") %>% 
    mutate(expectancy_type = factor(expectancy_type),
           sham_as_true  = factor(sham_as_true ),
           true_as_true = factor(true_as_true))

  bf_full = lmBF(expectancy ~ expectancy_type + sham_as_true + 
                   true_as_true, whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_reduced = lmBF(expectancy ~ sham_as_true + 
                      true_as_true, whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_diff = bf_full / bf_reduced
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
}

analysis_code_main_lmBF_factors_sham_as_true = function(data, rscale){
  data_with_ID = data %>% 
    mutate(ID = factor(paste0("ID_", 1:nrow(data))))
  
  data_long = data_with_ID %>% 
    gather(key = "expectancy_type", value = "expectancy", expectancy_sham_as_true:expectancy_true_as_control) %>% 
    filter(expectancy_type != "expectancy_sham_as_control") %>% 
    filter(expectancy_type != "expectancy_true_as_control") %>% 
    mutate(expectancy_type = factor(expectancy_type),
           sham_as_true  = factor(sham_as_true ),
           true_as_true = factor(true_as_true))

  bf_full = lmBF(expectancy ~ expectancy_type + sham_as_true, 
                 whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_reduced = lmBF(expectancy ~ sham_as_true, 
                    whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_diff = bf_full / bf_reduced
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
}



analysis_code_main_lmBF_factors_true_as_true = function(data, rscale){
  data_with_ID = data %>% 
    mutate(ID = factor(paste0("ID_", 1:nrow(data))))
  
  data_long = data_with_ID %>% 
    gather(key = "expectancy_type", value = "expectancy", expectancy_sham_as_true:expectancy_true_as_control) %>% 
    filter(expectancy_type != "expectancy_sham_as_control") %>% 
    filter(expectancy_type != "expectancy_true_as_control") %>% 
    mutate(expectancy_type = factor(expectancy_type),
           sham_as_true  = factor(sham_as_true ),
           true_as_true = factor(true_as_true))
  
  ### If the expectancy of sham techniques are comparable AND also the expectancy of true techniques are comparable
  
  bf_full = lmBF(expectancy ~ expectancy_type + true_as_true, 
                 whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_reduced = lmBF(expectancy ~ true_as_true, 
                    whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data_long)
  
  bf_diff = bf_full / bf_reduced
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
}

analysis_code_comparability_of_sham_techniques = function(data, rscale){
  # are expectancy of the different sham techniques comparable
  
  bf_ttest = ttestBF(formula = expectancy_sham_as_true ~ sham_as_true, data = data, rscale = rscale)
  bf = 1/matrix(bf_ttest)
  
  return(bf)
  }

analysis_code_comparability_of_true_techniques = function(data, rscale){
  # are expectancy of the different sham techniques comparable
  
  bf_ttest = ttestBF(formula = expectancy_true_as_true ~ true_as_true, data = data, rscale = rscale)
  bf = 1/matrix(bf_ttest)
  
  return(bf)
  }

analysis_compiler <- function(data, rscale){
  sham_comp = analysis_code_comparability_of_sham_techniques(data, rscale)
  true_comp = analysis_code_comparability_of_true_techniques(data, rscale)
  
  if((sham_comp >= 3) & (true_comp  >= 3)){
    bf = analysis_code_main_lmBF_factors_none(data, rscale)
  } else if((sham_comp >= 3) & (true_comp  < 3)){
    bf = analysis_code_main_lmBF_factors_true_as_true(data, rscale)
  } else if((sham_comp < 3) & (true_comp  >= 3)){
    bf = analysis_code_main_lmBF_factors_sham_as_true(data, rscale)
  } else if((sham_comp < 3) & (true_comp  < 3)){
    bf = analysis_code_main_lmBF_factors_sham_as_true_true_as_true(data, rscale)
  }
  
  return(bf)
}




simulation_plus_analysis <- function(
  mean_sham_as_true_embedded,
  mean_sham_as_true_whitenoise,
  mean_true_as_true_relaxation,
  mean_true_as_true_confusion,
  
  mean_sham_as_control_embedded,
  mean_sham_as_control_whitenoise,
  mean_true_as_control_relaxation,
  mean_true_as_control_confusion,
  
  sd,
  
  r_as_trues,
  r_as_controls,
  r_trues_with_controls,
  
  n_per_subgroup,
  
  rscale
){
  data = 
    simulate_data(
      mean_sham_as_true_embedded,
      mean_sham_as_true_whitenoise,
      mean_true_as_true_relaxation,
      mean_true_as_true_confusion,
      
      mean_sham_as_control_embedded,
      mean_sham_as_control_whitenoise,
      mean_true_as_control_relaxation,
      mean_true_as_control_confusion,
      
      sd,
      
      r_as_trues,
      r_as_controls,
      r_trues_with_controls,
      
      n_per_subgroup
    )
  
  bf = analysis_compiler(data = data, rscale = rscale)
  
 return(bf)
}











###########################################
#                                         #  
#            Power analysis               #
#                                         #
###########################################





iterations = 100

bfs = replicate(n = iterations, simulation_plus_analysis(
  mean_sham_as_true_embedded = 4.8,
  mean_sham_as_true_whitenoise = 4.8,
  mean_true_as_true_relaxation = 4.8,
  mean_true_as_true_confusion = 4.8,
  
  mean_sham_as_control_embedded = 2.3,
  mean_sham_as_control_whitenoise = 2.3,
  mean_true_as_control_relaxation = 2.3,
  mean_true_as_control_confusion = 2.3,
  
  sd = 2.35,
  
  r_as_trues = 0.63,
  r_as_controls = 0.63,
  r_trues_with_controls = 0,
  
  n_per_subgroup = 13,
  
  rscale = sqrt(2)/2 # sqrt(2)/2 or 1
))


sum(bfs >= 3)/length(bfs)
sum(bfs <= 1/3)/length(bfs)
