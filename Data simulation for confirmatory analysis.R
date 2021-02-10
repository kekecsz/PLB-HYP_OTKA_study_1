library(BayesFactor)
library(MASS)
library(tidyverse) 

simulate_data <- function(mean_sham_described_hypnosis_embedded,
                          mean_sham_described_hypnosis_whitenoise,
                          mean_true_described_hypnosis_relaxation,
                          mean_true_described_hypnosis_confusion,
                          
                          mean_sham_described_control_embedded,
                          mean_sham_described_control_whitenoise,
                          mean_true_described_control_relaxation,
                          mean_true_described_control_confusion,
                          
                          sd,
                          
                          r_described_hypnosiss,
                          r_described_controls,
                          r_hypnosiss_with_controls,
                          
                          total_N){
  
  Sigma = matrix(c(1, r_hypnosiss_with_controls, r_described_hypnosiss, r_hypnosiss_with_controls,
                   r_hypnosiss_with_controls, 1, r_hypnosiss_with_controls, r_described_controls,
                   r_described_hypnosiss, r_hypnosiss_with_controls, 1, r_hypnosiss_with_controls,
                   r_hypnosiss_with_controls, r_described_controls, r_hypnosiss_with_controls, 1), nrow = 4)
  
  my_data_pre = as.data.frame(
    mvrnorm(n = total_N,
            mu = c(0, 0, 0, 0),
            Sigma = Sigma)
  )
  
  my_data_pre2 = as.vector(as.matrix(my_data_pre))
  expectancy = my_data_pre2 * sd
  
  ID = rep(factor(paste0("ID_", 1:total_N)), 4)
  
  trial_type = rep(c("sham", "true"), each = total_N*2)
  procedure_type = c(rep(c("embedded", "whitenoise"), total_N/2), rep(c("whitenoise", "embedded"), total_N/2), rep(c("relaxation", "confusion"), total_N/2), rep(c("confusion", "relaxation"), total_N/2))
  description_type = c(rep(c("control", "hypnosis"), each = total_N), rep(c("control", "hypnosis"), each = total_N))
  
  my_data = data.frame(ID = ID, expectancy = expectancy, trial_type = trial_type, procedure_type = procedure_type, description_type = description_type)
  
  my_data[trial_type == "sham" & description_type == "hypnosis" & procedure_type == "embedded", "expectancy"] = my_data[trial_type == "sham" & description_type == "hypnosis" & procedure_type == "embedded", "expectancy"] + mean_sham_described_hypnosis_embedded
  my_data[trial_type == "sham" & description_type == "hypnosis" & procedure_type == "whitenoise", "expectancy"] = my_data[trial_type == "sham" & description_type == "hypnosis" & procedure_type == "whitenoise", "expectancy"] + mean_sham_described_hypnosis_whitenoise
  
  my_data[trial_type == "sham" & description_type == "control" & procedure_type == "embedded", "expectancy"] = my_data[trial_type == "sham" & description_type == "control" & procedure_type == "embedded", "expectancy"] + mean_sham_described_control_embedded
  my_data[trial_type == "sham" & description_type == "control" & procedure_type == "whitenoise", "expectancy"] = my_data[trial_type == "sham" & description_type == "control" & procedure_type == "whitenoise", "expectancy"] + mean_sham_described_control_whitenoise
  
  
  my_data[trial_type == "true" & description_type == "hypnosis" & procedure_type == "relaxation", "expectancy"] = my_data[trial_type == "true" & description_type == "hypnosis" & procedure_type == "relaxation", "expectancy"] + mean_true_described_hypnosis_relaxation
  my_data[trial_type == "true" & description_type == "hypnosis" & procedure_type == "confusion", "expectancy"] = my_data[trial_type == "true" & description_type == "hypnosis" & procedure_type == "confusion", "expectancy"] + mean_true_described_hypnosis_confusion
  
  my_data[trial_type == "true" & description_type == "control" & procedure_type == "relaxation", "expectancy"] = my_data[trial_type == "true" & description_type == "control" & procedure_type == "relaxation", "expectancy"] + mean_true_described_control_relaxation
  my_data[trial_type == "true" & description_type == "control" & procedure_type == "confusion", "expectancy"] = my_data[trial_type == "true" & description_type == "control" & procedure_type == "confusion", "expectancy"] + mean_true_described_control_confusion
  
  my_data = my_data[order(c(1:total_N, 1:total_N, 1:total_N, 1:total_N)),]
  
  return(my_data)
}

my_data = simulate_data(
  mean_sham_described_hypnosis_embedded = 4.8,
  mean_sham_described_hypnosis_whitenoise = 4.8,
  mean_true_described_hypnosis_relaxation = 4.8,
  mean_true_described_hypnosis_confusion = 4.8,
  mean_sham_described_control_embedded = 2.3,
  mean_sham_described_control_whitenoise = 2.3,
  mean_true_described_control_relaxation = 2.3,
  mean_true_described_control_confusion = 2.3,
  sd = 2.35,
  r_described_hypnosiss = 0.63,
  r_described_controls = 0.63,
  r_hypnosiss_with_controls = 0,
  total_N = 13*4
)
my_data