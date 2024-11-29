
####################################
###    Set working directory     ###
####################################
figure_location = "C:\\Users\\User\\Documents\\"

####################################
### Set parameters for analysis  ###
####################################
# scale parameter for the Bayesian analysis
rscale_to_use = 1
# set permuation test number of iterations (10000 will run for a several minutes on an average work computer)
perm_iter = 10000

##############################
###  Sensitivity analysis  ###
##############################

# Set this to 1 to run the analysis only on the subsample of
# participants who were deceived by the deception
sensitivity_analysis = 0


#######################
###  Load packages  ###
#######################

library(data.table)
library(tidyverse)
library(lme4)
library(BayesFactor) # for ttestBF
library(r2glmm) # for r2beta
library(ggdist) # for rain cloud plot
library(ggsci) # for rain cloud plot color palette
library(pbapply) # for progress bar on iterative rerun of the confirmatory analysis
library(boot) # for bootstrapping	
library(ggpubr) # for ggarrange
library(grid)
library(gsheet)

###############################
###  Load custom functions  ###
###############################

# compute mode of a distribution

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# correlation, to allow using tidy data in correlation functions

r_to_boot = function(data_for_r, indices){
  d <- data_for_r[indices,] # allows boot to select sample 	
  r = cor(d, method = "spearman", use = "complete.obs")[1,2]
  return(r)
}

corr_boot <- function(data_for_r){
  # bootstrapping with 10000 replications 	
  results.boot <- boot(data=data_for_r, statistic=r_to_boot, 	
                       R=10000)	
  r = results.boot$t0
  
  boot_CI = boot.ci(results.boot, type="bca", index=1)
  CI_lb = unlist(boot_CI)$bca4
  CI_ub = unlist(boot_CI)$bca5
  
  return(data.frame(r = r, CI_lb = CI_lb, CI_ub = CI_ub))
}

m_to_boot = function(data_for_r, indices){
  d <- data_for_r[indices] # allows boot to select sample 	
  m = mean(d, na.rm = T)
  return(m)
}

mean_boot <- function(data_for_r){
  # bootstrapping with 10000 replications 	
  results.boot <- boot(data=data_for_r, statistic=m_to_boot, 	
                       R=10000)	
  m = results.boot$t0
  
  boot_CI = boot.ci(results.boot, type="bca", index=1)
  CI_lb = unlist(boot_CI)$bca4
  CI_ub = unlist(boot_CI)$bca5
  
  return(data.frame(m = m, CI_lb = CI_lb, CI_ub = CI_ub))
}

##################
### Read data  ###
##################
data_pre = read.csv("https://raw.githubusercontent.com/kekecsz/PLB-HYP_OTKA_study_1/refs/heads/main/Data/PLB_HYP_S1_data.csv")

#####################################
### Data management - eligibility ###
#####################################

### Exclude participants who are not eligible (these exclusions were pre-registered)

# Exclude participants with technical errors during the experiments

ids_to_exclude = c("GJTG20", 
                   "BAHJ57",
                   "ISL7",
                   "VKSP98",
                   "SZESZI2",
                   "NRNB27",
                   "MEMT98",
                   "TCGV11",
                   "KPSZM41",
                   "2e611d51")

data_pre_noexperror = data_pre[!(data_pre$unique_id %in% ids_to_exclude), ]

# no participants were below 18 years of age, so no exclusions in this step

data_pre_noexperror$age_range

# no participants had incomplete sessions

# no participant had tried hypnosis before

data_pre_noexperror$tried_hypnosis_before
data_pre_noexperror$tried_hypnosis_before_for_pain

# exclude participants who have attended a lecture on hypnosis before

data_pre_noexperror_nohypnosisstudents = data_pre_noexperror[
  which(!(data_pre_noexperror$knowledge_source %like% "Attended a lecture about hypnosis")), ]


############################################
### Data management - recoding variables ###
############################################

# recoding age

data_pre_recoded = data_pre_noexperror_nohypnosisstudents %>% 
  mutate(age = as.numeric(recode(age_range,
                      "18 - 24" = 21,
                      "25 - 34" = 30,
                      "18 - 24" = 21,
                      "25 - 34" = 30,
                      "35 - 44" = 40,
                      "45 - 54" = 50,
                      "55 - 64" = 60,
                      "65 - 74" = 70,
                      "75 - 84" = 80,
                      "85 or older" = 90)))

# designate factors

data_pre_recoded = data_pre_recoded %>% 
  mutate(
          description_type_1 = factor(description_type_1),
          description_type_2 = factor(description_type_2),
          description_type_3 = factor(description_type_3),
          description_type_4 = factor(description_type_4),
          felt_discomfort = factor(felt_discomfort),
          gender = factor(gender),
          language = factor(language),
          prefered_procedure = factor(prefered_procedure),
          procedure_type_1 = factor(procedure_type_1),
          procedure_type_2 = factor(procedure_type_2),
          procedure_type_3 = factor(procedure_type_3),
          procedure_type_4 = factor(procedure_type_4),
          suspected_control = factor(suspected_control),
          c = factor(suspected_real_hypnosis),
          trial_type_1 = factor(trial_type_1),
          trial_type_2 = factor(trial_type_2),
          trial_type_3 = factor(trial_type_3),
          trial_type_4 = factor(trial_type_4),
          which_real_control = factor(which_real_control),
          which_real_hypnosis = factor(which_real_hypnosis)
          )


# adding hypnotizability category

data_pre_recoded = data_pre_recoded %>% 
  mutate(hypnotizability_category = case_when(hypnotizability_total %in% 0:4 ~ "low",
                                            hypnotizability_total %in% 5:8 ~ "medium",
                                            hypnotizability_total %in% 9:12 ~ "high"))


names(data_pre_recoded) <- gsub(x = names(data_pre_recoded), pattern = "\\-", replacement = "_")  

data = data_pre_recoded


# Create new variables: deceived_by_plb and deceived_by_fakecontrol

# we create a new variable called deceived_by_plb. This will take the value "Yes" if the person reported that they
# did not suspect that one of the techniques described to them as hypnosis was actually control, or
# if they designated either "embedded" or "whitenoise" as the real hypnosis out of the techniques described as hypnosis.
# Otherwise, this variable takes the value "No".

data$deceived_by_plb = "No"

data = data %>% 
  mutate(deceived_by_plb = "No")
data = data %>% 
  mutate(deceived_by_plb = replace(deceived_by_plb, suspected_control == "No", "Yes"))
data = data %>% 
  mutate(deceived_by_plb = replace(deceived_by_plb, which_real_hypnosis == "embedded", "Yes"))
data = data %>% 
  mutate(deceived_by_plb = replace(deceived_by_plb, which_real_hypnosis == "whitenoise", "Yes"))

# we create a new variable called deceived_by_fakecontrol. This will take the value "Yes" if the person reported that they
# did not suspect that one of the techniques described to them as control was actually hypnosis, or
# if they designated either "embedded" or "whitenoise" as the real hypnosis out of the techniques described as control.
# Otherwise, this variable takes the value "No".


data = data %>% 
  mutate(deceived_by_fakecontrol = "No")
data = data %>% 
  mutate(deceived_by_fakecontrol = replace(deceived_by_fakecontrol, suspected_real_hypnosis == "No", "Yes"))
data = data %>% 
  mutate(deceived_by_fakecontrol = replace(deceived_by_fakecontrol, which_real_control == "embedded", "Yes"))
data = data %>% 
  mutate(deceived_by_fakecontrol = replace(deceived_by_fakecontrol, which_real_control == "whitenoise", "Yes"))


### simplify the ID variable

# verifying that all IDs are unique in the dataset

length(data$unique_id)
length(unique(data$unique_id))

# now we give a simpler ID code to participants

data$ID = paste0("ID_", 1:length(data$unique_id))
  



### Exclude outliers, compute difference scores, and standardize (z-transform) EEG variables of interest

## Exclude outliers and z-transform spectral power data
power_measures_of_interest = c("FZ_alpha",
                               "FZ_theta",
                               "FZ_gamma",
                               "OZ_alpha",
                               "OZ_theta",
                               "OZ_gamma")


for(i in 1:length(power_measures_of_interest)){
  for(j in 1:4){
    experience_var_withnum = paste0(power_measures_of_interest[i], "_z_trans_experience", j)
    cases_to_exclude = which(abs(as.matrix(data[,experience_var_withnum])-mean(as.matrix(data[,experience_var_withnum]))) > 3*sd(as.matrix(data[,experience_var_withnum]))) ### Removing observations more than 3SD away from the mean in the EEG feature
    data[cases_to_exclude,experience_var_withnum] = NA
    print(paste0("outlier excluded from", experience_var_withnum, ": ", sum(is.na(data[,experience_var_withnum]))))
    
    output_var_withnum = paste0(power_measures_of_interest[i], "_ebds_", j) ### ebds stands for experience-baseline difference scaled
    data = data.frame(data, newCol = NA)
    names(data)[length(names(data))] = output_var_withnum
    data[,output_var_withnum] = as.numeric(scale(data[,experience_var_withnum]))
  }
}

## exclude outliers from connectivity data

connectivity_measures_of_interest = c("beta1_connect_OZ_PO4",
                             "alpha2_connect_PZ_area_FZ_area",
                             "theta_connect_O1_PZ")

for(i in 1:length(connectivity_measures_of_interest)){
  baseline_var = paste0(connectivity_measures_of_interest[i], "_baseline1")
  cases_to_exclude = which(abs(as.matrix(data[,baseline_var])-mean(as.matrix(data[,baseline_var]))) > 3*sd(as.matrix(data[,baseline_var]))) ### Removing observations more than 3SD away from the mean in the EEG feature
  data[cases_to_exclude,baseline_var] = NA
  print(paste0("outlier excluded from", baseline_var, ": ", sum(is.na(data[,baseline_var]))))
  for(j in 1:4){
    experience_var_withnum = paste0(connectivity_measures_of_interest[i], "_experience", j)
    cases_to_exclude = which(abs(as.matrix(data[,experience_var_withnum])-mean(as.matrix(data[,experience_var_withnum]))) > 3*sd(as.matrix(data[,experience_var_withnum]))) ### Removing observations more than 3SD away from the mean in the EEG feature
    data[cases_to_exclude,experience_var_withnum] = NA
    print(paste0("outlier excluded from", experience_var_withnum, ": ", sum(is.na(data[,experience_var_withnum]))))
  }
}

### compute the difference of pre-trial "Baseline" and "Experience" phase 
### of connectivity features of interest
### We also Z-transform (scale) variables to reduce inter-individual noise, 
### to avoid multicollinearity, and to aid comparability and visualization across features

for(i in 1:length(connectivity_measures_of_interest)){
  baseline_var = paste0(connectivity_measures_of_interest[i], "_baseline1")
  for(j in 1:4){
    experience_var_withnum = paste0(connectivity_measures_of_interest[i], "_experience", j)
    output_var_withnum = paste0(connectivity_measures_of_interest[i], "_ebds_", j) ### ebds stands for experience-baseline difference scaled
    data = data.frame(data, newCol = NA)
    names(data)[length(names(data))] = output_var_withnum
    data[,output_var_withnum] = data[,experience_var_withnum]-data[,baseline_var]
    data[,output_var_withnum] = as.numeric(scale(data[,output_var_withnum]))
  } 
}

connectivity_measures_without_trialnum = paste0(connectivity_measures_of_interest, "_ebds")
power_measures_of_without_trialnum = paste0(power_measures_of_interest, "_ebds")

EEG_change_vars_without_trialnum = c(power_measures_of_without_trialnum, connectivity_measures_without_trialnum)

### convert data to long format


data_long = melt(setDT(data), 
                 measure = patterns("description_type", "expectancy", "hypnosis_depth", "procedure_type", "trial_type", EEG_change_vars_without_trialnum),
                 variable.name = 'trial_number', value.name = c("description_type", "expectancy", "hypnosis_depth", "procedure_type", "trial_type", EEG_change_vars_without_trialnum))

data_long = data_long %>% 
  arrange(ID)



### determine order of factor levels

data_long = data_long %>% 
  mutate(procedure_type = factor(procedure_type, levels = c("relaxation", "confusion", "embedded", "whitenoise")))



# create a deceived_by_deception variable, which takes the value of deceived_by_plb 
# if the description type is hypnosis, or it takes the value of deceived_by_fakecontrol
# if the description type is control

data_long = data_long %>% 
  mutate(deceived_by_deception = case_when(
    (description_type == "hypnosis" & deceived_by_plb == "Yes") ~ "Yes",
    (description_type == "control" & deceived_by_fakecontrol == "Yes") ~ "Yes",
    TRUE ~ "No"))


# trial number converted to a numeric variable

data_long$trial_number = as.numeric(as.character(data_long$trial_number))


### Create data subsets for easier modeling


data_long_desc_control = data_long %>% 
  filter(description_type == "control")

data_long_desc_hypnosis = data_long %>% 
  filter(description_type == "hypnosis")





###############################################
###############################################
#########        Results       ################
###############################################
###############################################

###################################
### Participant characteristics ###
###################################

# gender
table(data$gender)
table(data$gender)[1]/sum(table(data$gender))

# age range
table(data$age_range)
table(data$age_range)[1]/sum(table(data$age_range))
table(data$age_range)[2]/sum(table(data$age_range))
table(data$age_range)[3]/sum(table(data$age_range))
table(data$age_range)[4]/sum(table(data$age_range))

# knowledge level
data %>% 
  summarize(mean = mean(knowledge_level_on_hypnosis),
            sd = sd(knowledge_level_on_hypnosis))
min(data$knowledge_level_on_hypnosis)
max(data$knowledge_level_on_hypnosis)

data %>% 
  ggplot()+
  aes(x = knowledge_level_on_hypnosis) +
  geom_density(fill = "lightblue")

# hypnotizability
sum(is.na(data$hypnotizability_total))

data %>% 
  drop_na(hypnotizability_total) %>% 
  summarize(mean = mean(hypnotizability_total),
            sd = sd(hypnotizability_total),
                    min = min(hypnotizability_total),
                    max = max(hypnotizability_total))
data %>% 
  ggplot()+
  aes(x = hypnotizability_total)+
  geom_density(fill = "lightblue")

# attitude toward hypnosis
data %>% 
  summarize(mean = mean(attitude_towards_hypnosis),
            sd = sd(attitude_towards_hypnosis),
            min = min(attitude_towards_hypnosis),
            max = max(attitude_towards_hypnosis))
data %>% 
  ggplot()+
  aes(x = attitude_towards_hypnosis)+
  geom_density(fill = "lightblue")

# belief about hypnoanalgesia effectiveness
data %>% 
  summarize(mean = mean(effective_of_hypnoanalgesia),
            sd = sd(effective_of_hypnoanalgesia),
            min = min(effective_of_hypnoanalgesia),
            max = max(effective_of_hypnoanalgesia))
data %>% 
  ggplot()+
  aes(x = effective_of_hypnoanalgesia)+
  geom_density(fill = "lightblue")

###############################
### Deception effectiveness ###
###############################

### How many were deceived by the placebo technique

data %>% 
  count(suspected_control) %>% 
  mutate(prop.table(n))

data %>% 
  filter(suspected_control == "Yes") %>%
  count(deceived_by_plb) %>% 
  mutate(prop.table(n))

data %>% 
  count(deceived_by_plb) %>% 
  mutate(prop.table(n))

### How many were deceived by the placebo technique

data %>% 
  count(suspected_real_hypnosis) %>% 
  mutate(prop.table(n))

data %>% 
  filter(suspected_real_hypnosis == "Yes") %>%
  count(deceived_by_fakecontrol) %>% 
  mutate(prop.table(n))

data %>% 
  count(deceived_by_fakecontrol) %>% 
  mutate(prop.table(n))


### Which was the more convincing placebo in terms of deception effectiveness

## embedded hypnosis

data_long %>%
  filter(procedure_type == "embedded", description_type == "control") %>% 
  count(which_real_control) %>% 
  mutate(prop.table(n))

data_long %>%
  filter(procedure_type == "embedded", description_type == "hypnosis") %>% 
  count(which_real_hypnosis) %>% 
  mutate(prop.table(n))

14/46

## white noise hypnosis
 
data_long %>%
  filter(procedure_type == "whitenoise", description_type == "control") %>% 
  count(which_real_control) %>% 
  mutate(prop.table(n))
  
data_long %>%
  filter(procedure_type == "whitenoise", description_type == "hypnosis") %>% 
  count(which_real_hypnosis) %>% 
  mutate(prop.table(n))

18/46

data_long = data_long %>% 
  mutate(trial_type = factor(trial_type),
         description_type = factor(description_type))




#####################################
####    Sensitivity analysis?    ####
#####################################

if(sensitivity_analysis == 1){
  data_long = data_long %>% 
    filter(deceived_by_deception == "Yes")
}


###############################################################################
####    The effect of induction procedure and description on expectancy    ####
###############################################################################


### Descriptive statistics of expectancy evoked by the different procedures,

#===============#
##  Table 1    ##

## when described as control
data_long %>% 
  filter(description_type == "control") %>% 
  group_by(procedure_type) %>% 
  summarize(mean = mean(expectancy),
            sd = sd(expectancy),
            n = n())

## when described as hypnosis
data_long %>% 
  filter(description_type == "hypnosis") %>% 
  group_by(procedure_type) %>% 
  summarize(mean = mean(expectancy),
            sd = sd(expectancy),
            n = n())


## average within subject effect of description type (collapsed across procedures) on expectancy
data_long %>% 
  group_by(unique_id, description_type) %>% 
  summarize(mean = mean(expectancy)) %>% 
  pivot_wider(names_from = description_type, values_from = mean) %>% 
  mutate(difference = hypnosis - control) %>% 
  ungroup() %>% 
  summarize(mean = mean(difference, na.rm = T),
            sd = sd(difference, na.rm = T))

## confidence interval for the effect of description_type

mod_mancheck = lmer(expectancy ~ description_type + (1|ID), data = data_long)
summary(mod_mancheck)

confint(mod_mancheck)
# Conclusion: the manipulation was successful

## moderating effect of trial order tested.
mod_Q1_a2 = lmer(expectancy ~ trial_number * description_type + (1|ID), data = data_long)

summary(mod_Q1_a2)
confint(mod_Q1_a2)


## average effect of procedure type (collapsed across description types)
data_long %>% 
  group_by(procedure_type) %>%
  summarize(mean = mean(expectancy),
            sd = sd(expectancy))


exp_output_table = data.frame(procedure_type = rep(unique(data_long$procedure_type), 2), 
                              description_type = rep(unique(data_long$description_type), each = 4),
                              mean = NA,
                              CI_lb = NA,
                              CI_ub = NA)
  

for(i in 1:length(unique(data_long$procedure_type))){
  for(j in 1:length(unique(data_long$description_type))){
    boot_results = data_long %>% 
      filter(procedure_type == unique(data_long$procedure_type)[i], 
             description_type == unique(data_long$description_type)[j]) %>% 
      select(expectancy) %>% 
      as.matrix() %>% 
      as.vector() %>% 
      mean_boot()
    
    exp_output_table[(exp_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (exp_output_table[, "description_type"] == unique(data_long$description_type)[j]), "mean"] = boot_results$m
    exp_output_table[(exp_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (exp_output_table[, "description_type"] == unique(data_long$description_type)[j]), "CI_lb"] = boot_results$CI_lb
    exp_output_table[(exp_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (exp_output_table[, "description_type"] == unique(data_long$description_type)[j]), "CI_ub"] = boot_results$CI_ub
    
      }
}


procedure_type_labels <- c("relaxation", "confusion", "embedded", "white noise")
names(procedure_type_labels) <- c("relaxation", "confusion", "embedded", "whitenoise")

library(grid)

grob_a <- grobTree(textGrob("Relaxation", x=0.35,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_1_a = data_long %>%
  filter(procedure_type == "relaxation") %>% 
  ggplot() +
  aes(x = description_type, y = expectancy, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = exp_output_table %>% 
                  filter(procedure_type == "relaxation"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = exp_output_table %>% 
               filter(procedure_type == "relaxation"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Expectancy",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_a) +
  scale_fill_aaas()

grob_b <- grobTree(textGrob("Confusion", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_1_b = data_long %>%
  filter(procedure_type == "confusion") %>% 
  ggplot() +
  aes(x = description_type, y = expectancy, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = exp_output_table %>% 
                  filter(procedure_type == "confusion"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = exp_output_table %>% 
               filter(procedure_type == "confusion"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Expectancy",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_b) +
  scale_fill_aaas()

grob_c <- grobTree(textGrob("Embedded", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_1_c = data_long %>%
  filter(procedure_type == "embedded") %>% 
  ggplot() +
  aes(x = description_type, y = expectancy, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = exp_output_table %>% 
                  filter(procedure_type == "embedded"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = exp_output_table %>% 
               filter(procedure_type == "embedded"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Expectancy",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_c) +
  scale_fill_aaas()

grob_d <- grobTree(textGrob("White noise", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_1_d = data_long %>%
  filter(procedure_type == "whitenoise") %>% 
  ggplot() +
  aes(x = description_type, y = expectancy, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = exp_output_table %>% 
                  filter(procedure_type == "whitenoise"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = exp_output_table %>% 
               filter(procedure_type == "whitenoise"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Expectancy",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_d) +
  scale_fill_aaas()

figure_1 <- ggarrange(figure_1_a, figure_1_b, figure_1_c, figure_1_d,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, 
                      common.legend = TRUE, legend="bottom")
figure_1





### create high resolution image of figure
#jpeg(paste0(figure_location, "figure_1.jpg"), width = 15, height = 15, units = 'cm', res = 300)
#figure_1 # Make plot
#dev.off()

# Rain cloud plot of the effect of procedure type and description type on expectancy

description_type_labels <- c("When described as control", "When described as hypnosis")
names(description_type_labels) <- c("control", "hypnosis")

figure_1_old = data_long %>% 
  ggplot() +
  aes(x = procedure_type, y = expectancy, fill = procedure_type) +
  
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # outliers should not show up for the boxplot
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    # ploting on left side
    side = "left",
    # adjusting position
    justification = 1.1,
    # adjust grouping (binning) of observations
    binwidth = 0.12
  ) +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Expectancy",
    fill = "Procedure type"
  ) +
  scale_y_continuous(breaks=0:10) +
  facet_grid(description_type ~ ., labeller = labeller(description_type = description_type_labels)) +
  scale_fill_aaas() +
  coord_flip()

suppressWarnings(print(figure_1_old))




###############################
####    Hypothesis test    ####
###############################

# We will use statistical inference to test the hypotheses. The confirmatory analysis will be conducted on the sample of participants fitting the eligibility criteria when data collection is successfully finished according to the stopping rules. 
# In the confirmatory analysis we will restrict our analysis to the trials with description type: hypnosis.
# We will contrast the likelihood of observing the data we observed under the following two models:
# Model 0 (M0) assumes that the sham hypnosis technique described as hypnosis evokes comparable expected hypnosis depth to the true hypnosis technique described as hypnosis. 
# Model 1 (M1) assumes that the sham hypnosis technique described as hypnosis evokes different expected hypnosis depth than the true hypnosis technique described as hypnosis.

## Variables

# Variables in the analysis:
# description_type: hypnosis vs. control procedure
# procedure_type: embedded, whitenoise, relaxation, or confusion
# trial_type: conventional hypnosis vs. sham hypnosis
# expectancy: an expectation by the recipient of the treatment that the treatment will elicit the desired response, it is a 11 point numerical rating scale ranging from 0 (Not Hypnotized at all) to 10 (Extremely Hypnotized)


## Comparability of procedure types within the same trial type

# We will start by exploring whether the different procedure types evoke comparable expectancy within each trial type. (That is, testing if there is evidence supporting the comparability of the expectancy evoked by the two sham hypnosis techniques, and the two true hypnosis techniques.)
# We will build two Bayesian mixed effect linear regression models. In the “full model” we will predict expectancy with trial type (sham vs. true) and procedure type (embedded, whitenoise, relaxation, confusion) as predictors, and with the random intercept of participant ID. 
# The reduced model will be the same with the exception that procedure type will not be included in the model. We will contrast the two models using Bayes factor to see whether there is evidence that the two models are comparable.

# Bayes factor (M0 vs M1) of the comparability of procedure type

analysis_code_comparability_of_procedure_types = function(data, rscale){
  
  data = data %>%
    filter(description_type == "hypnosis")
  
  bf_full = lmBF(expectancy ~ trial_type + procedure_type + ID, 
                 whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data) # in the “full model" we will predict expectancy with trial type, and procedure type as predictors, and the random intercept of participant ID
  
  bf_reduced = lmBF(expectancy ~ trial_type + ID, whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data) # in the reduced model the procedure type is excluded
  
  bf_diff = bf_full / bf_reduced # 
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
  
}

analysis_code_comparability_of_procedure_types(data= data_long, rscale=1)

## Testing the main hypothesis

# If in the above mentioned comparability analsysis we find evidence supporting the comparability of the full model with the reduced model, indicating that the procedure types are comparable with each other within the same trial type, we will use the reduced model to make statistical inference on the main hypothesis. Otherwise, we will use the full model for statistical inference.
# For this reason we have two functions below to perform the main statistical test, one using the full model (incuding procedure type), and one using the reduced model (without procedure type).

### main statistical test function - if procedures types are comparable within the same trial type

analysis_code_main_lmBF_comparable_procedure_type = function(data, rscale){
  data = data %>%
    filter(description_type == "hypnosis")
  
  bf_mod1 = lmBF(expectancy ~ trial_type + ID
                 , whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data)
  
  bf_mod2 = lmBF(expectancy ~ ID
                 , whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data)
  
  bf_diff = bf_mod1 / bf_mod2
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
}

### main statistical test function - if procedures types are not comparable within the same trial type

analysis_code_main_lmBF_not_comparable_procedure_type = function(data, rscale){
  
  data = data %>%
    filter(description_type == "hypnosis")
  bf_mod1 = lmBF(expectancy ~ trial_type + procedure_type + ID
                 , whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data)
  
  bf_mod2 = lmBF(expectancy ~ procedure_type + ID
                 , whichRandom="ID", rscaleFixed = rscale, rscaleCont = rscale, data = data)
  bf_diff = bf_mod1 / bf_mod2
  
  bf = 1/matrix(bf_diff)
  
  return(bf)
}

### Running the test

# This function performs the comparability analysis (shown above), followed by the appropriate main statistical test chosen based on the result of the comparability analysis. 

analysis_compiler <- function(data, rscale){
  comparability_of_procedure_types = analysis_code_comparability_of_procedure_types(data, rscale)
  
  if(comparability_of_procedure_types >= 3){
    bf = analysis_code_main_lmBF_comparable_procedure_type(data, rscale)
  } else if(comparability_of_procedure_types < 3){
    bf = analysis_code_main_lmBF_not_comparable_procedure_type(data, rscale)
  }
  
  return(bf)
}

# The function computes the Bayes Factor (M0 vs. M1) corresponding to the effect of trial type (sham vs. true) in the regression model.
# In case the Bayes factor is 3 or more, analysis support the Model 0.
# In case the Bayes factor us under 0.33 then the analysis supports Model 1.
# In case Bayes factor is between 3 and 0.33 then the analysis returned inconclusive report.

analysis_compiler(data=data_long, rscale= rscale_to_use)

# Rerunning the analysis 1000 times
# since the results were very close to the inference threshold, and due to
# randomness in the analysis, sometimes were above, sometimes below the threshold,
# we rerun the analysis 1000 times to see the patterm of results.

# This code runs for about 30 minutes.
# results_of_confirmatory_analysis = pbreplicate(1000, analysis_compiler(data=data_long, rscale= rscale_to_use))

# proportion of BFs that reach or exceed BF = 3
# sum(results_of_confirmatory_analysis >= 3)/length(results_of_confirmatory_analysis)

# plotting the results
# enframe(results_of_confirmatory_analysis) %>% 
#   ggplot() +
#   aes(x = value) +
#   geom_density() +
#   xlim(c(2.5, 3.5)) +
#   geom_vline(xintercept = 3, linetype="dashed")


### A better test is to compare each placebo procedure separately to each conventional procedure. 

data_long_embedded_vs_relaxation = data_long %>%
  filter(description_type == "hypnosis") %>% 
  filter(procedure_type == "embedded" | procedure_type == "relaxation")



bf_mod1_embedded_vs_relaxation = lmBF(expectancy ~ trial_type + ID
               , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_embedded_vs_relaxation)

bf_mod2_embedded_vs_relaxation = lmBF(expectancy ~ ID
               , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_embedded_vs_relaxation)

bf_diff_embedded_vs_relaxation = bf_mod1_embedded_vs_relaxation / bf_mod2_embedded_vs_relaxation

bf_embedded_vs_relaxation = 1/matrix(bf_diff_embedded_vs_relaxation)

bf_embedded_vs_relaxation


data_long_embedded_vs_confusion = data_long %>%
  filter(description_type == "hypnosis") %>% 
  filter(procedure_type == "embedded" | procedure_type == "confusion")

bf_mod1_embedded_vs_confusion = lmBF(expectancy ~ trial_type + ID
                                      , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_embedded_vs_confusion)

bf_mod2_embedded_vs_confusion = lmBF(expectancy ~ ID
                                      , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_embedded_vs_confusion)

bf_diff_embedded_vs_confusion = bf_mod1_embedded_vs_confusion / bf_mod2_embedded_vs_confusion

bf_embedded_vs_confusion = 1/matrix(bf_diff_embedded_vs_confusion)

bf_embedded_vs_confusion


data_long_whitenoise_vs_relaxation = data_long %>%
  filter(description_type == "hypnosis") %>% 
  filter(procedure_type == "whitenoise" | procedure_type == "relaxation")



bf_mod1_whitenoise_vs_relaxation = lmBF(expectancy ~ trial_type + ID
                                      , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_whitenoise_vs_relaxation)

bf_mod2_whitenoise_vs_relaxation = lmBF(expectancy ~ ID
                                      , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_whitenoise_vs_relaxation)

bf_diff_whitenoise_vs_relaxation = bf_mod1_whitenoise_vs_relaxation / bf_mod2_whitenoise_vs_relaxation

bf_whitenoise_vs_relaxation = 1/matrix(bf_diff_whitenoise_vs_relaxation)

bf_whitenoise_vs_relaxation


data_long_whitenoise_vs_confusion = data_long %>%
  filter(description_type == "hypnosis") %>% 
  filter(procedure_type == "whitenoise" | procedure_type == "confusion")

bf_mod1_whitenoise_vs_confusion = lmBF(expectancy ~ trial_type + ID
                                     , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_whitenoise_vs_confusion)

bf_mod2_whitenoise_vs_confusion = lmBF(expectancy ~ ID
                                     , whichRandom="ID", rscaleFixed = rscale_to_use, rscaleCont = rscale_to_use, data = data_long_whitenoise_vs_confusion)

bf_diff_whitenoise_vs_confusion = bf_mod1_whitenoise_vs_confusion / bf_mod2_whitenoise_vs_confusion

bf_whitenoise_vs_confusion = 1/matrix(bf_diff_whitenoise_vs_confusion)

bf_whitenoise_vs_confusion




##############################################################################################
####    The effect of induction procedure and description on subjective hypnosis depth    ####
##############################################################################################

### Descriptive statistics of hypnosis depth evoked by the different procedures,
## when described as control

data_long %>% 
  filter(description_type == "control") %>% 
  group_by(procedure_type) %>% 
  summarize(mean = mean(hypnosis_depth),
            sd = sd(hypnosis_depth),
            n = n())

## when described as hypnosis

data_long %>% 
  filter(description_type == "hypnosis") %>% 
  group_by(procedure_type) %>% 
  summarize(mean = mean(hypnosis_depth),
            sd = sd(hypnosis_depth),
            n = n())

## average within subject effect of description type (collapsed across procedures)
data_long %>% 
  group_by(unique_id, description_type) %>% 
  summarize(mean = mean(hypnosis_depth)) %>% 
  pivot_wider(names_from = description_type, values_from = mean) %>% 
  mutate(difference = hypnosis - control) %>% 
  ungroup() %>% 
  summarize(mean = mean(difference),
            sd = sd(difference))

## confidence interval for the effect of description_type

mod_description_type_on_hypnosis_depth = lmer(hypnosis_depth ~ description_type + (1|ID), data = data_long)
summary(mod_description_type_on_hypnosis_depth)
confint(mod_description_type_on_hypnosis_depth)


## moderating effect of trial order tested.
mod_Q1_a3 = lmer(hypnosis_depth ~ trial_number * description_type + (1|ID), data = data_long)

summary(mod_Q1_a3)
confint(mod_Q1_a3)


## average effect of procedure type (collapsed across description types)
data_long %>% 
  group_by(procedure_type) %>%
  summarize(mean = mean(hypnosis_depth),
            sd = sd(hypnosis_depth))


mod_description_type_and_procedure_type_on_hypnosis_depth = lmer(hypnosis_depth ~ description_type * procedure_type + (1|ID), data = data_long)
summary(mod_description_type_and_procedure_type_on_hypnosis_depth)
confint(mod_description_type_and_procedure_type_on_hypnosis_depth)



# Rain cloud plot of the effect of procedure type and description type on hypnosis_depth

depth_output_table = data.frame(procedure_type = rep(unique(data_long$procedure_type), 2), 
                              description_type = rep(unique(data_long$description_type), each = 4),
                              mean = NA,
                              CI_lb = NA,
                              CI_ub = NA)


for(i in 1:length(unique(data_long$procedure_type))){
  for(j in 1:length(unique(data_long$description_type))){
    boot_results = data_long %>% 
      filter(procedure_type == unique(data_long$procedure_type)[i], 
             description_type == unique(data_long$description_type)[j]) %>% 
      select(hypnosis_depth) %>% 
      as.matrix() %>% 
      as.vector() %>% 
      mean_boot()
    
    depth_output_table[(depth_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (depth_output_table[, "description_type"] == unique(data_long$description_type)[j]), "mean"] = boot_results$m
    depth_output_table[(depth_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (depth_output_table[, "description_type"] == unique(data_long$description_type)[j]), "CI_lb"] = boot_results$CI_lb
    depth_output_table[(depth_output_table[, "procedure_type"] == unique(data_long$procedure_type)[i]) & (depth_output_table[, "description_type"] == unique(data_long$description_type)[j]), "CI_ub"] = boot_results$CI_ub
    
  }
}



grob_a <- grobTree(textGrob("Relaxation", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_2_a = data_long %>%
  filter(procedure_type == "relaxation") %>% 
  ggplot() +
  aes(x = description_type, y = hypnosis_depth, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = depth_output_table %>% 
                  filter(procedure_type == "relaxation"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = depth_output_table %>% 
               filter(procedure_type == "relaxation"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Hypnosis depth",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_a) +
  scale_fill_aaas()

grob_b <- grobTree(textGrob("Confusion", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_2_b = data_long %>%
  filter(procedure_type == "confusion") %>% 
  ggplot() +
  aes(x = description_type, y = hypnosis_depth, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = depth_output_table %>% 
                  filter(procedure_type == "confusion"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = depth_output_table %>% 
               filter(procedure_type == "confusion"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Hypnosis depth",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_b) +
  scale_fill_aaas()

grob_c <- grobTree(textGrob("Embedded", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_2_c = data_long %>%
  filter(procedure_type == "embedded") %>% 
  ggplot() +
  aes(x = description_type, y = hypnosis_depth, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = depth_output_table %>% 
                  filter(procedure_type == "embedded"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = depth_output_table %>% 
               filter(procedure_type == "embedded"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Hypnosis depth",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_c) +
  scale_fill_aaas()

grob_d <- grobTree(textGrob("White noise", x=0.35,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=15, fontface="italic")))

figure_2_d = data_long %>%
  filter(procedure_type == "whitenoise") %>% 
  ggplot() +
  aes(x = description_type, y = hypnosis_depth, fill = description_type) +
  geom_violin() +
  geom_errorbar(data = depth_output_table %>% 
                  filter(procedure_type == "whitenoise"),
                aes(x = description_type, ymin = CI_lb, ymax = CI_ub, fill = description_type), width = 0.4, inherit.aes=F)+
  geom_point(data = depth_output_table %>% 
               filter(procedure_type == "whitenoise"), 
             aes(y = mean, x = description_type), shape = 21, size = 3, fill = "white") +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Procedure type",
    y = "Hypnosis depth",
    fill = "Label"
  ) +
  scale_y_continuous(breaks=0:10, limits = c(0, 10))+
  annotation_custom(grob_d) +
  scale_fill_aaas()

figure_2 <- ggarrange(figure_2_a, figure_2_b, figure_2_c, figure_2_d,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, 
                      common.legend = TRUE, legend="bottom")
figure_2





### create high resolution image of figure
#jpeg(paste0(figure_location, "figure_2.jpg"), width = 15, height = 15, units = 'cm', res = 300)
#figure_2 # Make plot
#dev.off()




### create high resolution image of figure
# jpeg(paste0(figure_location, "figure_2.jpg"), width = 15, height = 24, units = 'cm', res = 300)
# figure_2 # Make plot
# dev.off()


############################################################
####    The influence of expectancy on hypnosis depth   ####
############################################################
trial_type_labels <- c("unconventional", "conventional")
names(trial_type_labels) <- c("sham", "TRUE")

figure_3 = data_long %>% 
  ggplot() + 
  aes(y = hypnosis_depth, x = expectancy, color = trial_type) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x")+
  theme_bw() +
  theme(legend.position="bottom") +
  labs(
    x = "Expectancy",
    y = "Hypnosis depth",
    color = "Trial type"
  ) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
  facet_wrap(~ procedure_type) +
  scale_color_aaas(labels=trial_type_labels)

figure_3

### create high resolution image of figure
#jpeg(paste0(figure_location, "figure_3.jpg"), width = 15, height = 15, units = 'cm', res = 300)
#figure_3 # Make plot
#dev.off()


# Correlation 
corr_boot(data_long %>% 
            filter(procedure_type == "relaxation") %>% 
            select(hypnosis_depth, expectancy))
corr_boot(data_long %>% 
            filter(procedure_type == "confusion") %>% 
            select(hypnosis_depth, expectancy))
corr_boot(data_long %>% 
            filter(procedure_type == "embedded") %>% 
            select(hypnosis_depth, expectancy))
corr_boot(data_long %>% 
            filter(procedure_type == "whitenoise") %>% 
            select(hypnosis_depth, expectancy))


#######################
###   EEG analysis  ###
#######################

########################################
# In all of the analysis below, the variables are the 
# Z-transformed experience minus baseline difference values, denoted by ebds
#########################################


### Difference between conventional and placebo hypnosis 
### induction (when they are described as hypnosis) 
### in the change in the candidate EEG features

predictors = c("theta_ebds", "alpha_ebds", "gamma_ebds")
localizations = c("FZ", "OZ")
output_table_of_shamdiff = data.frame(localizations = NA, predictors = NA, mean_difference = NA, CIlb = NA, CIub = NA, n = NA) 
run_counter = 0
  
for(h in 1:length(predictors)){
  for(j in 1:length(localizations)){
    run_counter = run_counter+1
    shamdiff_mean_difference_table = data_long %>% 
      filter(description_type == "hypnosis") %>% 
      select(unique_id, trial_type, paste0(localizations[j], "_", predictors[h])) %>% 
      pivot_wider(names_from = trial_type, values_from = last_col()) %>% 
      mutate(difference = `TRUE` - sham)
    
    mean_boot_result = mean_boot(shamdiff_mean_difference_table$difference)
    mean_difference = mean_boot_result$m 
    mean_difference_CIlb = mean_boot_result$CI_lb
    mean_difference_CIub = mean_boot_result$CI_ub
    
    output_table_of_shamdiff[run_counter, "predictors"] = predictors[h]
    output_table_of_shamdiff[run_counter, "localizations"] = localizations[j]
    output_table_of_shamdiff[run_counter, "mean_difference"] = mean_difference
    output_table_of_shamdiff[run_counter, "CIlb"] = mean_difference_CIlb
    output_table_of_shamdiff[run_counter, "CIub"] = mean_difference_CIub
    output_table_of_shamdiff[run_counter, "n"] = sum(!is.na(shamdiff_mean_difference_table$difference))
  }
}

output_table_of_shamdiff





predictors_connect = c("beta1", "alpha2", "theta")
localizations_connect = c("connect_OZ_PO4_ebds", "connect_PZ_area_FZ_area_ebds", "connect_O1_PZ_ebds")
output_table_of_shamdiff_connect = data.frame(localizations = NA, predictors = NA, mean_difference = NA, CIlb = NA, CIub = NA, n = NA) 
run_counter = 0

for(h in 1:length(predictors_connect)){
    run_counter = run_counter+1
    shamdiff_mean_difference_table = data_long %>% 
      filter(description_type == "hypnosis") %>% 
      select(unique_id, trial_type, paste0(predictors_connect[h], "_", localizations_connect[h])) %>% 
      pivot_wider(names_from = trial_type, values_from = last_col()) %>% 
      mutate(difference = `TRUE` - sham)
    
    mean_boot_result = mean_boot(shamdiff_mean_difference_table$difference)
    mean_difference = mean_boot_result$m 
    mean_difference_CIlb = mean_boot_result$CI_lb
    mean_difference_CIub = mean_boot_result$CI_ub
    
    output_table_of_shamdiff_connect[run_counter, "predictors"] = predictors_connect[h]
    output_table_of_shamdiff_connect[run_counter, "localizations"] = localizations_connect[h]
    output_table_of_shamdiff_connect[run_counter, "mean_difference"] = mean_difference
    output_table_of_shamdiff_connect[run_counter, "CIlb"] = mean_difference_CIlb
    output_table_of_shamdiff_connect[run_counter, "CIub"] = mean_difference_CIub
    output_table_of_shamdiff_connect[run_counter, "n"] = sum(!is.na(shamdiff_mean_difference_table$difference))

}

output_table_of_shamdiff_connect




predictors = c("theta_ebds", "alpha_ebds", "gamma_ebds")
localizations = c("FZ", "OZ")
output_table_of_descmdiff = data.frame(localizations = NA, predictors = NA, mean_difference = NA, CIlb = NA, CIub = NA, n = NA) 
run_counter = 0

for(h in 1:length(predictors)){
  for(j in 1:length(localizations)){
    run_counter = run_counter+1
  
    descdiff_mean_difference_table_pre = data_long %>% 
      select(unique_id, trial_type, description_type, paste0(localizations[j], "_", predictors[h])) %>% 
      group_by(unique_id, description_type) %>% 
      summarise(mean = across(last_col(), mean)) %>% 
      ungroup() %>% 
      rename(mean = last_col()) %>% 
      pivot_wider(names_from = description_type, values_from = mean) %>% 
      mutate(difference = hypnosis - control)
    
    descdiff_mean_difference_table = data.frame(unique_id = descdiff_mean_difference_table_pre[,1], difference = as.matrix(descdiff_mean_difference_table_pre[,4]))

    mean_boot_result = mean_boot(descdiff_mean_difference_table$difference)
    mean_difference = mean_boot_result$m 
    mean_difference_CIlb = mean_boot_result$CI_lb
    mean_difference_CIub = mean_boot_result$CI_ub

    output_table_of_descmdiff[run_counter, "predictors"] = predictors[h]
    output_table_of_descmdiff[run_counter, "localizations"] = localizations[j]
    output_table_of_descmdiff[run_counter, "mean_difference"] = mean_difference
    output_table_of_descmdiff[run_counter, "CIlb"] = mean_difference_CIlb
    output_table_of_descmdiff[run_counter, "CIub"] = mean_difference_CIub
    output_table_of_descmdiff[run_counter, "n"] = sum(!is.na(descdiff_mean_difference_table$difference))
    
  }
}

output_table_of_descmdiff



predictors_connect = c("beta1", "alpha2", "theta")
localizations_connect = c("connect_OZ_PO4_ebds", "connect_PZ_area_FZ_area_ebds", "connect_O1_PZ_ebds")
output_table_of_descmdiff_connect = data.frame(localizations = NA, predictors = NA, mean_difference = NA, CIlb = NA, CIub = NA, n = NA) 
run_counter = 0

for(h in 1:length(predictors_connect)){
    run_counter = run_counter+1
    
    descdiff_mean_difference_table_pre = data_long %>% 
      select(unique_id, trial_type, description_type, paste0(predictors_connect[h], "_", localizations_connect[h])) %>% 
      group_by(unique_id, description_type) %>% 
      summarise(mean = across(last_col(), mean)) %>% 
      ungroup() %>% 
      rename(mean = last_col()) %>% 
      pivot_wider(names_from = description_type, values_from = mean) %>% 
      mutate(difference = hypnosis - control)
    
    descdiff_mean_difference_table = data.frame(unique_id = descdiff_mean_difference_table_pre[,1], difference = as.matrix(descdiff_mean_difference_table_pre[,4]))
    
    mean_boot_result = mean_boot(descdiff_mean_difference_table$difference)
    mean_difference = mean_boot_result$m 
    mean_difference_CIlb = mean_boot_result$CI_lb
    mean_difference_CIub = mean_boot_result$CI_ub
    
    output_table_of_descmdiff_connect[run_counter, "predictors"] = predictors_connect[h]
    output_table_of_descmdiff_connect[run_counter, "localizations"] = localizations_connect[h]
    output_table_of_descmdiff_connect[run_counter, "mean_difference"] = mean_difference
    output_table_of_descmdiff_connect[run_counter, "CIlb"] = mean_difference_CIlb
    output_table_of_descmdiff_connect[run_counter, "CIub"] = mean_difference_CIub
    output_table_of_descmdiff_connect[run_counter, "n"] = sum(!is.na(descdiff_mean_difference_table$difference))
    

}

output_table_of_descmdiff_connect

output_table_of_sham_and_desc_diff_pre = data.frame(manipulation = rep(c("trial_type", "description_type"), each = (nrow(output_table_of_shamdiff)+nrow(output_table_of_shamdiff_connect))))
  
output_table_of_sham_and_desc_diff = cbind(output_table_of_sham_and_desc_diff_pre, rbind(output_table_of_shamdiff, output_table_of_shamdiff_connect, output_table_of_descmdiff, output_table_of_descmdiff_connect))

output_table_of_sham_and_desc_diff = output_table_of_sham_and_desc_diff %>% 
  mutate(EEG_features = factor(paste0(predictors, "_", localizations), levels = c("theta_ebds_FZ", "alpha_ebds_FZ", "gamma_ebds_FZ", "theta_ebds_OZ", "alpha_ebds_OZ", "gamma_ebds_OZ", "theta_connect_O1_PZ_ebds", "alpha2_connect_PZ_area_FZ_area_ebds", "beta1_connect_OZ_PO4_ebds")))



manipulation_labels <- c("Description type", "Trial type")
names(manipulation_labels) <- c("description_type", "trial_type")

outcomes_labels <- c("Hypnosis depth", "Hypnotizability")
names(outcomes_labels) <- c("hypnosis_depth", "hypnotizability_total")

figure_4a = output_table_of_sham_and_desc_diff %>% 
  ggplot() +
  aes(y = mean_difference, x = EEG_features, color = manipulation) %>% 
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=CIlb, ymax=CIub, x = EEG_features, color = manipulation), width=.3) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(manipulation~., labeller = labeller(manipulation = outcomes_labels)) +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "EEG features",
    y = "Mean difference",
    fill = "Experimental manipulation"
  )+
  coord_flip()

figure_4a


figure_4b = data_long %>%
  select(FZ_alpha_ebds, description_type) %>% 
  drop_na() %>% 
  ggplot() +
  aes(x = description_type, y = FZ_alpha_ebds, fill = description_type) +
  
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    # outliers should not show up for the boxplot
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    # ploting on left side
    side = "left",
    # adjusting position
    justification = 1.1,
    # adjust grouping (binning) of observations
    binwidth = 0.08
  ) +
  # Themes and Labels
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "Label",
    y = "Alpha: FzArea",
    fill = "Label"
  ) +
  scale_fill_aaas() +
  coord_flip()
figure_4b



### Do these candidate EEG features help us predict whether the person underwent
### conventional or placebo hypnosis?

trial_type_mod_original = glmer(trial_type ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds + (1|ID),
                              data = data_long %>% 
                                filter(description_type == "hypnosis") %>% 
                                select(trial_type, FZ_alpha_ebds, FZ_theta_ebds, FZ_gamma_ebds, OZ_alpha_ebds, OZ_theta_ebds, OZ_gamma_ebds, beta1_connect_OZ_PO4_ebds, alpha2_connect_PZ_area_FZ_area_ebds, theta_connect_O1_PZ_ebds, ID) %>% 
                                drop_na(), family = "binomial")

glm_pred = NA
glm_pred[predict(trial_type_mod_original)<0] = "sham"
glm_pred[!(predict(trial_type_mod_original)<0)] = "TRUE"

trial_type_actual = as.vector(as.matrix(data_long %>% 
                                          filter(description_type == "hypnosis") %>% 
                                          select(trial_type, FZ_alpha_ebds, FZ_theta_ebds, FZ_gamma_ebds, OZ_alpha_ebds, OZ_theta_ebds, OZ_gamma_ebds, beta1_connect_OZ_PO4_ebds, alpha2_connect_PZ_area_FZ_area_ebds, theta_connect_O1_PZ_ebds, ID) %>% 
                                          drop_na() %>% 
                                          select(trial_type)))


# base rate of trial_type == "TRUE":
sum(trial_type_actual == "TRUE")/length(trial_type_actual)
# correct prediction rate:
correct_prediciton_rate_original = sum(glm_pred == trial_type_actual)/length(trial_type_actual)
correct_prediciton_rate_original

comparison_table = as.data.frame(cbind(glm_pred, trial_type_actual))

predictions_of_hypnosis = comparison_table %>% 
  filter(trial_type_actual == "TRUE")
true_positive_rate_original = sum(predictions_of_hypnosis$glm_pred == predictions_of_hypnosis$trial_type_actual)/nrow(predictions_of_hypnosis)
true_positive_rate_original

predictions_of_sham = comparison_table %>% 
  filter(trial_type_actual == "sham")
true_negative_rate_original = sum(predictions_of_sham$glm_pred == predictions_of_sham$trial_type_actual)/nrow(predictions_of_sham)
true_negative_rate_original



true_positive_rate_perm = NA
true_negative_rate_perm = NA

for(i in 1:perm_iter){
  print(paste0(i, "/", perm_iter))
  
  data_glm_perm = data_long %>% 
    filter(description_type == "hypnosis") %>% 
    select(trial_type, FZ_alpha_ebds, FZ_theta_ebds, FZ_gamma_ebds, OZ_alpha_ebds, OZ_theta_ebds, OZ_gamma_ebds, beta1_connect_OZ_PO4_ebds, alpha2_connect_PZ_area_FZ_area_ebds, theta_connect_O1_PZ_ebds, ID) %>% 
    drop_na() %>%
    mutate(trial_type_perm = sample(trial_type))
  
  trial_type_mod_perm = glmer(trial_type_perm ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds + (1|ID),
                                  data = data_glm_perm, family = "binomial")
  
  glm_pred = NA
  glm_pred[predict(trial_type_mod_perm)<0] = "sham"
  glm_pred[!(predict(trial_type_mod_perm)<0)] = "TRUE"
  
  trial_type_actual = as.vector(as.matrix(data_glm_perm %>% 
                                            select(trial_type_perm)))

  comparison_table = as.data.frame(cbind(glm_pred, trial_type_actual))
  
  predictions_of_hypnosis = comparison_table %>% 
    filter(trial_type_actual == "TRUE")
  true_positive_rate_perm[i] = sum(predictions_of_hypnosis$glm_pred == predictions_of_hypnosis$trial_type_actual)/nrow(predictions_of_hypnosis)

  predictions_of_sham = comparison_table %>% 
    filter(trial_type_actual == "sham")
  true_negative_rate_perm[i] = sum(predictions_of_sham$glm_pred == predictions_of_sham$trial_type_actual)/nrow(predictions_of_sham)

}

# true positive rate of the original model
true_positive_rate_original
# mean true positive rate in the permutation sample
mean(true_positive_rate_perm)
# percentage of permutation samples with as high or higher true positive rate as the original
sum(true_positive_rate_perm>=true_positive_rate_original)/length(true_positive_rate_perm)

# true negative rate of the original model
true_negative_rate_original
# mean true negative rate in the permutation sample
mean(true_negative_rate_perm)
# percentage of permutation samples with as high or higher true negative rate as the original
sum(true_negative_rate_perm>=true_negative_rate_original)/length(true_negative_rate_perm)




### How theta, alpha, and gamma power relate to hypnosis depth in conditions where
### conventional hypnosis was described as real hypnosis? We explore this in the frontal and occipital central (Z) position.

### Correlations

data_for_function = data_long %>% 
  filter(description_type == "hypnosis")
outcomes = c("hypnosis_depth", "hypnotizability_total")
predictors = c("theta_ebds", "alpha_ebds", "gamma_ebds")
localizations = c("FZ", "OZ")
trial_types = c("TRUE", "sham")

output_table_of_corfunc = data.frame(outcomes = NA, predictors = NA, localizations = NA, trial_types = NA, correlation = NA, CIlb = NA, CIub = NA)

run_counter = 0

for(h in 1:length(outcomes)){
  for(i in 1:length(predictors)){
    for(j in 1:length(localizations)){
      for(k in 1:length(trial_types)){
        run_counter = run_counter+1
        print(paste0(run_counter, "/", (length(outcomes)*length(predictors)*length(localizations)*length(trial_types))))
        
        data_r = data_for_function %>% 
          filter(trial_type == trial_types[k]) %>% 
          select(outcomes[h], paste0(localizations[j], "_", predictors[i]))
        
        cor_result = corr_boot(data_r)
        
        output_table_of_corfunc[run_counter, "outcomes"] = outcomes[h]
        output_table_of_corfunc[run_counter, "predictors"] = predictors[i]
        output_table_of_corfunc[run_counter, "localizations"] = localizations[j]
        output_table_of_corfunc[run_counter, "trial_types"] = trial_types[k]
        output_table_of_corfunc[run_counter, "correlation"] = cor_result$r
        output_table_of_corfunc[run_counter, "CIlb"] = cor_result$CI_lb
        output_table_of_corfunc[run_counter, "CIub"] = cor_result$CI_ub
      }
    }
  }
}

## table of correlations for all combinations of 
## outcomes (hypnosis_depth or hypnotizability_total),
## EEG features (alpha, theta, and gamma power),
## and localizations (FZ and OZ)

output_table_of_corfunc
# write.csv(output_table_of_corfunc, "C:\\Users\\User\\Documents\\output_table_of_corfunc.csv", row.names = F)

output_table_of_corfunc = output_table_of_corfunc %>% 
  mutate(predictors = factor(predictors, levels = c("theta_ebds", "alpha_ebds", "gamma_ebds")))
outcomes_labels <- c("Hypnosis depth", "Hypnotizability")
names(outcomes_labels) <- c("hypnosis_depth", "hypnotizability_total")
localizations_labels <- c("Region: FZ", "Region: OZ")
names(localizations_labels) <- c("FZ", "OZ")
trial_type_labels <- c("unconventional", "conventional")
names(trial_type_labels) <- c("sham", "TRUE")

figure_corr = output_table_of_corfunc %>% 
  ggplot() +
  aes(y = correlation, x = predictors, color = trial_types) %>% 
  geom_point(size = 3, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin=CIlb, ymax=CIub, x = predictors, color = trial_types), width=.3, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_discrete(labels= c("theta", "alpha", "gamma")) +
  scale_color_discrete(labels=trial_type_labels)+
  facet_grid(outcomes~localizations, labeller = labeller(outcomes = outcomes_labels, localizations = localizations_labels)) +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "EEG spectral component",
    y = "Spearman correlation (rho)",
    fill = "Outcomes",
    color = "Induction type"
  )

figure_corr


plot_a = data_for_function %>% 
  select(trial_type, hypnosis_depth, FZ_gamma_ebds) %>% 
  ggplot() + 
  aes(y = hypnosis_depth, x = FZ_gamma_ebds, color = trial_type) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se=F)+
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_discrete(labels=trial_type_labels)+
  labs(
    x = "gamma: FzArea",
    y = "Hypnosis depth",
    color = "Induction type"
  ) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10), limits = c(0,10))


plot_b = data_for_function %>% 
  select(trial_type, hypnotizability_total, FZ_gamma_ebds) %>% 
  ggplot() + 
  aes(y = hypnotizability_total, x = FZ_gamma_ebds, color = trial_type) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se=F)+
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_discrete(labels=trial_type_labels)+
  labs(
    x = "gamma: FzArea",
    y = "Hypnotizability",
    color = "Induction type"
  ) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12), limits = c(0,12))


plot_c = data_for_function %>% 
  select(trial_type, hypnotizability_total, OZ_theta_ebds) %>% 
  ggplot() + 
  aes(y = hypnotizability_total, x = OZ_theta_ebds, color = trial_type) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se=F)+
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_discrete(labels=trial_type_labels)+
  labs(
    x = "theta: OzArea",
    y = "Hypnotizability",
    color = "Induction type"
  ) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12), limits = c(0,12))

plot_d = data_for_function %>% 
  select(trial_type, hypnotizability_total, alpha2_connect_PZ_area_FZ_area_ebds) %>% 
  ggplot() + 
  aes(y = hypnotizability_total, x = alpha2_connect_PZ_area_FZ_area_ebds, color = trial_type) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se=F)+
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_discrete(labels=trial_type_labels)+
  labs(
    x = "alpha2: PzArea-FzArea connectivity",
    y = "Hypnotizability",
    color = "Induction type"
  ) +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10, 12), limits = c(0,12))


figure_8 <- ggarrange(plot_a,plot_b, plot_c, plot_d,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, 
                      common.legend = TRUE, legend="bottom")
figure_8



### create high resolution image of figure
#jpeg(paste0(figure_location, "figure_6.jpg"), width = 15, height = 15, units = 'cm', res = 300)
#figure_6 # Make plot
#dev.off()



### Connectivity correlations

predictors_connect = c("beta1_connect_OZ_PO4_ebds", "alpha2_connect_PZ_area_FZ_area_ebds", "theta_connect_O1_PZ_ebds")

output_table_of_corfunc_connect = data.frame(outcomes = NA, predictors = NA, trial_types = NA, correlation = NA, CIlb = NA, CIub = NA)

run_counter = 0

for(h in 1:length(outcomes)){
  for(i in 1:length(predictors_connect)){
    for(k in 1:length(trial_types)){
      run_counter = run_counter+1
      print(paste0(run_counter, "/", (length(outcomes)*length(predictors_connect)*length(trial_types))))
      
      data_r = data_for_function %>% 
        filter(trial_type == trial_types[k]) %>% 
        select(outcomes[h], predictors_connect[i])
      
      cor_result = corr_boot(data_r)
      
      output_table_of_corfunc_connect[run_counter, "outcomes"] = outcomes[h]
      output_table_of_corfunc_connect[run_counter, "predictors"] = predictors_connect[i]
      output_table_of_corfunc_connect[run_counter, "trial_types"] = trial_types[k]
      output_table_of_corfunc_connect[run_counter, "correlation"] = cor_result$r
      output_table_of_corfunc_connect[run_counter, "CIlb"] = cor_result$CI_lb
      output_table_of_corfunc_connect[run_counter, "CIub"] = cor_result$CI_ub
    }
  }
}

output_table_of_corfunc_connect

output_table_of_corfunc_connect = output_table_of_corfunc_connect %>% 
  mutate(predictors = factor(predictors, levels = c("theta_connect_O1_PZ_ebds", "alpha2_connect_PZ_area_FZ_area_ebds", "beta1_connect_OZ_PO4_ebds")))
outcomes_labels <- c("Hypnosis depth", "Hypnotizability")
names(outcomes_labels) <- c("hypnosis_depth", "hypnotizability_total")

figure_corr_connect = output_table_of_corfunc_connect %>% 
  ggplot() +
  aes(y = correlation, x = predictors, color = trial_types) %>% 
  geom_point(size = 3, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin=CIlb, ymax=CIub, x = predictors, color = trial_types), width=.3, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_discrete(labels= c("theta: O1-PZ", "alpha2: PZ-FZ", "beta1 OZ-PO4")) +
  scale_color_discrete(labels=trial_type_labels)+
  facet_grid(outcomes~., labeller = labeller(outcomes = outcomes_labels)) +
  theme_bw() +
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(
    x = "EEG spectral component",
    y = "Pearson correlation (r)",
    fill = "Outcomes",
    color = "Induction type"
  )

figure_corr_connect

### Predictive models with EEG predictors

## Predicting hypnosis depth with EEG predictors

# true hypnosis

mod_EEG_lm_hypnosis_depth_true = lm(hypnosis_depth ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                                    data = data_long %>% 
                                      filter(description_type == "hypnosis") %>% 
                                      filter(trial_type == "TRUE"))
original_r2_EEG_lm_hypnosis_depth_true = summary(mod_EEG_lm_hypnosis_depth_true)$r.squared

r2_perm_EEG_lm_hypnosis_depth_true = NA
for(i in 1:perm_iter){
  print(paste0(i, "/", perm_iter))
  mod_perm = lm(perm_hypnosis_depth ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                data = data_long %>% 
                  filter(description_type == "hypnosis") %>% 
                  filter(trial_type == "TRUE") %>%
                  mutate(perm_hypnosis_depth = sample(hypnosis_depth)))
  
  r2_perm_EEG_lm_hypnosis_depth_true[i] = summary(mod_perm)$r.squared
}

# R2 of the model
original_r2_EEG_lm_hypnosis_depth_true

# mean R2 in the permutation sample
mean(r2_perm_EEG_lm_hypnosis_depth_true)

# percentage of permutation samples with as high or higher R2 than the original R2
sum(r2_perm_EEG_lm_hypnosis_depth_true>=original_r2_EEG_lm_hypnosis_depth_true)/length(r2_perm_EEG_lm_hypnosis_depth_true)


# sham hypnosis

mod_EEG_lm_hypnosis_depth_sham = lm(hypnosis_depth ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                                    data = data_long %>% 
                                      filter(description_type == "hypnosis") %>% 
                                      filter(trial_type == "sham"))
original_r2_EEG_lm_hypnosis_depth_sham = summary(mod_EEG_lm_hypnosis_depth_sham)$r.squared

r2_perm_EEG_lm_hypnosis_depth_sham = NA
for(i in 1:perm_iter){
  print(paste0(i, "/", perm_iter))
  mod_perm = lm(perm_hypnosis_depth ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                data = data_long %>% 
                  filter(description_type == "hypnosis") %>% 
                  filter(trial_type == "TRUE") %>%
                  mutate(perm_hypnosis_depth = sample(hypnosis_depth)))
  
  r2_perm_EEG_lm_hypnosis_depth_sham[i] = summary(mod_perm)$r.squared
}


# R2 of the model
original_r2_EEG_lm_hypnosis_depth_sham

# mean R2 in the permutation sample
mean(r2_perm_EEG_lm_hypnosis_depth_sham)

# percentage of permutation samples with as high or higher R2 than the original R2
sum(r2_perm_EEG_lm_hypnosis_depth_sham>=original_r2_EEG_lm_hypnosis_depth_sham)/length(r2_perm_EEG_lm_hypnosis_depth_sham)




## Predicting hypnotizability with EEG predictors

# true hypnosis

mod_EEG_lm_hypnotizability_total_true = lm(hypnotizability_total ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                                           data = data_long %>% 
                                             filter(description_type == "hypnosis") %>% 
                                             filter(trial_type == "TRUE"))
original_r2_EEG_lm_hypnotizability_total_true = summary(mod_EEG_lm_hypnotizability_total_true)$r.squared

r2_perm_EEG_lm_hypnotizability_total_true = NA
for(i in 1:perm_iter){
  print(paste0(i, "/", perm_iter))
  mod_perm = lm(perm_hypnotizability_total ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                data = data_long %>% 
                  filter(description_type == "hypnosis") %>% 
                  filter(trial_type == "TRUE") %>%
                  mutate(perm_hypnotizability_total = sample(hypnotizability_total)))
  
  r2_perm_EEG_lm_hypnotizability_total_true[i] = summary(mod_perm)$r.squared
}


# R2 of the model
original_r2_EEG_lm_hypnotizability_total_true

# mean R2 in the permutation sample
mean(r2_perm_EEG_lm_hypnotizability_total_true)

# percentage of permutation samples with as high or higher R2 than the original R2
sum(r2_perm_EEG_lm_hypnotizability_total_true>=original_r2_EEG_lm_hypnotizability_total_true)/length(r2_perm_EEG_lm_hypnotizability_total_true)


# sham hypnosis

mod_EEG_lm_hypnotizability_total_sham = lm(hypnotizability_total ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                                           data = data_long %>% 
                                             filter(description_type == "hypnosis") %>% 
                                             filter(trial_type == "sham"))
original_r2_EEG_lm_hypnotizability_total_sham = summary(mod_EEG_lm_hypnotizability_total_sham)$r.squared

r2_perm_EEG_lm_hypnotizability_total_sham = NA
for(i in 1:perm_iter){
  print(paste0(i, "/", perm_iter))
  mod_perm = lm(perm_hypnotizability_total ~ FZ_alpha_ebds + FZ_theta_ebds + FZ_gamma_ebds + OZ_alpha_ebds + OZ_theta_ebds + OZ_gamma_ebds + beta1_connect_OZ_PO4_ebds + alpha2_connect_PZ_area_FZ_area_ebds + theta_connect_O1_PZ_ebds,
                data = data_long %>% 
                  filter(description_type == "hypnosis") %>% 
                  filter(trial_type == "TRUE") %>%
                  mutate(perm_hypnotizability_total = sample(hypnotizability_total)))
  
  r2_perm_EEG_lm_hypnotizability_total_sham[i] = summary(mod_perm)$r.squared
}

# R2 of the model
original_r2_EEG_lm_hypnotizability_total_sham

# mean R2 in the permutation sample
mean(r2_perm_EEG_lm_hypnotizability_total_sham)

# percentage of permutation samples with as high or higher R2 than the original R2
sum(r2_perm_EEG_lm_hypnotizability_total_sham>=original_r2_EEG_lm_hypnotizability_total_sham)/length(r2_perm_EEG_lm_hypnotizability_total_sham)










