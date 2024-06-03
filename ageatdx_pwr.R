library(tidyverse)
library(faux)
library(tidymodels)

# Source the supporting functions

source("./ageatdx_pwr_funs.R")

# Suppress dplyr summarise info

options(dplyr.summarise.inform = FALSE)

# Set seed for reproducibility

set.seed(134906)

# Run a power simulation for the analysis

# These are the options:

# n=1400 #Sample size
# edu=-0.2 #This is the standardised effect of confounder education on age at dx
# fun=0.4 #This is the standardised effect of confounder functioning on age at dx
# acc=-0.1 #This is the standardised effect of confounder access to services on age at dx
# paternal_eff = -0.1 #This is the standardised effect of age_at_dx on paternal wellbeing
# maternal_eff = -0.2 #This is the standardised effect of age_at_dx on maternal wellbeing
# sib_eff = -0.05 #This is the standardised effect of age_at_dx on sibling wellbeing
# ind_eff = -0.2 #This is the standardised effect of age_at_dx on individual wellbeing
# htest_level="family" #How is H3 specified? Alternatives are "person" and "outcome"
# n_replicates=100 #How many times should the scenario be simulated

# all apart from n_replicates can be specified as a vector of options to loop over

# Example 

sim1 <- pwr_sim(ind_eff=seq(0.01,0.1,0.01))

# Overall power

sim1 %>% select(ind_eff, overall_power)

# Power to detect effects in the individual

sim1$person_power %>%  
  purrr::reduce(bind_rows) %>% 
  filter(person == "self") %>% 
  bind_cols(sim1 %>% select(ind_eff)) %>% 
  select(ind_eff,power)




