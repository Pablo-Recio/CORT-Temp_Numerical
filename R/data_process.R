####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
# Clean the data
# Remove individuals who did not participate (more than 2 NAs) and split treatment into Temp and Cort
data_num <- read.csv(here("data/lizard_database_number.csv")) %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(choice)) <= 3) %>%
  ungroup()  %>%
  mutate(temp = gsub("[AB]_", "", trt),
         cort = gsub("_[2][38]", "", trt))  %>%
  mutate(temp = factor(temp,
    levels = c("23", "28"),
    labels = c("23" = "Cold", "28" = "Hot"))) %>%
  mutate(cort = factor(cort,
    levels = c("B", "A"),
    labels = c("B" = "CORT", "A" = "Control"))) %>%
  data.frame()


