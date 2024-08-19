####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats, MASS,
fitdistrplus)
# 
#### CLEAN THE MAIN DATABASE
# Remove individuals who did not participate (more than 2 NAs), split treatment into Temp and Cort, and relevel test_type
data_num <- read.csv(here("data/lizard_database_number.csv")) %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(choice)) <= 2) %>%
  ungroup()  %>%
  mutate(temp = gsub("[AB]_", "", trt),
         cort = gsub("_[2][38]", "", trt))  %>%
  mutate(temp = factor(temp,
    levels = c("23", "28"),
    labels = c("23" = "Cold", "28" = "Hot"))) %>%
  mutate(cort = factor(cort,
    levels = c("B", "A"),
    labels = c("B" = "CORT", "A" = "Control"))) %>%
  mutate(test_type = factor(test_type,
    levels = c("1 VS 4", "1 VS 3", "2 VS 4", "2 VS 3", "3 VS 4"))) %>%
  mutate(age = scale(age, center = TRUE, scale =FALSE)) %>%
data.frame()
#
#
#
#### CLEAN THE DATABASE FOR THE ORIENTATION TEST
data_orient <- read.csv(here("data/orientation_test.csv")) %>%
  mutate(temp = gsub("[AB]_", "", trt),
         cort = gsub("_[2][38]", "", trt))  %>%
  mutate(temp = factor(temp,
    levels = c("23", "28"),
    labels = c("23" = "Cold", "28" = "Hot"))) %>%
  mutate(cort = factor(cort,
    levels = c("B", "A"),
    labels = c("B" = "CORT", "A" = "Control"))) %>%
data.frame()

