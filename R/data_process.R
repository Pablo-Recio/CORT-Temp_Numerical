####################################
# 1_data_process
####################################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats)
# A) We extract the lizards we are using from the main data_base
meta_data <- read.csv("./data/lizard_database.csv") %>%
  filter(lizard_id %in% c("LG2859_23", "LG2893_23", "LG2882_23", "LG2915_23", "LG2946_23", "LG2917_23", "LG2811_23",
  "LG3182_23", "LG2975_23", "LG2869_23", "LG2898_23", "LG2897_23", "LG3205_23", "LG2812_23", "LG2934_23", "LG3208_23",
  "LG3035_23", "LG3151_23", "LG3076_23", "LG2866_23", "LG2831_23", "LG3204_23", "LG2976_23", "LG2894_23", "LG3186_23",
  "LG2810_23", "LG2947_23", "LG3022_23", "LG2941_23", "LG3150_23", "LG3161_23", "LG2991_23", "LG2880_23", "LG2914_23",
  "LG2977_23", "LG2923_23", "LG2943_23", "LG2935_23", "LG2844_23", "LG2932_23", "LG3126_23", "LG2846_23", "LG3121_23",
  "LG3119_23", "LG3168_23", "LG2813_23", "LG3075_23", "LG3175_23", "LG3209_23", "LG3117_23", "LG3134_23", "LG3146_23",
  "LG3167_23", "LG3077_23", "LG3172_23", "LG3128_23", "LG3002_23", "LG3130_23", "LG2998_23", "LG3059_23", "LG3073_23",
  "LG3133_23", "LG2817_23", "LG2988_23", "LG3000_23", "LG3019_23", "LG3061_23", "LG3004_23", "LG3058_23", "LG2972_23",
  "LG2978_23", "LG3020_23", "LG3037_23", "LG2933_23")) %>%
data.frame()
write.csv(meta_data, "./data/numerical.csv") 
#
# B) Clean the meta_data
# Remove individuals who did not participate (more than 15 NAs), remove trials [36-40] (Only in associative) and split treatment into Temp and Cort
data_asso <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_associative)) <= 15) %>%
    filter(trial_associative <= 35) %>%
  ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
    mutate( # Standarize data by trial (i.e. make the first trial where each individual participated their trial 1)
      first_non_na = min(which(!is.na(FC_associative))),
      Associative_Trial = ifelse(!is.na(first_non_na),trial_associative - first_non_na + 1, trial_associative))%>% 
      filter(Associative_Trial >= 1) %>%
  ungroup() %>% 
data.frame()


data_rev <- data %>%
  group_by(lizard_id) %>%
    filter(sum(is.na(FC_associative)) <= 15) %>%
    filter(sum(is.na(FC_reversal)) <= 15) %>%
   ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
    mutate(trial_reversal=as.numeric(trial_reversal)) %>%
  data.frame()


