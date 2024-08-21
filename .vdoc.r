#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats, MASS,
fitdistrplus, DHARMa)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-Methods
#| fig.cap: "Experimental design of early environment manipulation (**A**), numerical discrimination arena an measures(**B**), measures of the platform used for the experiments (**C**), and types of numerical trial (**D**), where each line indicates the position and orientation of the food items in each test."

knitr::include_graphics("./Others/NUMDIS_FIG_1.png")

#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: data_process
# Load the data and process it to get the main dfs
source(here("R", "data_process.R"))
#
#
#
#| label: sampleSize
# List with the sample sizes from the database (here("output/databases_clean/data_asso.csv") as the sample size per species and group is the same on each task. We used function sample (see func.R) to estimate the sample size per treatment and species.
source(here("R", "func.R"))
#
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
side <- c("Right", "Left")
#
n_list <- list()
#
for(k in 1:length(hormone)){
  for(l in 1:length(temperature)){
    for(i in 1:length(side)){
      n_list[[paste(hormone[k], temperature[l], side[i], sep = "_")]] <- sample(hormone[k], temperature[l], side [i])
    }
  }
}
#
#
#
#| label: models
# Fitting the model and extraction of posteriors 
refit = FALSE
if(refit){
  # Get the formula and family of variables
  lat <- bf(log(latency) ~ test_type*temp*cort + sex + age + (1|lizard_id) + (1|clutch), family = gaussian())
  choice <- bf(choice ~ test_type*temp*cort + sex + age + (1|lizard_id) + (1|clutch), family = bernoulli(link = "logit"))
  int <- bf(compared_interest ~ test_type*temp*cort + sex + age + (1|lizard_id) + (1|clutch),
    family = gaussian())
  # Fit the model
  model <- brm(lat + choice + int,
                chains = 4, cores = 4, iter = 3000, warmup = 1000, control = list(adapt_delta = 0.99, max_treedepth =12),
                data = data_num)
    # Write the model to a file
    saveRDS(model, file = here("output/models/model_num.rds"))
  } else {
      # Read the model from a file
      model <- readRDS(file = here("output/models/model_num.rds"))
  } 
# 
#
#
#
#| label: posteriors
# Code to extract posteriors from the model and get predicted values
source(here("R", "func.R"))
#
## Extract posteriors
posteriors <- as_draws_df(model)
# 
## Get predictions for every test and treatment based on the posteriors of the model. We used the function org_post (see func.R) to get the predictions.
#
# Latency
post_lat <- org_post("lat")
# Choice
post_choice <- org_post("choice")
prob_choice <- exp(post_choice)/(1 + exp(post_choice))
# Estimated interest for the higher amount of food
post_int <- org_post("int")
#
#
#
#
#
#
#
#
#
#| label: tbl-data
#| tbl-cap: "Estimates of the effects of each treatment and their interaction on the latency, choice, and interest shown for the higher amount of food in each of the numerical discrimination tests. The table shows the comparissons for each predictor and the p~mcmcm~ value of the comparissons. In bold, the estimates with p~mcmc~ < 0.05."
source(here("R", "func.R"))
#
# Making the df for each variable
#
## Initialize df_org as a list to store data frames for each variable
df_org <- list()
## Define variables and tests
variable <- c("Latency", "Choice", "Interest")
test <- c("1VS4", "1VS3", "2VS4", "2VS3", "3VS4")
## Loop over each variable
for(a in variable){
  # Select the original df to use depending on the variable
  if(a == "Latency"){
    post <- post_lat
  } else if(a == "Choice"){
    post <- prob_choice
  } else if(a == "Interest"){
    post <- post_int
  }
  df_tables <- data.frame()
  for(b in test){
    # Modify the df to get the values per each test
    post_mod <- post %>%
      dplyr::select(dplyr::matches(b))
    # Extract the required columns
    Cold_CORT <- post_mod[[paste0("X", b, "_Cold_CORT")]]
    Cold_Control <- post_mod[[paste0("X", b, "_Cold_Control")]]
    Hot_CORT <- post_mod[[paste0("X", b, "_Hot_CORT")]]
    Hot_Control <- post_mod[[paste0("X", b, "_Hot_Control")]]
    # Calculate effects and pmcmc values
    effect_temp <- format_dec(mean(c(Hot_CORT, Hot_Control)) - mean(c(Cold_CORT, Cold_Control)), 3)
    effect_hormone <- format_dec(mean(c(Hot_Control, Cold_Control)) - mean(c(Hot_CORT, Cold_CORT)), 3)
    interaction <- format_dec((mean(Hot_Control) - mean(Cold_Control)) - (mean(Hot_CORT) - mean(Cold_CORT)), 3)
    #
    pmcmc_temp <- format_p(pmcmc(c(Hot_CORT, Hot_Control) - c(Cold_CORT, Cold_Control)), 3)
    pmcmc_hormone <- format_p(pmcmc(c(Hot_Control, Cold_Control) - c(Hot_CORT, Cold_CORT)), 3)
    pmcmc_interaction <- format_p(pmcmc((Hot_Control - Cold_Control) - (Hot_CORT - Cold_CORT)), 3)
    # Define labels and values
    effects_predictors_names <- c("Temperature", 
                          "Hormone", 
                          "Interaction")
    values_predictors <- c(effect_temp, effect_hormone, interaction)
    values_pmcmc <- c(pmcmc_temp, pmcmc_hormone, pmcmc_interaction)
    # Create a temporary data frame for this loop iteration
    temporal_df <- data.frame(Predictor = effects_predictors_names,
                        Effect = values_predictors,
                        pmcmc = values_pmcmc,
                        Test = b,
                        stringsAsFactors = FALSE)
    # Merge the temporary data frame to the main df_tables
    df_tables <- rbind(df_tables, temporal_df)
  }
  df_tables$Variable <- a
  df_org[[a]] <- df_tables
}
# Combine all data frames in df_org into one data frame
table_df <- do.call(rbind, df_org)
# Convert the combined list into a data frame and modify for the table
table_df <- as.data.frame(table_df, stringsAsFactors = FALSE) %>%
  mutate(results = paste(Effect, "(p", pmcmc, ")")) %>%
  dplyr::select(Variable, Test, Predictor, results) %>%
  pivot_wider(names_from = Test, values_from = c(results), names_prefix = "")
#
#
# Making the table
#
## Table format
set_flextable_defaults(
 font.family = "Arial",
 fint.size = 10)
# Split the table_data df by task
real_table <- flextable(table_df) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    add_header_row(values = c("", "Tests"), colwidths = c(2, 5)) %>%
    flextable::compose(i = c(2,3,5,6,8,9), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    hline(i = c(3,6), j = c(1:7), part = "body") %>% # To make some horizontal lines
    vline(i = (1:9), j = c(2:6), part = "body") %>% # To make some vertical lines on body
    vline(i = 2, j=c(2:6), part = "header") %>% # To make some vertical lines on header
    autofit() 
real_table
#
#
#
#
#
#| label: fig-results
#| fig.cap: "" 
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#| label: tbl-summary
#| tbl-cap: "Summary of the model fitted to the data. The table shows the estimates of the effects of each treatment and their interaction on the latency, choice, and interest shown for the higher amount of food in each of the numerical discrimination tests. The table shows the comparissons for each predictor and the p~mcmcm~ value of the comparissons. In bold, the estimates with p~mcmc~ < 0.05."
```{r, results='asis', echo=FALSE}
cat("\\newpage")
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#| label: tbl-sides
source(here("R", "func.R"))
#
# 1) MAKING THE DF FOR THE TABLE
data_sides <- data_num %>%
  group_by(trt, test_type) %>% # Group by treatment and test
  summarise(
    right_choices = sum(side_choice == "Right", na.rm = TRUE),
    left_choices = sum(side_choice == "Left", na.rm = TRUE),
  ) %>% # Count the number of correct choices
  rowwise() %>%  # Ensure that mutate works row by row
  mutate(
    total_choices = right_choices + left_choices,  # Calculate total choices
    p_value = format_dec(binom.test(right_choices, total_choices)$p.value, 2)  # Perform binom test
  ) %>% # Test for differences between sides
  ungroup() %>%  # Remove rowwise grouping
  mutate(choices_Right_Left_p = paste0("R = ", right_choices, " | L = ", left_choices, " , p = ", p_value), # Create a column with the choices
        trt = recode(trt, 
                 "B_23" = "CORT-Cold (n = 20)", 
                 "A_23" = "Control-Cold (n = 20)", 
                 "B_28" = "CORT-Hot (n = 20)", 
                 "A_28" = "Control-Hot (n = 20)")) %>% # To order the treatments
  dplyr::select(trt, test_type, choices_Right_Left_p) %>% # Remove unused columns
  pivot_wider(names_from = test_type, values_from = c(choices_Right_Left_p)) %>% # To split between tests
data.frame()
#
# 2) MAKING THE TABLE
#
# Split the table_data df by task
sides_table <- flextable(data_sides) %>%
  compose(j = 1, value = as_paragraph(""), part = "header") %>%
  compose(j = 2, value = as_paragraph("1 VS 4"), part = "header") %>%
  compose(j = 3, value = as_paragraph("1 VS 3"), part = "header") %>%
  compose(j = 4, value = as_paragraph("2 VS 4"), part = "header") %>%
  compose(j = 5, value = as_paragraph("2 VS 3"), part = "header") %>%
  compose(j = 6, value = as_paragraph("3 VS 4"), part = "header") %>%
  set_table_properties( width = 0.9, layout = "autofit") %>%
  align(align = "center", part = "all")
sides_table
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#| label: tbl-sides
source(here("R", "func.R"))
#
# 1) MAKING THE DF FOR THE TABLE
data_or_res <- data_orient %>%
  mutate(choice_or = ifelse(choice_horizontal == 0, "Vertical", "Horizontal")) %>%
  dplyr::select(trt, choice_horizontal, choice_or) %>%
  group_by(trt) %>% # Group by treatment and test
  summarize(
    Horizontal_Count = sum(choice_or == "Horizontal"),
    Vertical_Count = sum(choice_or == "Vertical")
  ) %>% # Count the number of horizontal and vertical choices
  rowwise() %>%  # Ensure that mutate works row by row
  mutate(
    total_choices = Horizontal_Count + Vertical_Count,  # Calculate total choices
    p_value = format_dec(binom.test(Horizontal_Count, total_choices)$p.value, 3)  # Perform binom test
  ) %>% # Test for differences between sides
  ungroup() %>%  # Remove rowwise grouping
  mutate(trt = recode(trt, 
            "B_23" = "CORT-Cold (n = 20)", 
            "A_23" = "Control-Cold (n = 20)", 
            "B_28" = "CORT-Hot (n = 20)", 
            "A_28" = "Control-Hot (n = 20)")) %>% # To order the treatments
  dplyr::select(trt, Horizontal_Count, Vertical_Count, p_value) %>% # Remove unused columns
data.frame()
#
# 2) MAKING THE TABLE
#
# Split the table_data df by task
orient_table <- flextable(data_or_res) %>%
  compose(j = 1, value = as_paragraph(""), part = "header") %>%
  compose(j = 2, value = as_paragraph("Horizontal"), part = "header") %>%
  compose(j = 3, value = as_paragraph("Vertical"), part = "header") %>%
  compose(j = 4, value = as_paragraph("p-value"), part = "header") %>%
  set_table_properties( width = 0.9, layout = "autofit") %>%
  align(align = "center", part = "all")
orient_table
#
#
#
cat("\\newpage")
#
#
#
#
# Chunk for plotting the model
## Initial plot
plot(model)
#
## pp_checks
pp_check(model, resp = "loglatency")
pp_check(model, resp = "choice")
pp_check(model, resp = "comparedinterest")
#
#
#
#
