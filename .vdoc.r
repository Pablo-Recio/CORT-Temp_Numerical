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
#
## Extract posteriors
posteriors <- as_draws_df(model)
# 
## Get predictions for every test and treatment based on the posteriors of the model. We used the function org_post (see func.R) to get the predictions.

## A) Get predictions
post_lat <- org_post("lat")
post_choice <- org_post("choice")
post_int <- org_post("int")
#
#
#
#
#| label: tbl-data
#| tbl-cap: "Estimates of Associative learning slope for all the different treatments per each task, species and group. Mean shows the aritmetic mean of the estimates obtained from the posteriors of the model, and 95% CI indicates the 95% confidence interval of the mean. All p-values were obtained using pmcmc and test the hypothesis that the mean is equal to zero. In bold, those values that are significant (p-value <0.05)"
source(here("R", "func.R"))
#
############################## MAKING THE TABLE ##############################
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 fint.size = 10)
# Split the table_data df by task
real_table <- flextable(new_table_data) %>%
    bold(~ `p-value_Associative` < 0.05, ~ `p-value_Associative` + Mean_Associative + `95% CI_Associative`) %>%
    bold(~ `p-value_Reversal` < 0.05, ~ `p-value_Reversal` + Mean_Reversal + `95% CI_Reversal`) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    add_header_row(values = c("", "Associative task", "Reversal task"), colwidths = c(3, 3, 3)) %>%
    set_header_labels(Mean_Associative = "Mean",
                      `95% CI_Associative` = "95% CI",
                      `p-value_Associative` = "p-value",
                      Mean_Reversal = "Mean",
                      `95% CI_Reversal` = "95% CI",
                      `p-value_Reversal` = "p-value") %>%
    italic(j = 1, italic = TRUE, part = "body") %>% # To have names od species in italics
    flextable::compose(i = c(2:8,10:16), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::compose(i = c(2:4,6:8,10:12,14:16), j = 2, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the second column
    hline(i = c(4,12), j = c(2:9), part = "body") %>% # To make some horizontal lines
    hline(i = c(8), j = c(1:9), part = "body") %>% # To make some horizontal lines
    vline(i = (1:16), j = c(3,6), part = "body") %>% # To make some vertical lines on body
    vline(j=c(3,6), part = "header") %>% # To make some vertical lines on header
    autofit() 
real_table
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
cat("\\newpage")
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
    p_value = format_dec(binom.test(Horizontal_Count, total_choices)$p.value, 2)  # Perform binom test
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
