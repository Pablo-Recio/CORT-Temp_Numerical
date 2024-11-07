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
#| label: packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, forcats, MASS, nortest, fitdistrplus, ggh4x, PupillometryR, cowplot, png, grid, remotes, ggthemes, bayestestR, HDInterval)
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
#| fig.cap: "Experimental design. Panel (**A**) shows the early environment manipulation procedures. (**B**) illustrates the arena where the tests were performed, while panel (**C**) indicates the measurements of the platform used for the experiments. Finally, panel (**D**) displays the types of numerical tests used and the orientation of the crickets in each test."

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
                chains = 4, cores = 4, iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99, max_treedepth =12),
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
normality_results_lat <- check_normality(post_lat)
# Choice
post_choice <- org_post("choice")
prob_choice <- exp(post_choice)/(1 + exp(post_choice))
normality_results_choice <- check_normality(prob_choice)
write.csv(normality_results_choice, file = here("output/Checking/normality_results.csv"))
# Estimated interest for the higher amount of food
post_int <- org_post("int")
normality_results_int <- check_normality(post_int)
#
#
#
#
#
#
#| label: fig-results
#| fig.cap: "Estimates of log-latency (Latency), the probability of choosing the larger number of crickets first (Choice), and the estimated interest in the larger number of crickets (Interest) for each of the numerical tests performed across different developmental treatments. The x-axis represents the estimate, and the y-axis is the posterior density of the estimates. The different colours indicate the different treatments. Points and bars represent the median and 95% Higest Density Intervals (95% HPDI) of the estimates, respectively. Vertical dashed lines in Choice and Interest graphs values 0.5 and 0, respectively. Asterisks indicate values significantly different from 0."
#
source(here("R", "func.R"))
#
#### A) Organise the df used for both parts of the figure for each variable
#
# First the df for violin plots for each variable using df_violin_plot function from func.r
df_viol_lat <- df_violin_plot("lat")
df_viol_choice <- df_violin_plot("choice")
df_viol_int <- df_violin_plot("int")
#
# Second the df for the points and bars plot for each variable using df_bars_plot function from func.r
df_bars_lat <- df_bars_plot("lat")
df_bars_choice <- df_bars_plot("choice")
df_bars_int <- df_bars_plot("int")
#
#### B) Make the plots for each variable using plotting function from func.r
#
fig_lat <- plotting("lat")
fig_choice <- plotting("choice") 
fig_int <- plotting("int")
#
#### C) Combine the plots and make final figure
fig_results <- plot_grid(fig_lat, fig_choice, fig_int, ncol = 1) +
  annotate("text", x = 0.305, y = 0.56, label = "*", size = 7, color = "black") + 
  annotate("text", x = 0.61, y = 0.505, label = "*", size = 7, color = "black")
ggsave("./output/figures/fig_results.png", plot=fig_results, width = 25, height = 25, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_results.png")
#
#
#
#
#
#
#| label: tbl-data
#| tbl-cap: "Effects of temperature, CORT, and their interaction on Latency, Choice, and Interest in each of the numerical discrimination tests. The table shows the contrasts for each predictor (Temperature = [medianHot - medianCold]; Hormone = [medianControl - medianCORT]; and their Interaction = [(medianHot-Control - medianHot-CORT) - (medianCold-Control - medianCold-CORT)]). 95% Higest Density Intervals (95% HPDI) test the hypothesis that contrasts are different from zero."
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
    # Calculate effects and HDI values
    effect_hormone <- format_dec(median(c(Hot_Control, Cold_Control)) - median(c(Hot_CORT, Cold_CORT)), 1)
    effect_temp <- format_dec(median(c(Hot_CORT, Hot_Control)) - median(c(Cold_CORT, Cold_Control)), 1)
    interaction <- format_dec((median(Hot_Control) - median(Cold_Control)) - (median(Hot_CORT) - median(Cold_CORT)), 1)
    #
    hdi_hormone <- paste0("[", format_dec(hdi(c(Hot_Control, Cold_Control) - c(Hot_CORT, Cold_CORT))[1], 2), " , ",
                          format_dec(hdi(c(Hot_Control, Cold_Control) - c(Hot_CORT, Cold_CORT))[2], 2), "]")
    hdi_temp <- paste0("[", format_dec(hdi(c(Hot_CORT, Hot_Control) - c(Cold_CORT, Cold_Control))[1], 2), " , ",
                        format_dec(hdi(c(Hot_CORT, Hot_Control) - c(Cold_CORT, Cold_Control))[2], 2), "]")
    hdi_interaction <- paste0("[", format_dec(hdi((Hot_Control - Cold_Control) - (Hot_CORT - Cold_CORT))[1], 2), " , ",
                              format_dec(hdi((Hot_Control - Cold_Control) - (Hot_CORT - Cold_CORT))[2], 2), "]")
    # Define labels and values
    effects_predictors_names <- c("Hormone",
                          "Temperature",
                          "Interaction")
    values_predictors <- c(effect_hormone, effect_temp, interaction)
    values_hdi <- c(hdi_hormone, hdi_temp, hdi_interaction)
    # Create a temporary data frame for this loop iteration
    temporal_df <- data.frame(Predictor = effects_predictors_names,
                        Effect = values_predictors,
                        hdi = values_hdi,
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
  mutate(results = trimws(paste0(Effect, "\n", hdi))) %>%
  dplyr::select(Variable, Test, Predictor, results) %>%
  pivot_wider(names_from = Test, values_from = c(results), names_prefix = "")
#
#
# Making the table
#
## Table format
set_flextable_defaults(
 font.family = "Arial",
 font.size = 9)
# Create the table
real_table <- flextable(table_df) %>%
    align(align="center", part="all") %>% 
    add_header_row(values = c("", "Tests"), colwidths = c(2, 5)) %>%
    flextable::compose(i = c(3,6,8,9), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::compose(i = 2, j = 1, value = as_paragraph("log(latency)"), part = "body") %>%
    flextable::compose(i = 5, j = 1, value = as_paragraph("log(odds)"), part = "body") %>%
    hline(i = c(3,6), j = c(1:7), part = "body") %>% # To make some horizontal lines 
    set_table_properties(layout = "autofit", width = 1)
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
#| label: tbl-disc
#| tbl-cap: "Median and 95% HPDI of the estimated probability of choosing first the higher amount (Choice) and the estimated interest for the higher amount of food (Interest) per treatment group for each of the numerical tests performed. 95% HPDI test the hypothesis that Choice = 0.5, and Interest = 0, which would indicate a preference towards one of choices."
#
source(here("R", "func.R"))
#
# Making the df for each variable
#
## Initialize df_org as a list to store data frames for each variable
disc_list <- list()
## Define variables and tests
variable <- c("Choice", "Interest")
test <- c("1VS4", "1VS3", "2VS4", "2VS3", "3VS4")
## Loop over each variable
for(c in variable){
  # Select the original df to use depending on the variable
  if(c == "Choice"){
    post <- prob_choice
    val <- 0.5
  } else if(c == "Interest"){
    post <- post_int
    val <- 0
  }
  df_disc <- data.frame()
  for(d in test){
    # Modify the df to get the values per each test
    post_mod <- post %>%
      dplyr::select(dplyr::matches(d))
    # Extract the required columns
    Cold_CORT <- post_mod[[paste0("X", d, "_Cold_CORT")]]
    Cold_Control <- post_mod[[paste0("X", d, "_Cold_Control")]]
    Hot_CORT <- post_mod[[paste0("X", d, "_Hot_CORT")]]
    Hot_Control <- post_mod[[paste0("X", d, "_Hot_Control")]]
    # Calculate medians and 95% HPDI
    Cold_CORT_median <- format_dec(median(Cold_CORT), 2)
    Cold_Control_median <- format_dec(median(Cold_Control), 2)
    Hot_CORT_median <- format_dec(median(Hot_CORT), 2)
    Hot_Control_median <- format_dec(median(Hot_Control), 2)
    #
    int_Cold_CORT <- paste0("[", format_dec(hdi(Cold_CORT, credMass = 0.95)[1], 2), " , ",
                            format_dec(hdi(Cold_CORT, credMass = 0.95)[2], 2), "]")
    int_Cold_Control <- paste0("[", format_dec(hdi(Cold_Control, credMass = 0.95)[1], 2), " , ",
                            format_dec(hdi(Cold_Control, credMass = 0.95)[2], 2), "]")
    int_Hot_CORT <- paste0("[", format_dec(hdi(Hot_CORT, credMass = 0.95)[1], 2), " , ",
                            format_dec(hdi(Hot_CORT, credMass = 0.95)[2], 2), "]")
    int_Hot_Control <- paste0("[", format_dec(hdi(Hot_Control, credMass = 0.95)[1], 2), " , ",
                            format_dec(hdi(Hot_Control, credMass = 0.95)[2], 2), "]")
    # Define labels and values
    treatment_names <- c("Cold-CORT (n = 20)", "Cold-Control", "Hot-CORT", "Hot-Control")
    values_medians <- c(Cold_CORT_median, Cold_Control_median, Hot_CORT_median, Hot_Control_median)
    values_int <- c(int_Cold_CORT, int_Cold_Control, int_Hot_CORT, int_Hot_Control)
    # Create a temporary data frame for this loop iteration
    temporal_df <- data.frame(Treatment = treatment_names,
                        Effect = values_medians,
                        int = values_int,
                        Test = d,
                        stringsAsFactors = FALSE)
    # Merge the temporary data frame to the main df_tables
    df_disc <- rbind(df_disc, temporal_df)
  }
  df_disc$Variable <- c
  disc_list[[c]] <- df_disc
}
# Combine all data frames in df_org into one data frame
disc_df <- do.call(rbind, disc_list)
# Convert the combined list into a data frame and modify for the table
disc_df_final <- as.data.frame(disc_df, stringsAsFactors = FALSE) %>%
  mutate(int = gsub("^= ", "", int)) %>%
  pivot_wider(names_from = Variable, values_from = c(Effect, int), names_prefix = "") %>%
  dplyr::select(Test, Treatment, 
                Effect_Choice, int_Choice,
                Effect_Interest, int_Interest) %>%
  mutate(Treatment = factor(Treatment, levels = c("Hot-Control (n = 20)",
                                      "Hot-CORT (n = 20)",
                                      "Cold-Control (n = 20)",
                                      "Cold-CORT (n = 20)"))) %>%
data.frame()
#
#
# Making the table
#
## Table format
set_flextable_defaults(
 font.family = "Arial",
 font.size = 10)
# Create the table
disc_table <- flextable(disc_df_final) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    flextable::compose(i = c(2:4,6:8,10:12,14:16,18:20), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::add_header_row(values = c("", "Choice", "Interest"), colwidths = c(2, 2, 2)) %>%
    flextable::compose(i = 2, j = c(3,5), value = as_paragraph("Median"), part = "head") %>% 
    flextable::compose(i = 2, j = c(4,6), value = as_paragraph("95% HPDI"), part = "head") %>% 
    hline(i = c(4,8,12,16), j = c(1:6), part = "body") %>% # To make some horizontal lines
    vline(i = c(1:20), j = 2, part = "body") %>% # To make some vertical lines
    vline(i = c(1:20), j = 4, part = "body") %>% # To make some vertical lines
    vline(i = c(1,2), j = 2, part = "head") %>%
    vline(i = c(1,2), j = 4, part = "head") 
disc_table
#
#
#
#
#
#| label: tbl-discraw
#| tbl-cap: "Performance of each treatment in each of the numerical tests using the raw data. For the variable Latency and Interest, we show the median and the 95% CI. For the varaible Choice, we show the proportion of individuals that chose the higher number of crickets first."
#
data_table_raw <- data_num %>%
  dplyr::select(trt, test_type, latency, choice, compared_interest) %>%
  group_by(trt, test_type) %>%
  summarize(
    mean_latency = format_dec(mean(latency, na.rm = TRUE), 2),
    ci_latency_lower = format_dec(quantile(latency, 0.025, na.rm = TRUE), 2),
    ci_latency_upper = format_dec(quantile(latency, 0.975, na.rm = TRUE), 2),
    mean_interest = format_dec(mean(compared_interest, na.rm = TRUE), 2),
    ci_interest_lower = format_dec(quantile(compared_interest, 0.025, na.rm = TRUE), 2),
    ci_interest_upper = format_dec(quantile(compared_interest, 0.975, na.rm = TRUE), 2),
    Choice = as.character(sum(choice == 1, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    ci_latency = paste0("[", ci_latency_lower, " , ", ci_latency_upper, "]"),
    ci_interest = paste0("[", ci_interest_lower, " , ", ci_interest_upper, "]"),
  ) %>%
  mutate(
    Latency = trimws(paste0(mean_latency, "\n", ci_latency)),
    Interest = trimws(paste0(mean_interest, "\n", ci_interest))
  ) %>%
  mutate(trt = recode(trt, 
            "B_23" = "Cold-CORT (n = 20)", 
            "A_23" = "Cold-Control (n = 20)", 
            "B_28" = "Hot-CORT (n = 20)", 
            "A_28" = "Hot-Control (n = 20)")) %>% # To order the treatments
  dplyr::select(test_type, trt, Latency, Choice, Interest) %>% 
  pivot_longer(cols = c(Latency, Choice, Interest), 
               names_to = "variable", 
               values_to = "value") %>%
  pivot_wider(names_from = test_type, 
              values_from = value) %>%
  dplyr::select(variable, trt, "1 VS 4", "1 VS 3", "2 VS 4", "2 VS 3", "3 VS 4") %>%
  mutate(variable = factor(variable, levels = c("Latency", "Choice", "Interest")),
         trt = factor(trt, levels = c("Hot-Control (n = 20)",
                                      "Hot-CORT (n = 20)",
                                      "Cold-Control (n = 20)",
                                      "Cold-CORT (n = 20)"))) %>%
  arrange(variable, trt) %>%
data.frame()
data_table_raw_b <- data_table_raw %>%
  rename(
    Variable = variable, 
    Treatment = trt, 
    "1 VS 4" = X1.VS.4,
    "1 VS 3" = X1.VS.3,
    "2 VS 4" = X2.VS.4,
    "2 VS 3" = X2.VS.3,
    "3 VS 4" = X3.VS.4) 
#
#
# Table format
set_flextable_defaults(
 font.family = "Arial",
 font.size = 10)
# Create the table
table_raw <- flextable(data_table_raw_b) %>%
  # Add your custom header rows
  flextable::add_header_row(values = c("", "Tests", ""), colwidths = c(4, 1, 2)) %>%
  # Take away some text
  flextable::compose(i = c(2,3,4,6,7,8,10,11,12), j = 1, value = as_paragraph(""), part = "body") %>%
  hline(i = c(4,8), j = c(1:7), part = "body") # To make some horizontal lines
table_raw
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
#| label: summary_org
#
## Get the summary of the posteriors to get the estimates of each parameter of the model
sum <- summary(posteriors)
## Extract the values of the summary
sum_2 <- as.data.frame(sum) %>%
  dplyr::select(-mad) %>%
  filter(!str_detect(variable, "^r_"))
sum_df <- sum_2 %>%
  mutate(across(where(is.numeric), ~ sprintf(paste0("%.", 2, "f"), .)))  
#
#
#
#
#| label: tbl-summary_lat
#| tbl-cap: "Summary of the model fitted for the latency variable."
#
source(here("R", "func.R"))
## Get df only for latency
df_sum_lat <- sum_table_org(sum_df, "lat")
## Make table
table_sum_lat <- flextable(df_sum_lat)
table_sum_lat
#
#
#
#
#
#| label: tbl-summary_choice
#| tbl-cap: "Summary of the model fitted for the choice variable."
#
source(here("R", "func.R"))
## Get df only for choice
df_sum_choice <- sum_table_org(sum_df, "choice")
## Make table
table_sum_choice <- flextable(df_sum_choice)
table_sum_choice
#
#
#
#
#
#| label: tbl-summary_int
#| tbl-cap: "Summary of the model fitted for the interest variable."
#
source(here("R", "func.R"))
## Get df only for interest
df_sum_int <- sum_table_org(sum_df, "int")
## Make table
table_sum_int <- flextable(df_sum_int)
table_sum_int
#
#
#
#
#
#| label: tbl-summary_other
#| tbl-cap: "Summary of the model fitted for other effects."
#
source(here("R", "func.R"))
## Get df only for interest
df_sum_other <- sum_df %>%
  filter(!str_detect(variable, "^b_"))
## Make table
table_sum_other <- flextable(df_sum_other)
table_sum_other
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
#| tbl-cap: "Number of individuals per treatment that chose the right (R) or left (L) side in each of the numerical tests. p-value indicates the result of the binomial test comparing the number of choices between sides."
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
#| label: tbl-orientation
#| tbl-cap: "Number of individuals per treatment that chose the cricket oriented horizontally (Horizontal) or vertically (Vertical) in our control tests. p-value indicates the result of the binomial test comparing the number of choices between both choices."
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
            "B_23" = "Cold-CORT (n = 20)", 
            "A_23" = "Cold-Control (n = 20)", 
            "B_28" = "Hot-CORT (n = 20)", 
            "A_28" = "Hot-Control (n = 20)"),
          trt = factor(trt, levels = c("Hot-Control (n = 20)",
                                  "Hot-CORT (n = 20)",
                                  "Cold-Control (n = 20)",
                                  "Cold-CORT (n = 20)"))) %>% # To order the treatments
  dplyr::select(trt, Horizontal_Count, Vertical_Count, p_value) %>% # Remove unused columns
data.frame()
#
# 2) MAKING THE TABLE
#
# Split the table_data df by task
orient_table <- flextable(data_or_res) %>%
  flextable::compose(j = 1, value = as_paragraph(""), part = "header") %>%
  flextable::compose(j = 2, value = as_paragraph("Horizontal"), part = "header") %>%
  flextable::compose(j = 3, value = as_paragraph("Vertical"), part = "header") %>%
  flextable::compose(j = 4, value = as_paragraph("p-value"), part = "header") %>%
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
#
#
#
#
#| label: fig-searchWOS
#| fig.cap: "Search query on Web of Science."

knitr::include_graphics("./Others/Web of Science Search results.png")

#
#
#
#| label: fig-searchScopus
#| fig.cap: "Search query on Scopus."

knitr::include_graphics("./Others/Scopus - Document search results.png")

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
