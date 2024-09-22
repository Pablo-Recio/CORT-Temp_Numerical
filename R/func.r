#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)
#
####################
####################  
# Get predictions based on posteriors
#' @title org_post
#' @param var The response variable we are interested in ("lat", "choice", "int")
org_post <- function(var){
  df <- data.frame() #create empty df
  sex <- unique(data_num$sex)
  test <- data_num$test_type %>%
    unique() %>%
    levels() %>%
    gsub(" ", "", .) %>%
    factor()
  therm <- unique(data_num$temp)
  horm <- unique(data_num$cort) # create the vectors employed for the loops
  # create names for each variable and each df:
  if(var == "lat"){
    name <- "loglatency"
  } else if(var == "choice"){
    name <- "choice"
  } else if(var == "int"){
    name <- "comparedinterest"
  } else{
    stop("Variable not valid")
  }
  # Modify the main df depending on the selected variable
  data <- posteriors %>% 
    dplyr::select(contains(paste0("b_", name))) %>%
    dplyr::select(-matches("age"))
  # Create function to select the columns
  choose_col <- function(initial_df, test_type, temperature, hormone){
    if(test_type == "1VS4"){
      col_df_1 <- initial_df %>%
        dplyr::select(-dplyr::matches("1VS3|2VS4|2VS3|3VS4"))
    } else if (test_type == "1VS3"){
      col_df_1 <- initial_df %>%
        dplyr::select(-dplyr::matches("2VS4|2VS3|3VS4"))
    } else if(test_type == "2VS4"){
      col_df_1 <- initial_df %>%
        dplyr::select(-dplyr::matches("1VS3|2VS3|3VS4"))
    } else if(test_type == "2VS3"){
      col_df_1 <- initial_df %>%
        dplyr::select(-dplyr::matches("1VS3|2VS4|3VS4"))
    } else if(test_type == "3VS4"){
      col_df_1 <- initial_df %>%
        dplyr::select(-dplyr::matches("1VS3|2VS4|2VS3"))
    }
    if(temperature == "Hot"){
      col_df_2 <- col_df_1
    } else if(temperature == "Cold"){
      col_df_2 <- col_df_1 %>%
        dplyr::select(-dplyr::matches("Hot")) %>%
      data.frame()
    }
    if(hormone == "Control"){
      col_df_final <- col_df_2
    } else if(hormone == "CORT"){
      col_df_final <- col_df_2 %>%
        dplyr::select(-dplyr::matches("Control")) %>%
      data.frame()
    }
  return(col_df_final)
  }
  # Create empty lists for females and males
  males <- list()
  females <- list()
  for(i in sex){
    if(i == "m"){
      df_m <- data %>%
        mutate(!!paste0("b_", name, "_Intercept") := 
                 !!sym(paste0("b_", name, "_Intercept")) + 
                 !!sym(paste0("b_", name, "_sexm"))) %>%
        dplyr::select(-contains(paste0("b_", name, "_sexm")))
      # Loop over test types, temperatures, and hormone levels
      for(j in test){
        for(k in therm){
          for(l in horm){
            vector_m_df <- choose_col(df_m, j, k, l)
            vector_m <- rowSums(vector_m_df)
            sum_vector_m <- base::sample(vector_m, size = 6000, replace = FALSE)
            column_name_m <- paste0(j, "_", k, "_", l)
            assign(column_name_m, sum_vector_m)
            males[[column_name_m]] <- get(column_name_m)
          }
        }
      }
    } else if(i == "f"){
      df_f <- data %>% 
        dplyr::select(-contains(paste0("b_", name, "_sexm")))
      # Loop over test types, temperatures, and hormone levels
      for(j in test){
        for(k in therm){
          for(l in horm){
            vector_f_df <- choose_col(df_f, j, k, l)
            vector_f <- rowSums(vector_f_df)
            sum_vector_f <- base::sample(vector_f, size = 6000, replace = FALSE)
            column_name_f <- paste0(j, "_", k, "_", l)
            assign(column_name_f, sum_vector_f)
            females[[column_name_f]] <- get(column_name_f)
          }
        }
      }
    }
  }
# Convert lists to data frames
males_df <- as.data.frame(males)
females_df <- as.data.frame(females)
# Get final df
df <- rbind(males_df, females_df) 
return(df)
}
####################
####################
# Estimate p-values using pmcm
#' @title pMCMC Function
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
#' @param dir The direction of the one-tail test (<) or (>)
pmcmc <- function(x, null = 0, twotail = TRUE, dir){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    if(dir == "<"){
    (max(sum(x>=null) / length(x)))
    } else{
      if(dir == ">"){
        (max(sum(x<=null) / length(x)))
      } else{
        stop("dir not valid")
      }
    }
  }
}
####################
####################
# Function to format numbers with 2 decimal places
#' @title format_dec
#' @param x The object
#' @param n The number of decimals
format_dec <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  return(as.numeric(z))
}
####################
####################
# Function to format p_values with n decimal places
#' @title format_p
#' @param x The object
#' @param n The number of decimals
format_p <- function(x, n) {
  z <- sprintf(paste0("%.",n,"f"), x)
  tmp <- ifelse(as.numeric(z) <= 0.001, "< 0.001",
         ifelse(as.numeric(z) <= 0.05 & as.numeric(z) > 0.001, "< 0.05",
                paste0("= ", as.character(z))))
  return(tmp)
}
####################
####################
# Function to organise the posteriors for the summary (Suppl) tables
#' @title sum_table_org
#' @param data_post the df employed
#' @param beh the variable selected for each table ("lat", "choice", "int")
sum_table_org <- function(data_post, beh){
  if(beh == "lat"){
    name_var <- "b_loglatency"
  } else if(beh == "choice"){
    name_var <- "b_choice"
  } else if(beh == "int"){
    name_var <- "b_comparedinterest"
  } else {
    stop("Incorrect variable")
  }
  data_sum <- data_post %>%
    filter(str_detect(variable, name_var))
return(data_sum)
}
####################
####################
# Function to create the df for the violin plots
#' @title df_violin_plot
#' @param data_fig for choosing df based on the variable selected ("lat", "choice", "int")
df_violin_plot <- function(data_fig){
  test <- c("1VS4", "1VS3", "2VS4", "2VS3", "3VS4")
  fig_df_1 <- data.frame()
  if(data_fig == "lat"){
    fig_df_0 <- post_lat
  } else if(data_fig == "choice"){
    fig_df_0 <- prob_choice
  } else if(data_fig == "int"){
    fig_df_0 <- post_int
  } else {
    stop("Incorrect data")
  }
  for(m in test){
    # Modify the df to get the values per each test
    fig_df_sub <- fig_df_0 %>%
      dplyr::select(dplyr::matches(m))
    # Extract the required columns
    Cold_CORT <- fig_df_sub[[paste0("X", m, "_Cold_CORT")]]
    Cold_Control <- fig_df_sub[[paste0("X", m, "_Cold_Control")]]
    Hot_CORT <- fig_df_sub[[paste0("X", m, "_Hot_CORT")]]
    Hot_Control <- fig_df_sub[[paste0("X", m, "_Hot_Control")]]
    # Create a temporary data frame for this loop iteration
    temporal_df <- data.frame(cbind(Test = m,
                                    Cold_CORT,
                                    Cold_Control,
                                    Hot_CORT,
                                    Hot_Control))
    # Merge the temporary data frame to the main df_tables
    fig_df_1 <- rbind(fig_df_1, temporal_df)
  }
  # Modify the obtain df to make it adequate for the plot
  df_viol_plot <- fig_df_1 %>%
    pivot_longer(cols = c(Cold_CORT, Cold_Control, Hot_CORT, Hot_Control),
               names_to = "Treatment",
               values_to = "Estimate") %>%
    mutate(Treatment = factor(Treatment, levels = c("Cold_CORT", "Cold_Control", "Hot_CORT", "Hot_Control"), 
                                      labels = c("Cold_CORT" = "Cold-CORT (n = 20)",
                                                  "Cold_Control" = "Cold-Control (n = 20)",
                                                  "Hot_CORT" = "Hot-CORT (n = 20)",
                                                  "Hot_Control" = "Hot-Control (n = 20)")),
            Test = factor(Test, levels = c("1VS4", "1VS3", "2VS4", "2VS3", "3VS4")),
            Estimate = as.numeric(Estimate)) %>%
  data.frame()
#
return(df_viol_plot)
}
####################
####################
# Function to create the df for the bar plots
#' @title df_bar_plot
#' @param data_fig for choosing df based on the variable selected ("lat", "choice", "int")
df_bars_plot <- function(data_fig){
  if(data_fig == "lat"){
    fig_df_2 <- df_viol_lat
  } else if(data_fig == "choice"){
    fig_df_2 <- df_viol_choice
  } else if(data_fig == "int"){
    fig_df_2 <- df_viol_int
  } else {
    stop("Incorrect data")
  }
  # Modify the df for the bar plots
  df_bars_plot <- fig_df_2 %>%
    group_by(Test, Treatment) %>%
    summarize(mean_estimate = mean(Estimate, na.rm = TRUE),
              sd_estimate = sd(Estimate, na.rm = TRUE)) %>%
    ungroup()
  data.frame()
#
return(df_bars_plot)
}
####################
####################
# Function to create the figure for the bar plots
#' @title plotting
#' @param type_var for choosing df based on the variable selected ("lat", "choice", "int")
plotting <- function(type_var){
  # Establish the df depending on the variable
  df_violin_name <- paste0("df_viol_", type_var)
  df_violin <- get(df_violin_name)
  #
  df_points_name <- paste0("df_bars_", type_var)
  df_points <- get(df_points_name)
  # Plot the predicted estimates for each treatment and test
  plot <- ggplot(df_violin, aes(x = Treatment, y = Estimate, fill = Treatment)) +
  geom_flat_violin(alpha = 0.6) + 
  scale_fill_manual(values = c("Cold-CORT (n = 20)" = "#00008B",
                              "Cold-Control (n = 20)" = "#68bde1",
                              "Hot-CORT (n = 20)" = "#b50101",
                              "Hot-Control (n = 20)" = "#fa927d"),
                    breaks = c("Hot-Control (n = 20)",
                              "Hot-CORT (n = 20)",
                              "Cold-Control (n = 20)",
                              "Cold-CORT (n = 20)")) +
  facet_wrap(~ Test, scales = "free_x", nrow = 1) +  # Facet by Treatment
  geom_point(data = df_points, aes(y = mean_estimate, x = Treatment), position = position_dodge(width = 0.75), 
              color = "black", fill = "black", size = 3) +
  geom_segment(data = df_points, aes(y = mean_estimate - sd_estimate, yend = mean_estimate + sd_estimate, 
                x = Treatment, xend = Treatment), size = 1.5, color = "black") +
  ylim(min(df_violin$Estimate), max(df_violin$Estimate)) +
  coord_flip() +
  theme_few() +
  theme(plot.margin = margin(2, 3, 5, 3, "mm"),
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.text = element_text(size = 10, family = "sans"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, family = "sans",
                                margin = margin(t = 8)),
    legend.position = "right",
    legend.title = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 11, family = "sans")
    )
  # Specify details depending on variable
  if(type_var == "lat"){
    final_plot <- plot +
    labs(y = "Latency")
  } else if(type_var == "choice"){
    final_plot <- plot +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    scale_y_continuous(
      breaks = c(0.25, 0.5, 0.75),
      labels = c("0.25", "0.5", "0.75")) +  # Specify which labels to show
    labs(y = "Choice")
  } else if(type_var == "int"){
    final_plot <- plot +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(y = "Interest")
  } else {
    stop("Incorrect var_type")
  }
  #
return(final_plot)
}

