#################### 
####################
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes)
#
####################
####################
# Extract sample size per group
#' @title sample
#' @description Estract sample size per group
#' @param corti To select cort treatment ("CORT"/"Control")
#' @param therm To select temp treatment ("Cold"/"Hot")
#' @param bias To select the side selected first ("Left"/"Right")
sample <- function(corti, therm, bias){
  #Specify database
  data <- data_num
  # Count sample
  sample_size <- data %>%
                filter(cort == corti, temp == therm, side_choice == bias, test_type == "1 VS 4") %>%
                group_by(lizard_id) %>%
                summarise(n = n()) %>%
                summarise(total_count = sum(n)) %>%
                pull(total_count)
  return(sample_size)
}
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
            sum_vector_m <- base::sample(vector_m, size = 2500, replace = FALSE)
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
            sum_vector_f <- base::sample(vector_f, size = 2500, replace = FALSE)
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
