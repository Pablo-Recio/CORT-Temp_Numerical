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
  test <- unique(data_num$test_type)
  temp <- unique(data_num$temp)
  horm <- unique(data_num$cort) # create the vectors employed for the loops
  # create names for each variable and each df:
  if(var == "lat"){
    name <- "loglatency"
  } else if(var == "choice"){
    name <- "choice"
  } else if(var == "int"){
    name <- comparedinterest
  } else{
    stop("Variable not valid")
  } 
  data <- posteriors %>% dplyr::select(contains(paste0("b_", name)))
  for(i in sex){
    if(i == "m"){
      df_m <- data %>%
        mutate(!!paste0("b_", name, "_Intercept") := 
                 !!sym(paste0("b_", name, "_Intercept")) + 
                 !!sym(paste0("b_", name, "_sexm"))) %>%
        dplyr::select(-contains(paste0("b_", name, "_sexm")))
    } else if(i == "f"){
      df_f <- data %>% 
        dplyr::select(-contains(paste0("b_", name, "_sexm")))
    }
  }
df <- rbind(df_m, df_f) 
return(df)
}
print(paste0("b_", name, "_Intercept"))
print(paste0("b_", name, "sexm"))
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
