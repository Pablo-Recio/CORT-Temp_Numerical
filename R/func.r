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
