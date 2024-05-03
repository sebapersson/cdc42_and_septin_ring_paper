library(tidyverse)
library(stringr)
library(ggthemes)
library(readxl)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_theme <- theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                                plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=21))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


# Function that converts a string number, where the dot is replaced by 
# P, to a numeric double. 
# Args:
#     number_str, a number on string format where the . is replaced by P
# Returns:
#     the input converted to numeric 
str_to_number <- function(number_str)
{
  number_str <- str_replace(number_str, "P", ".")
  return(as.numeric(number_str))
}


# Function that will calculate the scale-factor for a parameter 
# given the directory tag. E.g if provided with fn_bzeroR4P05 and 
# asked to extract R scale the function return 4.05, if provided 
# with k2 it will return 1 (as the standard value is used for k2). 
# Args:
#     param c("k2", "k_2", "R"), the parameter to extract
#     param_tag, the parameter tag for the experiment 
# Returns:
#    the numberic scale factor or 1 if standa rd value 
calc_scale_fact_param <- function(param, param_tag)
{
  # Standard value provided 
  if(is.na(param_tag)) return(1)
  
  param_tag <- str_replace(param_tag, param, "")
  scale_factor <- str_to_number(param_tag)
  return(scale_factor)
}


process_files <- function(dir_experiment, file_name="/Pole_data.csv") 
{

  # Process all the directories 
  dirs_exp <- list.files(dir_experiment)
  
  # Patterns for finding scaling for parameters 
  pattern_k2b <- "k2b(\\d|P)+"
  pattern_k3 <- "k3(\\d|P)+"
  pattern_Dm <- "Dm(\\d|P)+"
  pattern_Bem <- "BemGefTot(\\d|P)+"
  pattern_k8h <- "k8h(\\d|P)+"
  pattern_k9h <- "k9h(\\d|P)+"
  pattern_k8max <- "k8max(\\d|P)+"
  pattern_k8max2 <- "k8max2(\\d|P)+"
  pattern_k8max_2 <- "k8max_2(\\d|P)+"
  pattern_cdc42 <- "Cdc42Tot(\\d|P)+"
  pattern_r <- "r(\\d|P)+"
  pattern_f <- "F(\\d|P)+"
  data_tot <- tibble()
  data_time <- tibble()
  for(param_case in dirs_exp){
    # Compute scaling factor
    print(param_case)
    param_tag_k2b <- str_extract(param_case, pattern_k2b)
    scale_fact_k2b <- calc_scale_fact_param("k2b", param_tag_k2b)
    
    param_tag_k8h <- str_extract(param_case, pattern_k8h)
    scale_fact_k8h <- calc_scale_fact_param("k8h", param_tag_k8h)
    
    param_tag_k9h <- str_extract(param_case, pattern_k9h)
    scale_fact_k9h <- calc_scale_fact_param("k9h", param_tag_k9h)    
    
    param_tag_k8max <- str_extract(param_case, pattern_k8max)
    scale_fact_k8max <- calc_scale_fact_param("k8max", param_tag_k8max)
    
    param_tag_k8max2 <- str_extract(param_case, pattern_k8max2)
    scale_fact_k8max2 <- calc_scale_fact_param("k8max", param_tag_k8max2)    
    
    param_tag_k8max_2 <- str_extract(param_case, pattern_k8max_2)
    scale_fact_k8max_2 <- calc_scale_fact_param("k8max_2", param_tag_k8max_2)    
    
    param_tag_k3 <- str_extract(param_case, pattern_k3)
    scale_fact_k3 <- calc_scale_fact_param("k3", param_tag_k3)
    
    param_tag_Dm <- str_extract(param_case, pattern_Dm)
    scale_fact_Dm <- calc_scale_fact_param("Dm", param_tag_Dm)
    
    param_tag_Bem <- str_extract(param_case, pattern_Bem)
    scale_fact_Bem <- calc_scale_fact_param("BemGefTot", param_tag_Bem)
    
    param_tag_r <- str_extract(param_case, pattern_r)
    scale_fact_r <- calc_scale_fact_param("r", param_tag_r)
    
    param_tag_f <- str_extract(param_case, pattern_f)
    scale_fact_f <- calc_scale_fact_param("F", param_tag_f)
    
    param_tag_cdc42 <- str_extract(param_case, pattern_cdc42)
    scale_fact_cdc42 <- calc_scale_fact_param("Cdc42Tot", param_tag_cdc42)
    
    # Process file if exists 
    file_path <- file.path(dir_experiment, param_case, file_name)
    if(file.exists(file_path)){
      # For making it real time 
      data_read <- read_csv(file_path, col_types = cols()) %>%
        mutate(scale_k2b = scale_fact_k2b) %>%
        mutate(scale_bem = scale_fact_Bem) %>%
        mutate(scale_Dm = scale_fact_Dm) %>%
        mutate(scale_r = scale_fact_r) %>%
        mutate(scale_cdc42 = scale_fact_cdc42) %>% 
        mutate(scale_f = scale_fact_f) %>%
        mutate(scale_k3 = scale_fact_k3) |> 
        mutate(scale_k8h = scale_fact_k8h) |> 
        mutate(scale_k9h = scale_fact_k9h) |> 
        mutate(scale_k8max = scale_fact_k8max) |> 
        mutate(scale_k8max2 = scale_fact_k8max2) |> 
        mutate(scale_k8max_2 = scale_fact_k8max_2)
      
    }else{
      next
    }
    data_tot <- data_tot %>% bind_rows(data_read)
  }
  
  
  # Post-process master plotting object 
  if(file_name != "/Pole_dist_end.csv"){
    data_tot <- data_tot %>%
      mutate(scale_Dm = as.factor(scale_Dm)) %>%
      mutate(scale_k2b = as.factor(scale_k2b)) %>%
      mutate(scale_bem = as.factor(scale_bem)) %>%
      mutate(scale_cdc42 = as.factor(scale_cdc42)) %>%
      mutate(scale_k3 = as.factor(scale_k3)) %>%
      mutate(r_num = scale_r * 2.5) %>%
      mutate(pole_size = pole_ratio * r_num^2 * 4 * pi) %>%
      mutate(scale_r = as.factor(scale_r)) %>%
      mutate(scale_f = as.factor(scale_f)) %>%
      mutate(scale_k8h = as.factor(scale_k8h)) |> 
      mutate(scale_k9h = as.factor(scale_k9h)) |>       
      mutate(scale_k8max = as.factor(scale_k8max)) |> 
      mutate(scale_k8max2 = as.factor(scale_k8max2)) |> 
      mutate(scale_k8max_2 = scale_fact_k8max_2) |> 
      mutate(term_criteria = as.factor(term_crit))
  }else{
      data_tot <- data_tot %>%
        mutate(scale_Dm = as.factor(scale_Dm)) %>%
        mutate(scale_k2b = as.factor(scale_k2b)) %>%
        mutate(scale_bem = as.factor(scale_bem)) %>%
        mutate(scale_k3 = as.factor(scale_k3)) %>%
        mutate(r_num = scale_r * 2.5) %>%
        mutate(scale_k8h = as.factor(scale_k8h)) |> 
        mutate(scale_k9h = as.factor(scale_k9h)) |>         
        mutate(scale_k8max = as.factor(scale_k8max)) |> 
        mutate(scale_k8max2 = as.factor(scale_k8max2)) |> 
        mutate(scale_k8max_2 = scale_fact_k8max_2) |> 
        mutate(scale_r = as.factor(scale_r)) 
  }
  
  return(data_tot)
}


# We extract at max conc.
get_cdc42max <- function(dir_experiment) 
{

  # Process all the directories 
  dirs_exp <- list.files(dir_experiment)
  
  # Patterns for finding scaling for parameters 
  pattern_k2b <- "k2b(\\d|P)+"
  pattern_k3 <- "k3(\\d|P)+"
  pattern_Dm <- "Dm(\\d|P)+"
  pattern_Bem <- "BemGefTot(\\d|P)+"
  pattern_k8h <- "k8h(\\d|P)+"
  pattern_k9h <- "k9h(\\d|P)+"
  pattern_k8max <- "k8max(\\d|P)+"
  pattern_k8max2 <- "k8max2(\\d|P)+"
  pattern_k8max_2 <- "k8max_2(\\d|P)+"
  pattern_cdc42 <- "Cdc42Tot(\\d|P)+"
  pattern_r <- "r(\\d|P)+"
  pattern_f <- "F(\\d|P)+"
  data_tot <- tibble()
  for(param_case in dirs_exp){
    # Compute scaling factor
    print(param_case)
    param_tag_k2b <- str_extract(param_case, pattern_k2b)
    scale_fact_k2b <- calc_scale_fact_param("k2b", param_tag_k2b)
    
    param_tag_k8h <- str_extract(param_case, pattern_k8h)
    scale_fact_k8h <- calc_scale_fact_param("k8h", param_tag_k8h)
    
    param_tag_k9h <- str_extract(param_case, pattern_k9h)
    scale_fact_k9h <- calc_scale_fact_param("k9h", param_tag_k9h)    
    
    param_tag_k8max <- str_extract(param_case, pattern_k8max)
    scale_fact_k8max <- calc_scale_fact_param("k8max", param_tag_k8max)
    
    param_tag_k8max2 <- str_extract(param_case, pattern_k8max2)
    scale_fact_k8max2 <- calc_scale_fact_param("k8max", param_tag_k8max2)    
    
    param_tag_k8max_2 <- str_extract(param_case, pattern_k8max_2)
    scale_fact_k8max_2 <- calc_scale_fact_param("k8max_2", param_tag_k8max_2)    
    
    param_tag_k3 <- str_extract(param_case, pattern_k3)
    scale_fact_k3 <- calc_scale_fact_param("k3", param_tag_k3)
    
    param_tag_Dm <- str_extract(param_case, pattern_Dm)
    scale_fact_Dm <- calc_scale_fact_param("Dm", param_tag_Dm)
    
    param_tag_Bem <- str_extract(param_case, pattern_Bem)
    scale_fact_Bem <- calc_scale_fact_param("BemGefTot", param_tag_Bem)
    
    param_tag_r <- str_extract(param_case, pattern_r)
    scale_fact_r <- calc_scale_fact_param("r", param_tag_r)
    
    param_tag_f <- str_extract(param_case, pattern_f)
    scale_fact_f <- calc_scale_fact_param("F", param_tag_f)
    
    param_tag_cdc42 <- str_extract(param_case, pattern_cdc42)
    scale_fact_cdc42 <- calc_scale_fact_param("Cdc42Tot", param_tag_cdc42)
    
    # Process file if exists 
    file_path <- file.path(dir_experiment, param_case, "Pole_dist_t.csv")
    if(file.exists(file_path)){
      # For making it real time 
      data_read <- read_csv(file_path, col_types = cols())
      for(i in unique(data_read$index)){
        data_index <- data_read |> 
          filter(index == i) 
        tmax = data_index$t[which.max(data_index$Cdc42T)]
        data_tmax <- data_index |> filter(abs(t - tmax) < 1e-8)
        data_max_val <- data_tmax |> 
          group_by(t) |> 
          summarise(max_val = max(Cdc42T + BemGEF42), 
                    cdc42max = max(Cdc42T))
        data_ret <- data_max_val |> 
          select(max_val, cdc42max) |> 
          mutate(scale_k2b = factor(scale_fact_k2b)) %>%
          mutate(scale_bem = factor(scale_fact_Bem)) %>%
          mutate(scale_Dm = scale_fact_Dm) %>%
          mutate(scale_r = scale_fact_r) %>%
          mutate(scale_cdc42 = scale_fact_cdc42) %>% 
          mutate(scale_f = scale_fact_f) %>%
          mutate(scale_k3 = scale_fact_k3) |> 
          mutate(scale_k8h = scale_fact_k8h) |> 
          mutate(scale_k9h = scale_fact_k9h) |> 
          mutate(scale_k8max = scale_fact_k8max) |> 
          mutate(scale_k8max2 = scale_fact_k8max2) |> 
          mutate(scale_k8max_2 = scale_fact_k8max_2) |> 
          mutate(index = i) |> 
          mutate(r_num = as.numeric(scale_r) * 2.5)
        data_tot <- bind_rows(data_tot, data_ret)
      }
    }else{
      next
    }
  }
  return(data_tot)
}

