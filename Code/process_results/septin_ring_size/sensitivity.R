library(tidyverse)
library(stringr)
library(readxl)

source("ring_size_functions.R")

# General plotting parameters
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_theme <- theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=21))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0

# -----------------------------------------------------------------------------------------------------------
# Functions for helping to process the results
# -----------------------------------------------------------------------------------------------------------
str_to_number <- function(number_str)
{
  number_str <- str_replace(number_str, "P", ".")
  return(as.numeric(number_str))
}

calc_scale_fact_param <- function(param, param_tag)
{
  # Standard value provided 
  if(is.na(param_tag)) return(1)
  
  param_tag <- str_replace(param_tag, param, "")
  scale_factor <- str_to_number(param_tag)
  return(scale_factor)
}

process_file_name <- function(file_name) 
{
  # Patterns for finding scaling for parameters 
  pattern_k12b <- "k12b(\\d|P)+"
  pattern_k13 <- "k13(\\d|P)+"
  pattern_k15 <- "k15(\\d|P)+"
  pattern_k16 <- "k16(\\d|P)+"
  pattern_k17 <- "k17(\\d|P)+"
  pattern_k18 <- "k18(\\d|P)+"
  pattern_k19 <- "k19(\\d|P)+"
  pattern_k20 <- "k20(\\d|P)+"
  pattern_k21 <- "k21(\\d|P)+"
  pattern_k22 <- "k22(\\d|P)+"
  pattern_k23 <- "k23(\\d|P)+"
  pattern_k24 <- "k24(\\d|P)+"
  pattern_k25 <- "k25(\\d|P)+"
  
  param_tag_k12b <- str_extract(file_name, pattern_k12b)
  scale_fact_k12b <- calc_scale_fact_param("k12b", param_tag_k12b)
  
  param_tag_k13 <- str_extract(file_name, pattern_k13)
  scale_fact_k13 <- calc_scale_fact_param("k13", param_tag_k13)
  
  param_tag_k15 <- str_extract(file_name, pattern_k15)
  scale_fact_k15 <- calc_scale_fact_param("k15", param_tag_k15)
  
  param_tag_k16 <- str_extract(file_name, pattern_k16)
  scale_fact_k16 <- calc_scale_fact_param("k16", param_tag_k16)
  
  param_tag_k17 <- str_extract(file_name, pattern_k17)
  scale_fact_k17 <- calc_scale_fact_param("k17", param_tag_k17)
  
  param_tag_k18 <- str_extract(file_name, pattern_k18)
  scale_fact_k18 <- calc_scale_fact_param("k18", param_tag_k18)
  
  param_tag_k18 <- str_extract(file_name, pattern_k18)
  scale_fact_k18 <- calc_scale_fact_param("k18", param_tag_k18)
  
  param_tag_k19 <- str_extract(file_name, pattern_k19)
  scale_fact_k19 <- calc_scale_fact_param("k19", param_tag_k19)
  
  param_tag_k20 <- str_extract(file_name, pattern_k20)
  scale_fact_k20 <- calc_scale_fact_param("k20", param_tag_k20)
  
  param_tag_k21 <- str_extract(file_name, pattern_k21)
  scale_fact_k21 <- calc_scale_fact_param("k21", param_tag_k21)
  
  param_tag_k22 <- str_extract(file_name, pattern_k22)
  scale_fact_k22 <- calc_scale_fact_param("k22", param_tag_k22)
  
  param_tag_k23 <- str_extract(file_name, pattern_k23)
  scale_fact_k23 <- calc_scale_fact_param("k23", param_tag_k23)
  
  param_tag_k24 <- str_extract(file_name, pattern_k24)
  scale_fact_k24 <- calc_scale_fact_param("k24", param_tag_k24)
  
  param_tag_k25 <- str_extract(file_name, pattern_k25)
  scale_fact_k25 <- calc_scale_fact_param("k25", param_tag_k25)
  
  data_scale <- tibble(k12b = scale_fact_k12b, k13 = scale_fact_k13, k15 = scale_fact_k15,
                       k16 = scale_fact_k16, k17 = scale_fact_k17, k18 = scale_fact_k18, 
                       k19 = scale_fact_k19, k20 = scale_fact_k20, k21 = scale_fact_k21,
                       k22 = scale_fact_k22, k23 = scale_fact_k23, k24 = scale_fact_k24,
                       k25 = scale_fact_k25)
  return(data_scale)
}

get_p_change <- function(ref_experiment, test_experiment)
{
  data_ref <- process_file_name(ref_experiment)
  data_test <- process_file_name(test_experiment)
  data_change <- tibble()
  k <- 0
  for(i in 1:ncol(data_ref)){
    pname <- colnames(data_ref)[i]
    val_ref <- data_ref[1, i] |> as.numeric()
    val_test <- data_test[1, i] |> as.numeric()
    if(val_ref != val_test){
      data_change <- data_change |> 
        bind_rows(tibble(param = pname, 
                         fact = val_test / val_ref))
      k <- k + 1
    }
  }
  if(k != 1){
    print("Something wrong with parameter names")
  }
  return(data_change)
}

# -----------------------------------------------------------------------------------------------------------
# Process for wide
# -----------------------------------------------------------------------------------------------------------
ref_experiment <- "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding"
dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Sense_wide_standard")

# Processing ring diameter takes quite some time, hence results are written to disk in order to avoid 
# recomputing
path_save_tmp <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Wide_sense.csv")
if(!file.exists(path_save_tmp)){
  test_experiments <- list.files(dir_results)
  data_res <- tibble()
  for(test_experiment in test_experiments){
    p_change <- get_p_change(ref_experiment, test_experiment)
    if(!(p_change$fact[1] %in% c(0.5, 2.0))) next
    print(sprintf("Parameter %s with factor %.1f", p_change$param[1], p_change$fact[1]))
    
    ring_diameter <- process_no_cable(file.path(dir_results, test_experiment), c(1, 2, 3), "0.4") |> 
      group_by(tag, index) |> 
      summarise(median_inner = mean(d_min), 
                median_outer = mean(d_max), 
                median_mean = mean((d_min+d_max)*0.5)) |> 
      mutate(parameter = p_change$param[1],
             parameter_factor = p_change$fact[1])
    data_res <- bind_rows(data_res, ring_diameter)
  }
  write_csv(data_res, path_save_tmp)
  data_res_wide <- data_res |> 
    mutate(tag = "wide")

}else{
  data_res_wide <- read_csv(path_save_tmp, col_types = cols()) |> 
    mutate(tag = "wide")
}

# -----------------------------------------------------------------------------------------------------------
# Process for not-wide
# -----------------------------------------------------------------------------------------------------------
ref_experiment <- "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding"
dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Sense_conc_standard")

path_save_tmp <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Not_ide_sense.csv")
if(!file.exists(path_save_tmp)){
  test_experiments <- list.files(dir_results)
  data_res <- tibble()
  for(test_experiment in test_experiments){
    p_change <- get_p_change(ref_experiment, test_experiment)
    
    # Combinations where simulations fail. As this is a sensitivity analysis 
    # not every parameter combination is expected to be robust
    if(p_change$fact[1] == 2.0 && p_change$param == "k15") next
    
    if(!(p_change$fact[1] %in% c(0.5, 2.0))) next
    print(sprintf("Parameter %s with factor %.1f", p_change$param[1], p_change$fact[1]))
    
    ring_diameter <- process_no_cable(file.path(dir_results, test_experiment), c(1, 2, 3), "0.4") |> 
      group_by(tag, index) |> 
      summarise(median_inner = mean(d_min), 
                median_outer = mean(d_max), 
                median_mean = mean((d_min+d_max)*0.5)) |> 
      mutate(parameter = p_change$param[1],
             parameter_factor = p_change$fact[1])
    data_res <- bind_rows(data_res, ring_diameter)
  }
  write_csv(data_res, path_save_tmp)
  data_res_not_wide <- data_res |> 
    mutate(tag = "not_wide")

}else{
  data_res_not_wide <- read_csv(path_save_tmp, col_types = cols()) |> 
    mutate(tag = "not_wide")
}

# -----------------------------------------------------------------------------------------------------------
# Putting it all together in a way I can plot the data
# -----------------------------------------------------------------------------------------------------------
## Factor 0.5
factor_check <- 0.5
parameters_check <- c("k12b", "k13", "k15", "k16", "k17", "k18", "k19", "k20", "k21", 
                      "k22", "k23", "k24", "k25")
data_plot <- tibble()
for(parameter_check in parameters_check){
  data_wide <- data_res_wide |> 
    filter(parameter_factor == factor_check & parameter == parameter_check)
  data_not_wide <- data_res_not_wide |> 
    filter(parameter_factor == factor_check & parameter == parameter_check)
  mean_not_wide <- mean(data_not_wide$median_mean)
  data_plot_tmp <- bind_rows(data_wide, data_not_wide) |> 
    mutate(diameter_norm = median_mean /mean_not_wide)
  data_plot <- bind_rows(data_plot, data_plot_tmp)
}
data_plot <- data_plot |> mutate(parameter = factor(parameter))

p1 <- ggplot(data_plot, aes(tag, diameter_norm)) + 
  geom_jitter(aes(group = parameter_check), width = 0.1) + 
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", linewidth= 0.3, geom = "crossbar") + 
  facet_wrap(vars(parameter), nrow = 3, ncol = 5) + 
  scale_x_discrete(labels = c("Conentrated exocytosis", "Wide exocytosis")) + 
  labs(x = "", y = "Normalised septin ring diameter", 
       title = "Parameter decrease with factor 0.5") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

## Factor 2.0
factor_check <- 2.0
parameters_check <- c("k12b", "k13", "k16", "k17", "k18", "k19", "k20", "k21", 
                      "k22", "k23", "k24", "k25")
data_plot <- tibble()
for(parameter_check in parameters_check){
  data_wide <- data_res_wide |> 
    filter(parameter_factor == factor_check & parameter == parameter_check)
  data_not_wide <- data_res_not_wide |> 
    filter(parameter_factor == factor_check & parameter == parameter_check)
  mean_not_wide <- mean(data_not_wide$median_mean)
  data_plot_tmp <- bind_rows(data_wide, data_not_wide) |> 
    mutate(diameter_norm = median_mean /mean_not_wide)
  data_plot <- bind_rows(data_plot, data_plot_tmp)
}
data_plot <- data_plot |> mutate(parameter = factor(parameter))

p2 <- ggplot(data_plot, aes(tag, diameter_norm)) + 
  geom_jitter(aes(group = parameter_check), width = 0.1) + 
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", linewidth= 0.3, geom = "crossbar") + 
  facet_wrap(vars(parameter), nrow = 3, ncol = 5) + 
  scale_x_discrete(labels = c("Conentrated exocytosis", "Wide exocytosis")) + 
  labs(x = "", y = "Normalised septin ring diameter", 
       title = "Parameter increase with factor 2.0") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

dir_save <- file.path("..", "..", "..", "Results", "Septin_ring_v2", "Sensitivity")
if(!dir.exists(dir_save)) dir.create(dir_save)
ggsave(file.path(dir_save, "factor_0P5.svg"), p1, width = BASE_WIDTH*3, height = BASE_HEIGHT*2)
ggsave(file.path(dir_save, "factor_2.svg"), p2, width = BASE_WIDTH*3, height = BASE_HEIGHT*2)
