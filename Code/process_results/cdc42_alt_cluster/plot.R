library(tidyverse)
library(stringr)
library(ggthemes)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_theme <- theme_tufte(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                                plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0

# -----------------------------------------------------------------------------------------------------------
# Helper functions to process the results
# -----------------------------------------------------------------------------------------------------------

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


# Function that outputs a data-frame (tibble)
# of the output from a specific experiment. 
# More specifically, the function aggregates 
# the Simulation_data.csv into a single 
# tibble which then can be plotted, or 
# procssed further. 
# Args:
#     dir_experiment, path to where in 
#        intermeidate the result is stored
# Returns:
#     tibble of aggregated Simulated_data.csv 
#        files 
process_simulated_data_files <- function(dir_experiment)
{
  
  # Sanity check input 
  if(!dir.exists(dir_experiment)){
    print("Error: Provided experiment directory does not exist")
    return(1)
  }
  
  # Process all the directories 
  dirs_exp <- list.files(dir_experiment)
  
  # Patterns for finding scaling for parameters 
  pattern_k2 <- "k2(\\d|P)+"
  pattern_k_2 <- "k_2(\\d|P)+"
  pattern_R <- "R(\\d|P)+"
  pattern_G0 <- "G0(\\d|P)+"
  data_tot <- tibble()
  for(param_case in dirs_exp){
    # Compute scaling factors
    param_tag_k2 <- str_extract(param_case, pattern_k2)
    param_tag_k_2 <- str_extract(param_case, pattern_k_2)
    param_tag_R <- str_extract(param_case, pattern_R)
    param_tag_G0 <- str_extract(param_case, pattern_G0)
    scale_fact_k2 <- calc_scale_fact_param("k2", param_tag_k2)
    scale_fact_k_2 <- calc_scale_fact_param("k_2", param_tag_k_2)
    scale_fact_R <- calc_scale_fact_param("R", param_tag_R)
    scale_fact_G0 <- calc_scale_fact_param("G0", param_tag_G0)
    
    # Process file if exists 
    file_path <- file.path(dir_experiment, param_case, "Simulation_data.csv")
    if(file.exists(file_path)){
      # For making it real time 
      DA <- 0.011
      R_val <- 3.95 * scale_fact_R
      time_converstion <- (R_val^2) / (DA * 60)
      
      data_read <- read_csv(file_path, col_types = cols()) %>%
        mutate(scale_R = scale_fact_R) %>%
        mutate(scale_k2 = scale_fact_k2) %>%
        mutate(scale_k_2 = scale_fact_k_2) |> 
        mutate(scale_G0 = scale_fact_G0)
      
      if("X1" %in% colnames(data_read)) data_read <- data_read |> select(-X1)
      
    }else{
      next
    }
    data_tot <- data_tot %>% bind_rows(data_read)
  }
  # Post-process master plotting object 
  data_tot <- data_tot %>%
    mutate(scale_R = as.factor(scale_R)) %>%
    mutate(scale_k2 = as.factor(scale_k2)) %>%
    mutate(scale_k_2 = as.factor(scale_k_2)) %>%
    mutate(scale_G0 = as.factor(scale_G0)) |> 
    mutate(term_criteria = as.factor(term_criteria))
  
  return(data_tot)
}


# -----------------------------------------------------------------------------------------------------------
# Produce the plots
# -----------------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Old_model", "pole_size_pos")
data_tot <- process_simulated_data_files(dir_experiment) |> 
  mutate(scale_R = as.numeric(as.character(scale_R))) |> 
  mutate(R = 3.95 * scale_R) |> 
  mutate(t_end = t_end * R^2 /(0.011*60)) 

data_pole <- data_tot |> 
  mutate(pole_size = ratio * 4 * pi * (R)^2) |> 
  filter(ratio < 20) |> 
  filter(R > 3.0) |> 
  filter(n_poles == 1) |> 
  filter(scale_G0 == 1) |> 
  mutate(V = R^3 * 4 * pi/3, 
         pole_diameter = 2 *sqrt(pole_size / pi))

p1 <- ggplot(data_pole, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, method="lm", color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  labs(x = "Cell radius [µm]", y = "Cdc42-GTP pole-size", title = "Model : k8h") + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]", 
       title = "Positive feedback : Pole diameter increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

p2 <- ggplot(data_pole, aes(V, u_max)) + 
  geom_smooth(linewidth=2.0, method="lm", color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  labs(x = "Cell radius [µm]", y = "Cdc42-GTP pole-size", title = "Model : k8h") + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]", 
       title = "Positive feedback : Max Cdc42 increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

mm1 <- lm(log10(pole_diameter) ~ log10(V), data = data_pole)
mm2 <- lm(log10(u_max) ~ log10(V), data = data_pole)   
summary(mm1)
summary(mm2)

dir_save <- file.path("..", "..", "..", "Results", "Cluster_area", "Old_model")
if(!dir.exists(dir_save)) dir.create(dir_save)

ggsave(file.path(dir_save, "Size_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Conc_vol.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
