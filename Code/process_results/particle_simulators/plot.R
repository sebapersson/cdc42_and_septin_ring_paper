library(tidyverse)
library(stringr)
library(ggthemes)
library(readxl)

# General plotting parameters 
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_theme <- theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=21))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0

# -----------------------------------------------------------------------------------------------------------------
# Function for processing particle simulator without and with diffusion 
# -----------------------------------------------------------------------------------------------------------------
# Process particles without diffusion
process_particle_results <- function(tag_check)
{
  
  dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Particle_no_diff")
  files_in_dir <- list.files(dir_results)
  
  files_result <- files_in_dir[which(!is.na(str_match(files_in_dir, str_c("^", tag_check))))]
  
  data_tot <- tibble()
  for(i in 1:length(files_result)){
    pattern_cover <- "cover(\\d|\\.)+"
    cover_value <- as.numeric(str_replace(str_extract(files_result[i], pattern_cover), "cover", ""))
    pattern_alpha <- "alpha(\\d|\\.)+"
    alpha_value <- as.numeric(str_sub(str_replace(str_extract(files_result[i], pattern_alpha), "alpha", ""), end=-2))
    
    
    data_tmp <- read_csv(file.path(dir_results, files_result[i]), col_types = cols()) |> 
      mutate(tag = tag_check, 
             alpha = alpha_value, 
             cover = cover_value)
    data_tot <- bind_rows(data_tot, data_tmp)
  }
  
  dir_save <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", tag_check)
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
  
  data_plot <- data_tot |> 
    mutate(alpha = factor(alpha), 
           cover = factor(cover))
  data_test <- data_plot |> filter(alpha == 0.4) |> drop_na()
  
  for(alpha_value in unique(data_plot$alpha)){
    data_plot__ <- data_plot |> filter(alpha == alpha_value) |> drop_na()
    
    file_save1 <- file.path(dir_save, str_c("Inner_alpha", as.character(alpha_value), ".svg"))
    file_save2 <- file.path(dir_save, str_c("Outer_alpha", as.character(alpha_value), ".svg"))
    
    data_inner <- data_plot__ |> filter(point_type == "inner")
    data_outer <- data_plot__ |> filter(point_type == "outer")
    
    if(nrow(data_inner) > 0){
      p1 <- ggplot(data_inner, aes(strength, distance, color = point_type)) + 
        geom_point(size=0.5, alpha=0.4) + 
        geom_smooth() + 
        scale_color_manual(values = cbPalette[c(3)], name = "Diameter measured") + 
        facet_wrap(~cover, scale="free") + 
        scale_y_log10() + 
        my_theme +
        labs(y = "Ring diameter", x = "Degree exocytes are concentrated in the centre", 
             title = "Ring diameter when 'pole' covers 0.01 - 0.05 membrane area") +
        theme(legend.position = "bottom")
      ggsave(file_save1, p1, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.0, dpi=300)
    }
    p2 <- ggplot(data_outer, aes(strength, distance, color = point_type)) + 
      geom_point(size=0.5, alpha=0.4) + 
      geom_smooth() + 
      scale_color_manual(values = cbPalette[c(4)], name = "Diameter measured") + 
      facet_wrap(~cover, scale="free") + 
      my_theme +
      labs(y = "Ring diameter", x = "Degree exocytes are concentrated in the centre", 
           title = "Ring diameter when 'pole' covers 0.01 - 0.05 membrane area") +
      theme(legend.position = "bottom")
    ggsave(file_save2, p2, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.0, dpi=300)
  }
  
  return(data_tot)
}

# Process particles with diffusion 
process_diffusion <- function(tag_check, with_exo=T, low_diffusion=F, check_recruited=F, check_cables=F)
{
  if(check_recruited == T && check_cables == F){
    dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Particle_with_diff", "With_diffusion")
  }else if(check_recruited == F && check_cables == T){
    dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Particle_with_diff", "With_diffusion_cables")
  }else{
    dir_results <- file.path("..", "..", "..", "Intermediate", "Simulations", "Particle_with_diff", "With_diffusion_reac_rec")
  }
  
  files_in_dir <- list.files(dir_results)
  if(tag_check == "amount"){ # Spelling is hard
    files_result <- files_in_dir[which(!is.na(str_match(files_in_dir, str_c("^ammount"))))]
  }else{
    files_result <- files_in_dir[which(!is.na(str_match(files_in_dir, str_c("^", tag_check))))]
  }
  
  if(with_exo == T){
    iuse <- which(is.na(str_match(files_result, str_c("noexo"))))
    files_result <- files_result[iuse]
  }else{
    iuse <- which(!is.na(str_match(files_result, str_c("noexo"))))
    files_result <- files_result[iuse]
  }
  if(low_diffusion == F){
    iuse <- which(is.na(str_match(files_result, str_c("low_diffusion"))))
    files_result <- files_result[iuse]    
  }else{
    iuse <- which(!is.na(str_match(files_result, str_c("low_diffusion"))))
    files_result <- files_result[iuse]    
  }
  
  
  data_tot <- tibble()
  for(i in 1:length(files_result)){
    pattern_cover <- "frac_cover(\\d|\\.)+"
    cover_value <- as.numeric(str_sub(str_replace(str_extract(files_result[i], pattern_cover), "frac_cover", ""), end=-2))
    pattern_Dm <- "Dm(\\d|\\.)+"
    Dm_value <- as.numeric(str_replace(str_extract(files_result[i], pattern_Dm), "Dm", ""))
    Dm_value <- ifelse(is.na(Dm_value), 0.0045, Dm_value)
    pattern_koff <- "koff(\\d|\\.)+"
    koff_value <- as.numeric(str_replace(str_extract(files_result[i], pattern_koff), "koff", ""))
    koff_value <- ifelse(is.na(koff_value), 0.1, koff_value)
    pattern_amount <- "ammount(\\d|\\.)+"
    amount_value <- as.numeric(str_replace(str_extract(files_result[i], pattern_amount), "ammount", ""))
    amount_value <- ifelse(is.na(amount_value), 100, amount_value)
    pattern_alpha <- "alpha(\\d|\\.)+"
    alpha_value <- as.numeric(str_replace(str_extract(files_result[i], pattern_alpha), "alpha", ""))
    alpha_value <- ifelse(is.na(alpha_value), 0.5, alpha_value)
    
    data_tmp <- read_csv(file.path(dir_results, files_result[i]), col_types = cols()) |> 
      mutate(tag = tag_check, 
             Dm = Dm_value, 
             koff = koff_value,
             amount = amount_value,
             alpha=alpha_value,
             cover = cover_value)
    data_tot <- bind_rows(data_tot, data_tmp)
  }
  
  colnames(data_tot)[colnames(data_tot) == tag_check] <- "what_plot"
  
  data_line <- data_tot |> 
    group_by(what_plot, weight_centre, cover, runindex) |> 
    summarise(median = median(dist), 
              quant90 = quantile(dist, 0.9), 
              quant10 = quantile(dist, 0.1))
  
  
  if(check_recruited == F && check_cables == F){
    dir_save <- file.path("..", "..", "..", "Results", "Particle_simulator_diff", "Diffusion", 
                          str_c(tag_check, ifelse(with_exo, "", "noexo"), ifelse(!low_diffusion, "", "lowdiff")))
  }else if(check_cables == T){
    dir_save <- file.path("..", "..", "..", "Results", "Particle_simulator_diff", "Diffusion_cables", 
                          str_c(tag_check, ifelse(with_exo, "", "noexo"), ifelse(!low_diffusion, "", "lowdiff")))
  }else{
    dir_save <- file.path("..", "..", "..", "Results", "Particle_simulator_diff", "Diffusion_rec", 
                          str_c(tag_check, ifelse(with_exo, "", "noexo"), ifelse(!low_diffusion, "", "lowdiff")))
  }
  
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  
  for(i in 1:length(unique(data_tot$cover))){
    what_cover <- unique(data_tot$cover)[i]
    
    data_violin <- data_tot |> 
      filter(cover == what_cover) |> 
      mutate(weight_centre_fct = as.factor(weight_centre)) |> 
      mutate(what_plot = as.factor(what_plot))
    data_line_ <- data_line |> 
      filter(cover == what_cover) |> 
      mutate(what_plot = as.factor(what_plot))
    
    p1 <- ggplot(data_violin, aes(weight_centre_fct, dist, fill = what_plot)) + 
      geom_violin(draw_quantiles = c(0.5, 0.9)) + 
      facet_wrap(~what_plot, scale="free") +
      scale_fill_manual(values = cbPalette[-1], name = tag_check) +
      my_theme
    
    data_line_ <- data_line |> filter(cover == what_cover) |> mutate(what_plot = as.factor(what_plot))
    data_line__ <- tibble() |> 
      bind_rows(data_line_ |> select(what_plot, weight_centre, cover, runindex, median) |> 
                  rename("dist" = "median") |> 
                  mutate(crit_plot = "median")) |> 
      bind_rows(data_line_ |> select(what_plot, weight_centre, cover, runindex, quant10) |> 
                  rename("dist" = "quant10") |> 
                  mutate(crit_plot = "quant10")) |> 
      bind_rows(data_line_ |> select(what_plot, weight_centre, cover, runindex, quant90) |> 
                  rename("dist" = "quant90") |> 
                  mutate(crit_plot = "quant90")) |> 
      mutate(what_plot = as.factor(what_plot))
    data_mean <- data_line_ |> 
      group_by(what_plot) |> 
      summarise(mean_median = median(median), 
                mean_quant90 = median(quant90), 
                mean_quant10 = median(quant10))
    data_mean_ <- tibble() |> 
      bind_rows(data_mean |> 
                  select(what_plot, mean_median) |> 
                  rename("mean_val" = "mean_median") |> 
                  mutate(crit_plot = "median")) |> 
      bind_rows(data_mean |> 
                  select(what_plot, mean_quant10) |> 
                  rename("mean_val" = "mean_quant10") |> 
                  mutate(crit_plot = "quant10")) |> 
      bind_rows(data_mean |> 
                  select(what_plot, mean_quant90) |> 
                  rename("mean_val" = "mean_quant90") |> 
                  mutate(crit_plot = "quant90"))
    data_plot_line <- inner_join(data_line__, data_mean_, by=c("what_plot", "crit_plot")) |> 
      filter(crit_plot != "quant10") 
    
    p2 <- ggplot(data_plot_line, aes(weight_centre, dist / mean_val, color=crit_plot)) + 
      geom_point(size=3.0) + 
      geom_smooth(se=F) +
      facet_wrap(~what_plot) + 
      geom_hline(yintercept = 1.0) + 
      scale_color_manual(values = cbPalette[-1]) +
      scale_y_continuous(breaks=seq(from=0.7, by=0.1, to=1.5)) + 
      my_theme
    
    path_save1 <- file.path(dir_save, str_c("Violin_", tag_check, "_cover", as.character(what_cover), ".svg"))
    path_save2 <- file.path(dir_save, str_c("Line_", tag_check, "_cover", as.character(what_cover), ".svg"))
    ggsave(path_save1, p1, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*2.0)
    ggsave(path_save2, p2, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*2.0)
    
  }
  return(data_tot)
}

# -----------------------------------------------------------------------------------------------------------------
# Illustration of the beta-parameter
# -----------------------------------------------------------------------------------------------------------------
path <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Dist_data_weight0.1.csv")
data1 <- read_csv(path, col_types = cols()) |> 
  mutate(weight = "0.1")
path <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Dist_data_weight0.5.csv")
data2 <- read_csv(path, col_types = cols()) |> 
  mutate(weight = "0.5")
path <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Dist_data_weight1.0.csv")
data3 <- read_csv(path, col_types = cols()) |> 
  mutate(weight = "1.0")
path <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Dist_data_weight1.5.csv")
data4 <- read_csv(path, col_types = cols()) |> 
  mutate(weight = "1.5")
path <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Dist_data_weight2.0.csv")
data5 <- read_csv(path, col_types = cols()) |> 
  mutate(weight = "2.0")
data_plot <- bind_rows(data1, data2, data3, data4, data5)
data_mean <- data_plot |> 
  group_by(weight) |> 
  summarise(mean = mean(dist_centre))
p <- ggplot(data_plot, aes(weight, dist_centre, fill=weight)) + 
  geom_violin(bw=0.03, draw_quantiles=0.5, linewidth=1.0) + 
  scale_y_continuous(limits=c(0, 0.53), 
                     breaks = c(0, 0.25, 0.5)) + 
  scale_fill_brewer(palette = "Oranges", name = "β") + 
  labs(y = "Distance to pole centre", x = "β", 
       title = "With increased β the exocytosis distribution becomes more centred", 
       subtitle = "horisontal lines correspond to the median for each β") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))
path_save <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Hit_dist", "Illustration.svg")
ggsave(path_save, p, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------------
# Processing the simple particle model (without diffusion)
# -----------------------------------------------------------------------------------------------------------------
tag_check <- "Frac_cover"
data_tot <- process_particle_results("Frac_cover")
data_tot_ <- data_tot |> 
  group_by(id, point_type, alpha, cover, strength) |> 
  summarise(dist = mean(distance, drop.na=T)) |> 
  filter(point_type == "outer")
data_plot1 <- data_tot_ |> filter(cover == 0.02 & alpha == 0.4) |> filter(strength < 3) |> drop_na()
data_plot2 <- data_tot_ |> filter(cover == 0.02 & alpha == 0.8) |> filter(strength < 3) |> drop_na()
data_plot <- bind_rows(data_plot1, data_plot2) |> 
  rename("α" = "alpha")
p <- ggplot(data_plot, aes(strength, dist)) + 
  geom_smooth(linewidth=2.0, method="lm", color="grey30") + 
  geom_point(size=3.0, color="grey20", alpha=0.8) + 
  facet_wrap(~α, scale = "free_y") + 
  theme_bw(base_size = 14) +
  labs(x = "β", y = "Ring diameter", 
       title = "Particle simulator : Ring diameter does not increase with wide-spread exocytosis (smaller β)", 
       subtitle = "Results robust to tuning parameters (e.g. α=0.4, 0.8). Each dot = 1 simulation") + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
path_save <- file.path("..", "..", "..", "Results", "Particle_simulator_no_diff", "Diameter_const.svg")
ggsave(path_save, p, width = BASE_WIDTH*1.7, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------------
# Process results for the more complex particle simulator (with diffusion)
# -----------------------------------------------------------------------------------------------------------------
data_tot <- process_diffusion("Dm")
data_plot <- data_tot |> 
  filter(what_plot %in% c(0.00045, 0.0045)) |> 
  mutate(what_plot = factor(what_plot)) |> 
  filter(cover == 0.03) |> 
  group_by(what_plot, runindex, weight_centre) |> 
  summarise(dist_mean = mean(dist), 
            dist80 = quantile(dist, 0.8))
p <- ggplot(data_plot, aes(weight_centre, dist_mean)) + 
  geom_smooth(linewidth=2.0, color="grey30", method="lm") + 
  geom_point(size=3.0, color="grey20", alpha=0.8) + 
  facet_wrap(~what_plot) + 
  theme_bw(base_size = 14) +
  ylim(0.7, 2.4) +
  labs(x = "β", y = "Ring diameter", 
       title = "Particle simulator : Ring diameter increase with less wide-spread exocytosis (larger β)", 
       subtitle = "Effect bigger with lower diffusion dm=0.00045. Each dot = 1 simulation and koff=0.1") + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
mm <- lm(dist_mean ~ weight_centre, data = filter(data_plot, what_plot == 0.00045))
summary(mm)
mm <- lm(dist_mean ~ weight_centre, data = filter(data_plot, what_plot == 0.0045))
summary(mm)
path_save <- file.path("..", "..", "..", "Results", "Particle_simulator_diff", "Diameter_dm.svg")
ggsave(path_save, p, width = BASE_WIDTH*1.7, height = BASE_HEIGHT)

data_tot <- process_diffusion("koff")
data_plot <- data_tot |> 
  filter(what_plot %in% c(0.01, 0.05)) |> 
  mutate(what_plot = factor(what_plot)) |> 
  filter(cover == 0.03) |> 
  group_by(what_plot, runindex, weight_centre) |> 
  summarise(dist_mean = mean(dist), 
            dist80 = quantile(dist, 0.8))
p <- ggplot(data_plot, aes(weight_centre, dist_mean)) + 
  geom_smooth(linewidth=2.0, color="grey30", method="lm") + 
  geom_point(size=3.0, color="grey20", alpha=0.8) + 
  facet_wrap(~what_plot) + 
  ylim(0.7, 2.4) +
  theme_bw(base_size = 14) +
  labs(x = "β", y = "Ring diameter", 
       title = "Particle simulator : Ring diameter increase with less wide-spread exocytosis (larger β)", 
       subtitle = "Effect bigger with stronger koff = 0.05. Dm=0.00045") + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
mm <- lm(dist_mean ~ weight_centre, data = filter(data_plot, what_plot == 0.01))
summary(mm)
mm <- lm(dist_mean ~ weight_centre, data = filter(data_plot, what_plot == 0.05))
summary(mm)
path_save <- file.path("..", "..", "..", "Results", "Particle_simulator_diff", "Diameter_koff.svg")
ggsave(path_save, p, width = BASE_WIDTH*1.7, height = BASE_HEIGHT)
