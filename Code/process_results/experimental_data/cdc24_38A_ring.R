setwd(file.path("Code", "process_results", "experimental_data"))
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
SCALE_AREA <- 0.0865202^2
SCALE_DIAMETER <- 0.0865202

# -----------------------------------------------------------------------------------------------------------------
# Diploid results
# -----------------------------------------------------------------------------------------------------------------
# Septin ring
file_path <- file.path("..", "..", "..", "Data", "v0_cdc24_38a_cdc10_plotted_data.csv")
data_b <- read_csv(file_path, col_types = cols()) |> 
  mutate(logvolume = log10(cell_vol_at_max_ring_diam_fl), 
         log_diameter = log10(max_cdc10_ring_diameter * SCALE_DIAMETER))
mybins <- cut(data_b$logvolume, breaks=7)
data_foo <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(log_diameter), 
            sd_val = sd(log_diameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 14)
valx <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(logvolume))
data_bins <- inner_join(data_foo, valx)
my_col <- c(cbPalette[6], cbPalette[2])
p1 <- ggplot(data_b, aes(logvolume, log_diameter, color=strain_type)) + 
  geom_point(size=2.0, alpha=0.3) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=0.05) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Max Cdc10 diameter [µm]", 
       x = "Cell volume at max ring diameter [fL]", 
       title = "Title in progress") + 
  scale_x_continuous(limits = c(1.47, 2.50)) +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_b |> 
  group_by(strain_type) |> 
  summarise(n = n())
data_foo <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(cdc10_ring_diameter_median_in_S), 
            sd_val = sd(cdc10_ring_diameter_median_in_S) / sqrt(n()), 
            n = n()) 
valx <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(cell_vol_median_in_S))
data_bins <- inner_join(data_foo, valx)
data_wt <- (data_bins |> filter(strain_type != "mutant"))[1:5, ]
data_mutant <- (data_bins |> filter(strain_type == "mutant"))[1:5, ]
ratio_val <- data_mutant$mean_val / data_wt$mean_val
print(sprintf("Mean val = %.3f +/- %.3f", mean(ratio_val), sd(ratio_val)))

# Cluster area
file_path <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_cdc42_strain-DLY1657_max_cluster_size.csv")
data_g <- read_csv(file_path, col_types = cols()) |>  
  mutate(logvolume = log10(cell_vol_fl), 
         logarea = log10(max_cdc42_cluster_size * SCALE_AREA))
mybins <- cut(data_g$logvolume, breaks=9)
data_foo <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(logarea), 
            sd_val = sd(logarea) / sqrt(n()), 
            n = n()) |> 
  filter(n > 13)
valx <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(logvolume))
data_bins <- inner_join(data_foo, valx)

my_col <- c(cbPalette[6], cbPalette[2])
p6 <- ggplot(data_g, aes(logvolume, logarea, color=strain_type)) + 
  geom_point(size=2.0, alpha=0.2) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=0.05) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Max Cdc42 cluster area [µm²]", 
       x = "Cell volume at max ring diameter [fL]", 
       title = "Title in progress") + 
  scale_x_continuous(limits = c(1.47, 2.50)) +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_g |> 
  group_by(strain_type) |> 
  summarise(n = n())
data_foo <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(max_cdc42_cluster_size), 
            sd_val = sd(max_cdc42_cluster_size) / sqrt(n()), 
            n = n()) 
valx <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(cell_vol_fl))
data_bins <- inner_join(data_foo, valx)
data_wt <- (data_bins |> filter(strain_type != "mutant"))[2:8, ]
data_mutant <- (data_bins |> filter(strain_type == "mutant"))[1:7, ]
ratio_val <-sqrt(data_mutant$mean_val) / sqrt(data_wt$mean_val)
print(sprintf("Mean val = %.3f +/- %.3f", mean(ratio_val), sd(ratio_val)))


dir_save <- file.path("..", "..", "..", "Results", "Cdc24_38A", "Experimental_data")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

ggsave(file.path(dir_save, "Max_cdc10_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Max_cluster_vol_HOMO.svg"), p6, width = BASE_WIDTH, height = BASE_HEIGHT)

data_WT <- data_g |> filter(strain_type == "WT")
mm_WT <- lm(logarea ~ logvolume, data_WT)
summary(mm_WT)
data_mutant <- data_g |> filter(strain_type != "WT")
mm_mutant <- lm(logarea ~ logvolume, data_mutant)
summary(mm_mutant)

ggplot(data_g, aes(logvolume, logarea)) + 
  geom_point(alpha = 0.4, size=1.0) + 
  geom_smooth(method = "lm") + 
  facet_wrap(~strain_type) + 
  theme_bw(base_size = 16)

# Over time
path_file <- file.path("..", "..", "..", "Data", "v0_cdc24_38a_cdc10_ring_over_time.csv")
data_s <- read_csv(path_file, col_types = cols()) |> 
  mutate(ring_diameter_mean = ring_diameter_mean * 0.0865202, 
         standard_error_mean = standard_error_mean * 0.0865202) # Pixel -> µm

p1 <- ggplot(data_s, aes(frame_i, ring_diameter_mean, color = strain_type, fill = strain_type)) + 
  geom_line(linewidth=1.5) +
  geom_ribbon(aes(ymin = ring_diameter_mean - standard_error_mean, ymax = ring_diameter_mean + standard_error_mean), 
              alpha=0.5, color=NA) + 
  labs(y = "Cdc10 ring diameter [µm]", 
       x = "Frame number since bud emergence", 
       title = "Title in progress") + 
  scale_x_continuous(breaks = seq(from = -10, by = 10, to = 50)) + 
  scale_color_manual(values = cbPalette[c(6, 2)]) +
  scale_fill_manual(values = cbPalette[c(6, 2)]) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_s |> 
  group_by(strain_type) |> 
  summarise(n_max = max(sample_size), 
            n_mean = mean(sample_size))

data_wt <- data_s |> 
  filter(strain_type == "WT") |> 
  select(frame_i, ring_diameter_mean) |> 
  rename("ring_diameter_wt" = "ring_diameter_mean")
data_mutant <- data_s |> 
  filter(strain_type != "WT") |> 
  filter(frame_i %in% unique(data_wt$frame_i)) |> 
  select(frame_i, ring_diameter_mean) |> 
  rename("ring_diameter_mutant" = "ring_diameter_mean")
data_plot <- inner_join(data_wt, data_mutant, by = c("frame_i"))

p2 <- ggplot(data_plot, aes(frame_i, ring_diameter_mutant / ring_diameter_wt)) + 
  geom_line(linewidth=1.0) +
  scale_y_continuous(breaks = c(0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10), limits = c(0.8, 1.10)) +
  labs(x = "Frame since bud emergence", 
       y = "Ring diameter mutant normalised with WT", 
       title = "Ring diameter normalised with WT") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "Diameter_time.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Diameter_time_norm.svg"), p6, width = BASE_WIDTH, height = BASE_HEIGHT)
