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

# -----------------------------------------------------------------------------------------------------------------
# Experimental results for the cdc24_kk mutant
# -----------------------------------------------------------------------------------------------------------------
file_path <- file.path("..", "..", "..", "Data", "v0_cdc24kk_cdc11_plotted_data.csv")
data_b <- read_csv(file_path, col_types = cols()) |> 
  mutate(logvolume = log10(cell_vol_median_in_S), 
         logdiameter = log10(max_cdc11_ring_diameter * 0.0865202)) # Data given in pixels
mybins <- cut(data_b$cell_vol_at_max_ring_diam_fl, breaks=14)
data_foo <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(logdiameter), 
            sd_val = sd(logdiameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 9)
valx <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(logvolume))
data_bins <- inner_join(data_foo, valx)

my_col <- cbPalette[c(3, 2)]
p1 <- ggplot(data_b, aes(logvolume, logdiameter, color=strain_type)) + 
  geom_point(size=2.0, alpha=0.2) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=0.05) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Max Cdc11 diameter [µm]", 
       x = "Mean cell volume three frames after bud emergence [fL]", 
       title = "Title in progress") + 
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

file_path <- file.path("..", "..", "..", "Data", "v0_cdc24kk_cdc11_ring_over_time.csv")
data_c <- read_csv(file_path, col_types = cols()) |> 
  mutate(ring_diameter_mean = ring_diameter_mean * 0.0865202, 
         standard_error_mean = standard_error_mean * 0.0865202)
p2 <- ggplot(data_c, aes(frame_i, ring_diameter_mean, color = strain_type, fill = strain_type)) + 
  geom_line(linewidth=1.5) +
  geom_ribbon(aes(ymin = ring_diameter_mean - standard_error_mean, ymax = ring_diameter_mean + standard_error_mean), 
              alpha=0.5, color=NA) + 
  labs(y = "Cdc11 ring diameter [µm]", 
       x = "Frame number since bud emergence", 
       title = "Title in progress") + 
  scale_x_continuous(breaks = seq(from = -10, by = 10, to = 50)) + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

dir_save <- file.path("..", "..", "..", "Results", "Cdc24_kk", "Experimental_data")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

ggsave(file.path(dir_save, "Cdc11_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Cdc11_time.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
