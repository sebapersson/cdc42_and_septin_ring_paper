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
# Experimental results for the exo84 data
# -----------------------------------------------------------------------------------------------------------------
file_path <- file.path("..", "..", "..", "Data", "v0_bni1del_cdc10_plotted_data.csv")
data_f <- read_csv(file_path, col_types = cols()) |> 
  mutate(volume = cell_vol_at_max_ring_diam_fl, 
         ring_diameter = max_cdc10_ring_diameter * 0.0865202)
mybins <- cut(data_f$volume, breaks=8)
data_foo <- data_f |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(ring_diameter), 
            sd_val = sd(ring_diameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 8)
valx <- data_f |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(volume))
data_bins <- inner_join(data_foo, valx)

my_col <- cbPalette[c(4, 2)]
p1 <- ggplot(data_f, aes(volume, ring_diameter, color=strain_type)) + 
  geom_point(size=2.0, alpha=0.2) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=10) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Max Cdc10 diameter [µm]", 
       x = "Cell volume at max ring diameter [fL]", 
       title = "Title in progress") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), limits = c(0, 500)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_f |> 
  group_by(strain_type) |> 
  summarise(n = n())

path_data <- file.path("..", "..", "..", "Data", "v0_bni1del_exo84_plotted_data.csv")
data_e <- read_csv(path_data, col_types = cols()) |> 
  mutate(volume = cell_vol_at_max_exo84_diam_fl, 
         diameter = max_exo84_cluster_diameter * 0.0865202)
mybins <- cut(data_e$volume, breaks=8)
data_foo <- data_e |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(diameter), 
            sd_val = sd(diameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 10)
valx <- data_e |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(volume))
data_bins <- inner_join(data_foo, valx)
p2 <- ggplot(data_e, aes(volume, diameter, color=strain_type)) + 
  geom_point(size=2.0, alpha=0.2) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=10) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Max Exo84 cluster diameter [µm]", 
       x = "Cell volume at max Exo84 diameter [fL]", 
       title = "Title in progress") + 
  labs(y = "Mean Exo84 cluster diameter three frames before [µm]", 
       x = "Cell volume mean three frames before [fL]", 
       title = "Three frames before mean (vol and diameter)") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), limits = c(0, 500)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_e |> 
  group_by(strain_type) |> 
  summarise(n = n())

path_data <- file.path("..", "..", "..", "Data", "v0_bni1del_exo84_vs_cdc10.csv")
data_g <- read_csv(path_data, col_types = cols()) |> 
  mutate(max_cdc10_diameter = max_cdc10_diameter * 0.0865202, 
         max_exo84_diameter = max_exo84_diameter * 0.0865202)
mybins <- cut(data_g$max_cdc10_diameter, breaks=8)
data_foo <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(max_exo84_diameter), 
            sd_val = sd(max_exo84_diameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 5)
valx <- data_g |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(max_cdc10_diameter))
data_bins <- inner_join(data_foo, valx)
p3 <- ggplot(data_g, aes(max_cdc10_diameter, max_exo84_diameter, color=strain_type)) +
  geom_point(size=2.0, alpha=0.2) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=0.1) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Mean Exo84 diameter three before [µm]", 
       x = "Max Cdc10 diameter [µm]", 
       title = "Title in progress") + 
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

dir_save <- file.path("..", "..", "..", "Results", "Exo84", "Experimental_data")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
ggsave(file.path(dir_save, "Cdc11_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Exo84_vol.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Exo84_cdc10.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------------
# Actin data
# -----------------------------------------------------------------------------------------------------------------
library(R.matlab)
path_file <- file.path("..", "..", "..", "Data", "WT_and_Bni1del_data_only_Both_Exps_20211102.mat")
data_matlab <- readMat(path_file)

abp140_cluster_area <- data_matlab$WT.F.apb140AreaAll |> as.numeric()
abp140_cell_volume <- data_matlab$WT.F.apb140VolCell |> as.numeric()
bni1del_cluster_area <- data_matlab$Bni1del.F.apb140AreaAll |> as.numeric()
bni1del_cell_volume <- data_matlab$Bni1del.F.apb140VolCell |> as.numeric()

data_b <- tibble(cluster_area = abp140_cluster_area, 
                  cell_volume = abp140_cell_volume, 
                  strain_type = "WT") |> 
  bind_rows(tibble(cluster_area = bni1del_cluster_area, 
                   cell_volume = bni1del_cell_volume, 
                   strain_type = "Mutant"))
mybins <- cut(data_b$cell_volume, breaks=8)
data_foo <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(cluster_area), 
            sd_val = sd(cluster_area) / sqrt(n()), 
            n = n()) |> 
  filter(n > 10) |> 
  drop_na()
valx <- data_b |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(cell_volume))
data_bins <- inner_join(data_foo, valx)

p1 <- ggplot(data_b, aes(cell_volume, cluster_area, color=strain_type)) +
  geom_point(size=2.0, alpha=0.4) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=10) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(y = "Abp140 cluster area [µm²]", 
       x = "Cell volume [fL]", 
       title = "Title in progress") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), limits = c(0, 500)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_sum <- data_b |> 
  group_by(strain_type) |> 
  summarise(n = n())

abp140_ring_diameter <- data_matlab$WT.F.mNeptuneLengthHalfHeight |> as.numeric()
abp140_cell_volume <- data_matlab$WT.F.mNeptuneVolCell |> as.numeric()
bni1del_ring_diameter <- data_matlab$Bni1del.F.mNeptuneLengthHalfHeight |> as.numeric()
bni1del_cell_volume <- data_matlab$Bni1del.F.mNeptuneVolCell |> as.numeric()
data_c <- tibble(cluster_area = abp140_cluster_area, 
                 ring_diameter = abp140_ring_diameter, 
                 strain_type = "WT") |> 
  bind_rows(tibble(cluster_area = bni1del_cluster_area, 
                   ring_diameter = bni1del_ring_diameter, 
                   strain_type = "Mutant")) |> 
  filter(ring_diameter < 5)
mybins <- cut(data_c$cluster_area, breaks=7)
data_foo <- data_c |> 
  mutate(mybins = mybins) |> 
  group_by(strain_type, mybins) |> 
  summarise(mean_val = mean(ring_diameter), 
            sd_val = sd(ring_diameter) / sqrt(n()), 
            n = n()) |> 
  filter(n > 10)
valx <- data_c |> 
  mutate(mybins = mybins) |> 
  group_by(mybins) |> 
  summarise(mean_x = mean(cluster_area))
data_bins <- inner_join(data_foo, valx) |> 
  filter(strain_type == "WT")

p2 <- ggplot(data_c, aes(cluster_area, ring_diameter, color=strain_type)) +
  geom_point(size=2.0, alpha=0.4) + 
  geom_errorbar(data=data_bins, mapping=aes(mean_x, mean_val, ymax=mean_val + sd_val, ymin=mean_val-sd_val), 
                linewidth=1.8, width=0.1) + 
  geom_line(data=data_bins, mapping=aes(mean_x, mean_val), linewidth=1.5) + 
  labs(x = "Abp140 cluster area [µm²]", 
       y = "Cdc10 ring diameter [µm]", 
       title = "Title in progress") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
# Check percentage close to zero
data_sum <- data_c |> 
  group_by(strain_type) |> 
  summarise(close_to_zero = sum(cluster_area < 0.1) / n(),
            less_than_one = sum(cluster_area < 1.0) / n(),
            n = n())

ggsave(file.path(dir_save, "Abp140_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Cdc11_abp140.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
