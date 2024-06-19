library(tidyverse)
library(stringr)
library(ggthemes)
library(readxl)
library(ggExtra)

# General plotting parameters
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
my_theme <- theme_bw(base_size = 16) + theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  plot.subtitle = element_text(hjust = 0.5)
) +
  theme(axis.title = element_text(size = 21))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0
SCALE_AREA <- 0.0865202^2

# -----------------------------------------------------------------------------------------------------------
# Cluster area (with hormone conc.) - Fig1
# -----------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_cdc42_max_cluster_size.csv")
data <- read_csv(path_file, col_types = cols()) |>
  mutate(strain_type = factor(strain_type, levels = c("WT", "mutant"))) |>
  mutate(ring_diameter = 2 * sqrt(max_cdc42_cluster_size / pi)) |>
  filter(log10(max_cdc42_cluster_size*SCALE_AREA) > -0.50)

# Check ring diameter vs cell volume
data_WT <- data |>
  filter(strain_type == "WT")
data_mutant <- data |> filter(strain_type == "mutant")

p1_WT <- ggplot(data_WT, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size*SCALE_AREA))) +
  geom_point(alpha = 0.7, size = 2.0, color = cbPalette[2]) +
  geom_smooth(linewidth = 2.0, method = "lm", color = "grey50", linewidth = 2.0) +
  scale_color_manual(values = cbPalette[-1]) +
  labs(
    y = "log10(Maximum Cdc42-GTP cluster size) [µm]", x = "log10(Cell volume) [fL]",
    title = "Cdc42 pole diameter increases with cell volume in WT",
    subtitle = "log10 transformed x- and y-axis"
  ) +
  scale_x_continuous(limits = c(1.5, 2.6), breaks = c(1.5, 1.8, 2.1, 2.4)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, -0.25, 0.0, 0.25, 0.5)) + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

p1_mutant <- ggplot(data_mutant, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size*SCALE_AREA))) +
  geom_point(alpha = 0.7, size = 2.0, color = cbPalette[6]) +
  geom_smooth(linewidth = 2.0, method = "lm", color = "grey50", linewidth = 2.0) +
  scale_color_manual(values = cbPalette[-1]) +
  labs(
    y = "log10(Maximum Cdc42-GTP cluster size) [µm]", x = "log10(Cell volume) [fL]",
    title = "Cdc42 pole diameter increases with cell volume in Cdc24_38A",
    subtitle = "log10 transformed x- and y-axis"
  ) +
  scale_x_continuous(limits = c(1.5, 2.6), breaks = c(1.5, 1.8, 2.1, 2.4)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, -0.25, 0.0, 0.25, 0.5)) + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

mm_WT <- lm(log10(max_cdc42_cluster_size) ~ log10(cell_vol_fl), data = data_WT)
summary(mm_WT) # slope=0.45 and p-value <2e-16
mm_mutant <- lm(log10(max_cdc42_cluster_size) ~ log10(cell_vol_fl), data = data_mutant)
summary(mm_mutant) # slope 0.39 and p-value <2e-16

data <- read_csv(path_file, col_types = cols())
data_sum <- data |>
  group_by(strain_type) |>
  summarise(n = n())

dir_save <- file.path("..", "..", "..", "Results", "Cluster_area", "Experimental_data")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
ggsave(file.path(dir_save, "Pole_diameter_vol_WT.svg"), p1_WT, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_vol_mutant.svg"), p1_mutant, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------
# Cluster area (with hormone conc.) - Fig3
# -----------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_cdc42_max_cluster_size.csv")
data <- read_csv(path_file, col_types = cols()) |> 
    filter(log10(max_cdc42_cluster_size*SCALE_AREA) > -0.5)

data_sum <- data |>
  group_by(strain_type) |>
  summarise(n = n())
my_col <- c(cbPalette[c(6, 2)])
p1 <- ggplot(data, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size*SCALE_AREA), color = strain_type)) +
  geom_point(size = 2.0, alpha = 0.2) +
  geom_smooth(method = "lm") +
  labs(
    y = "Max Cdc42-GTP cluster area [µm²]",
    x = "Cell volume [fL]",
    title = "v1_cdc24_38a_cdc42_max_cluster_size.csv"
  ) +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  ylim(-0.5, 0.5) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
p2 <- ggplot(filter(data, max_cdc42_cluster_size < 500), aes(max_cdc42_cluster_size*SCALE_AREA, cdc42_cluster_max, color = strain_type)) +
  geom_point(size = 1.0, alpha = 0.5) +
  geom_density2d() +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  labs(x = "Max Cdc42-GTP cluster area [µm²]", y = "Max Cdc42-GTP cluster intensity [A.U]") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
p2 <- ggMarginal(p2, type = "density", groupColour = TRUE, groupFill = TRUE)
sum(data$max_cdc42_cluster_size > 500) # 6 extreme cells to make the plot interpretable

ggsave(file.path(dir_save, "Pole_diameter_vol_both.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_intensity.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------
# Cluster area over time (Fig. S6)
# -----------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_whi5_cdc42_max_cluster_size.csv")
data <- read_csv(path_file, col_types = cols()) |> 
  filter(frame_i > -21 & frame_i < 21)

p1 <- ggplot(data, aes(frame_i*3, cdc42_cluster_size_um2_mean, color=strain_type, fill=strain_type)) + 
  geom_line(linewidth=1.5) + 
  geom_ribbon(aes(ymin=cdc42_cluster_size_um2_mean-cdc42_cluster_size_um2_standard_error_mean, 
                  ymax=cdc42_cluster_size_um2_mean+cdc42_cluster_size_um2_standard_error_mean), 
              alpha = 0.5) +
  labs(x = "Time since bud emergence [min]", y = "Cdc42-GTP cluster area [µm^2]") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = cbPalette[c(6, 2)], name = "") +
  scale_fill_manual(values = cbPalette[c(6, 2)], name = "") +
  theme(
    legend.position = "bottom",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
ggsave(file.path(dir_save, "Cluster_area_time.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------
# Cluster area (with hormone conc.) - FigS3
# -----------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_whi5_cdc42_max_cluster_size.csv")
data <- read_csv(path_file, col_types = cols()) |>
  filter(log10(max_cdc42_cluster_size*SCALE_AREA) > -0.5)

my_col <- c(cbPalette[c(6, 2)])
p1 <- ggplot(data, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size*SCALE_AREA), color = strain_type)) +
  geom_point(size = 2.0, alpha = 0.2) +
  geom_smooth(method = "lm") +
  labs(
    y = "Max Cdc42-GTP cluster area [µm²]",
    x = "Cell volume [fL]",
    title = "v1_cdc24_38a_cdc42_max_cluster_size.csv"
  ) +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
p2 <- ggplot(data, aes(max_cdc42_cluster_size*SCALE_AREA, cdc42_cluster_max, color = strain_type)) +
  geom_point(size = 1.0, alpha = 0.5) +
  geom_density2d() +
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  labs(x = "Max Cdc42-GTP cluster area [µm²]", y = "Max Cdc42-GTP cluster intensity [A.U]") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
p2 <- ggMarginal(p2, type = "density", groupColour = TRUE, groupFill = TRUE)

data$hormone_conc <- factor(data$hormone_conc, levels = c("0nM", "40nM", "100nM"))
coluse <- c("#efedf5", "#bcbddc", "#756bb1")
p3 <- ggplot(data, aes(strain_type, cell_vol_fl, fill = hormone_conc)) +
  geom_boxplot() +
  labs(x = "", y = "Cell volume [fL]") +
  scale_fill_manual(values = coluse) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

mm_wt <- lm(log10(max_cdc42_cluster_size*SCALE_AREA) ~ log10(cell_vol_fl),  filter(data, strain_type != "mutant"))
summary(mm_wt)

data <- read_csv(path_file, col_types = cols())
data_sum <- data |>
  group_by(strain_type, hormone_conc) |>
  summarise(n = n())

ggsave(file.path(dir_save, "Pole_diameter_vol_both_whi5.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_intensity_whi5.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_volumes_whi5.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)

data_sum <- data |>
  group_by(hormone_conc, strain_type) |>
  summarise(n = n())
