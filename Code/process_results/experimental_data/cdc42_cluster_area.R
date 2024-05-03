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

# -----------------------------------------------------------------------------------------------------------
# Cluster area (with hormone conc.) - Fig1
# -----------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "Pole_data.csv")
data <- read_csv(path_file, col_types = cols()) |> 
  mutate(hormone_conc = factor(hormone_conc, levels=c("0nM", "40nM", "100nM"))) |> 
  mutate(strain_type = factor(strain_type, levels = c("WT", "mutant"))) |> 
  mutate(ring_diameter = 2 * sqrt(max_cdc42_cluster_size / pi)) |> 
  filter(cdc42_cluster_max > 4)

# Check ring diameter vs cell volume
my_col <- c("#c7e9b4", "#41b6c4", "#225ea8")
data_WT <- data |> filter(strain_type == "WT")
data_mutant <- data |> filter(strain_type == "mutant")

p1_WT <- ggplot(data_WT, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size))) +
  geom_point(aes(color=hormone_conc), alpha=0.7, size=2.0) + 
  geom_smooth(linewidth=2.0, color="grey50", method="lm") + 
  scale_color_manual(values = cbPalette[-1]) + 
  labs(y = "log10(Maximum Cdc42-GTP cluster size) [µm]", x = "log10(Cell volume) [fL]", 
       title = "Cdc42 pole diameter increases with cell volume in WT", 
       subtitle = "log10 transformed x- and y-axis") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  scale_y_continuous(limits = c(0.3, 2.9), breaks = c(0.5, 1.5, 2.5)) + 
  scale_x_continuous(limits=c(1.3, 2.5), breaks = c(1.3, 1.6, 1.9, 2.2, 2.5)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

p1_mutant <- ggplot(data_mutant, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size))) +
  geom_point(aes(color=hormone_conc), alpha=0.7, size=2.0) + 
  geom_smooth(linewidth=2.0, color="grey50", method="lm") + 
  scale_color_manual(values = cbPalette[-1]) + 
  labs(y = "log10(Maximum Cdc42-GTP cluster size) [µm]", x = "log10(Cell volume) [fL]", 
       title = "Cdc42 pole diameter increases with cell volume in Cdc24_38A", 
       subtitle = "log10 transformed x- and y-axis") + 
  scale_y_continuous(limits = c(0.3, 2.9), breaks = c(0.5, 1.5, 2.5)) + 
  scale_x_continuous(limits=c(1.3, 2.5), breaks = c(1.3, 1.6, 1.9, 2.2, 2.5)) + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

mm_WT <- lm(log10(max_cdc42_cluster_size) ~ log10(cell_vol_fl) + hormone_conc, data=data_WT)
summary(mm_WT) # slope=0.52 and p-value <2e-16
mm_mutant <- lm(log10(max_cdc42_cluster_size) ~ log10(cell_vol_fl) + hormone_conc, data=data_mutant)
summary(mm_mutant) # slope 0.28 and p-value 1.3e-10

data_wt_sum <- data_WT |> 
  group_by(hormone_conc) |> 
  summarise(n = n())
data_mutant_sum <- data_mutant |> 
  group_by(hormone_conc) |> 
  summarise(n = n())
# Check Cdc42 max intensity
my_col <- c("#c7e9b4", "#41b6c4", "#225ea8")
p2 <- ggplot(data, aes(cell_vol_fl, cdc42_cluster_max, color=hormone_conc)) + 
  geom_point(alpha=0.8, size=2.0) + 
  scale_color_manual(values = cbPalette[-1]) + 
  facet_wrap(~strain_type) + 
  labs(y = "Maximum Cdc42-GTP intensity [A.U]", x = "Cell volume [fL]", 
       title = "Cdc42 pole intensity might increase with cell size in WT", 
       subtitle = "log10 transformed x- and y-axis") + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  geom_smooth(method="lm") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))


# Check Cdc42 max intensity and ring diameter 
my_col <- c("#c7e9b4", "#41b6c4", "#225ea8")
p3 <- ggplot(data, aes(cdc42_cluster_max, ring_diameter, color=hormone_conc)) + 
  geom_point(alpha=0.8, size=2.0) + 
  scale_color_manual(values = cbPalette[-1]) + 
  facet_wrap(~strain_type) + 
  labs(y = "Maximum Cdc42-GTP intensity [A.U]", x = "Maximum Cdc42-GTP diameter [µm]", 
       title = "Pole diameter increases with maximum Cdc42-GTP intensity", 
       subtitle = "log10 transformed x- and y-axis") + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  geom_smooth(method="lm") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

data_WT <- data |> filter(strain_type == "WT")
data_mutant <- data |> filter(strain_type == "mutant")
mm1 <- lm(log10(ring_diameter) ~ log10(cell_vol_fl) + hormone_conc, data=data_WT)
summary(mm1)
mm2 <- lm(log10(ring_diameter) ~ log10(cell_vol_fl) + hormone_conc, data=data_mutant)
summary(mm2)

dir_save <- file.path("..", "..", "..", "Results", "Cluster_area", "Experimental_data")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

ggsave(file.path(dir_save, "Pole_diameter_vol_WT.svg"), p1_WT, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_vol_mutant.svg"), p1_mutant, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Cdc42_max_vol.svg"), p2, width = BASE_WIDTH*1.5, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_cdc42_max.svg"), p3, width = BASE_WIDTH*1.5, height = BASE_HEIGHT)


# -----------------------------------------------------------------------------------------------------------
# Cluster area (with hormone conc.) - Fig3
# -----------------------------------------------------------------------------------------------------------
library(ggExtra)
path_file <- file.path("..", "..", "..", "Data", "v1_cdc24_38a_whi5_cdc42_max_cluster_size.csv")
data <- read_csv(path_file, col_types = cols()) |> 
  filter(log10(max_cdc42_cluster_size) > 1.62)

data_sum <- data |> 
  group_by(strain_type) |> 
  summarise(n = n())
my_col <- c(cbPalette[c(6, 2)])
p1 <- ggplot(data, aes(log10(cell_vol_fl), log10(max_cdc42_cluster_size), color=strain_type)) + 
  geom_point(size=2.0, alpha=0.2) + 
  geom_smooth(method="lm") + 
  labs(y = "Max Cdc42-GTP cluster area [µm²]", 
       x = "Cell volume [fL]", 
       title = "v1_cdc24_38a_cdc42_max_cluster_size.csv") + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  ylim(1, 3) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggplot(data, aes(max_cdc42_cluster_size, cdc42_cluster_max, color = strain_type)) + 
  geom_point(size=1.0, alpha=0.5) + 
  geom_density2d() + 
  scale_color_manual(values = my_col, name = "Hormone conc.") +
  scale_fill_manual(values = my_col, name = "Hormone conc.") +
  labs(x = "Max Cdc42-GTP cluster area [µm²]", y = "Max Cdc42-GTP cluster intensity [A.U]") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggMarginal(p2, type = "density",  groupColour = TRUE, groupFill = TRUE) 

ggsave(file.path(dir_save, "Pole_diameter_vol_both.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pole_diameter_intensity.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
