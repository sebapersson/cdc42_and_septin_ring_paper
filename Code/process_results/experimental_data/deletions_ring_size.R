setwd(file.path("Code", "process_results", "experimental_data"))

library(tidyverse)
library(stringr)
library(ggthemes)
library(readxl)
library(R.matlab)

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
# Experimental results for the effect on ring size for various deletions
# -----------------------------------------------------------------------------------------------------------------

# Redo the reference fitting (how ring diameter relates to volume)
file_ref <- file.path("..", "..", "..", "Data", "SCGE_all_importwithlongaxis_20190815.mat")
data_ref <- readMat(file_ref)

data_fit <- tibble(volumes_all = data_ref$volumes.all |> as.numeric(), 
                   ring_diameters_all = data_ref$ringdiameters.all |> as.numeric()) |> 
  mutate(logvolumes_all = log((0.086666666^3)*volumes_all) / log(10) |> as.numeric(), 
         logringdiameters_all = log((0.086666666)*ring_diameters_all) / log(10))
mm <- lm(logringdiameters_all ~ logvolumes_all, data=data_fit)  
summary(mm)
ggplot(data_fit, aes(logvolumes_all, logringdiameters_all)) + 
  geom_point() + 
  geom_smooth(method = "lm")
p_haploid_all <- c(as.numeric(mm$coefficients[2]), as.numeric(mm$coefficients[1]))

# Produce the plots in the MS
path_file <- file.path("..", "..", "..", "Data", "SCGE_all_importwithlongaxis_20200512_newexps.mat")
data_matlab <- readMat(path_file)

bnr1_del <- tibble(volume = data_matlab$volumes.Bnr1del |> as.numeric(), 
                   longaxis = data_matlab$longaxis.Bnr1del |> as.numeric(), 
                   ring_diameters = data_matlab$ringdiameters.Bnr1del |> as.numeric(),
                   label = "Bnr1_del")
bni1_del <- tibble(volume = data_ref$volumes.Bni1del |> as.numeric(), 
                   longaxis = data_ref$longaxis.Bni1del |> as.numeric(), 
                   ring_diameters = data_ref$ringdiameters.Bni1del |> as.numeric(),
                   label = "Bni1_del")
spa2_del <- tibble(volume = data_matlab$volumes.Spa2del |> as.numeric(), 
                   longaxis = data_matlab$longaxis.Spa2del |> as.numeric(),
                   ring_diameters = data_matlab$ringdiameters.Spa2del |> as.numeric(),
                   label = "Spa2_del")
bud6_del <- tibble(volume = data_matlab$volumes.Bud6del |> as.numeric(), 
                   longaxis = data_matlab$longaxis.Bud6del |> as.numeric(), 
                   ring_diameters = data_matlab$ringdiameters.Bud6del |> as.numeric(),
                   label = "Bud6_del")
wt <- tibble(volume = data_matlab$volumes.haploid.WTnew |> as.numeric(), 
             longaxis = data_matlab$longaxis.haploid.WTnew |> as.numeric(), 
             ring_diameters = data_matlab$ringdiameters.haploid.WTnew |> as.numeric(),
             label = "WT_new")

data_plot <- bind_rows(bnr1_del, bni1_del,  spa2_del, bud6_del, wt) |> 
  mutate(label = factor(label, levels = c("WT_new", "Bni1_del", "Spa2_del", "Bud6_del", "Bnr1_del"))) |> 
  mutate(volume_fl = volume *  0.086666666^3,
         elongation = (pi/6)^(1/3)*0.086666666*longaxis/((0.086666666^3)*volume)^(1/3), 
         ring_diameter_plot = ring_diameters * 0.086666666, 
         logvolume = log((0.086666666^3)*volume) / log(10), 
         logringdiameters = log((0.086666666)*ring_diameters) / log(10))
# Deviation from model fit
data_plot <- data_plot |> 
  mutate(diff = logringdiameters - (p_haploid_all[2] + p_haploid_all[1]*logvolume), 
         ratio = 10^(diff))

col_use <- c(cbPalette[2], cbPalette[4], cbPalette[4], cbPalette[4], cbPalette[1])
p1 <- ggplot(data_plot, aes(label, volume_fl, color=label, fill=label)) + 
  geom_boxplot(color="black") + 
  scale_fill_manual(values = col_use) + 
  scale_color_manual(values = col_use) + 
  labs(x = "", y = "Cell volume [fL]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

p2 <- ggplot(data_plot, aes(label, elongation, color=label, fill=label)) + 
  geom_boxplot(color="black") + 
  scale_fill_manual(values = col_use) + 
  scale_color_manual(values = col_use) + 
  labs(x = "", y = "Elongation [µm]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

p3 <- ggplot(data_plot, aes(label, ring_diameter_plot, color=label, fill=label)) + 
  geom_boxplot(color="black") + 
  scale_fill_manual(values = col_use) + 
  scale_color_manual(values = col_use) + 
  labs(x = "", y = "Ring diameter [µm]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

p4 <- ggplot(data_plot, aes(label, ratio, color=label, fill=label)) + 
  geom_boxplot(color="black") + 
  scale_fill_manual(values = col_use) + 
  scale_color_manual(values = col_use) + 
  labs(x = "", y = "Ring diameter / model prediction") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

dir_save <- file.path("..", "..", "..", "Results", "Various_deletions", "Experimental_data")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

data_n <- data_plot |> 
  group_by(label) |> 
  summarise(n = n())

ggsave(file.path(dir_save, "Volume.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Elongation.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Ring_diameter.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Prediction.svg"), p4, width = BASE_WIDTH, height = BASE_HEIGHT)


# -----------------------------------------------------------------------------------------------------------------
# Hemizygote experiments
# -----------------------------------------------------------------------------------------------------------------
path_file <- file.path("..", "..", "..", "Data", "UTF-8v0_cdc24_Bem1_Bem2_Cdc42_Gic2_hemi_cdc10_plotted_data.csv")
data <- read_csv(path_file, col_types = cols()) |>
  mutate(ring_diameter = max_cdc10_ring_diameter * 0.0865202,
         cell_vol_fl = cell_vol_at_max_ring_diam_fl) |>
  filter(ring_diameter < 4) # 2 data-points

# We use WT median +- 25fL as filtering criteria for the plot
diff <- 25
wt_median <- data |>
  filter(strain_type == "WT") |>
  group_by(strain_type) |>
  summarise(median = median(cell_vol_fl)) |>
  pull(median)
data_plot <- data |>
  filter(cell_vol_fl > wt_median - diff) |>
  filter(cell_vol_fl < wt_median + diff) |>
  mutate(type_col = case_when(strain_type == "WT" ~ "WT",
                              strain_type == "BEM2/bem2Δ" ~ "BEM2",
                              T ~ "other"))

col_use <- cbPalette[c(8, 1, 2)]
p1 <- ggplot(data_plot, aes(strain_type, ring_diameter, fill = type_col)) +
  geom_boxplot() +
  labs(x = "", y = "Cdc10 ring diameter [µm]") +
  scale_x_discrete(limits = c("WT",  "CDC42/cdc42Δ", "CDC24/cdc24Δ", "BEM1/bem1Δ", "GIC2/gic2Δ", "BEM2/bem2Δ")) +
  scale_fill_manual(values = col_use) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

p2 <- ggplot(data_plot, aes(strain_type, cell_vol_fl, fill = type_col)) +
  geom_boxplot() +
  labs(x = "", y = "Cell volume [fL]") +
  scale_x_discrete(limits = c("WT",  "CDC42/cdc42Δ", "CDC24/cdc24Δ", "BEM1/bem1Δ", "GIC2/gic2Δ", "BEM2/bem2Δ")) +
  scale_fill_manual(values = col_use) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

data_plot_p3 <- data_plot |>
  filter(strain_type %in% c("WT", "BEM2/bem2Δ"))

p3 <- ggplot(data_plot_p3, aes(strain_type, ring_diameter, fill = type_col)) +
  geom_boxplot(width = 0.5) +
  labs(x = "", y = "Ring diameter [µm]") +
  scale_x_discrete(limits = c("WT", "BEM2/bem2Δ")) +
  scale_fill_manual(values = cbPalette[c(8, 2)]) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

data_n <- data_plot |> 
  group_by(strain_type) |>
  summarise(n = n())

ggsave(file.path(dir_save, "Homozygote_ring_diameter.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Homozygote_volume.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Homozygote_ring_bem2.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------------
# Ts data
# -----------------------------------------------------------------------------------------------------------------
coluse <- c("#E69F00", "#D55E00")
diff <- 25

# Cdc10
path_file <- file.path("..", "..", "..", "Data", "v0_ts_mutants_cdc10_plotted_data.csv")
data <- read_csv(path_file, col_types = cols()) |>
  mutate(ring_diameter = max_cdc10_ring_diameter * 0.0865202,
         cell_vol_fl = cell_vol_at_max_ring_diam_fl)
wt_median <- data |>
  filter(strain_type == "WT") |>
  group_by(strain_type) |>
  summarise(median = median(cell_vol_fl)) |>
  pull(median)
data_plot <- data |>
  filter(cell_vol_fl > wt_median - diff) |>
  filter(cell_vol_fl < wt_median + diff) |>
  mutate(strain_type = factor(strain_type, levels = c("WT", "TS_mutant")))

p1 <- ggplot(data_plot, aes(strain_type, ring_diameter, fill = strain_type)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = coluse) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
ggsave(file.path(dir_save, "TS_ring_diameter.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)

p2 <- ggplot(data_plot, aes(strain_type, cell_vol_fl, fill = strain_type)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = coluse) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
data_n <- data_plot |>
  group_by(strain_type) |>
  summarise(n = n())

# Exo84
path_file <- file.path("..", "..", "..", "Data", "v0_ts_mutants_exo84_plotted_data.csv")
data <- read_csv(path_file, col_types = cols()) |>
  mutate(exo84_diameter = max_exo84_cluster_diameter * 0.0865202,
         cell_vol_fl = cell_vol_at_max_exo84_diam_fl)
wt_median <- data |>
  filter(strain_type == "WT") |>
  group_by(strain_type) |>
  summarise(median = median(cell_vol_fl)) |>
  pull(median)
data_plot <- data |>
  filter(cell_vol_fl > wt_median - diff) |>
  filter(cell_vol_fl < wt_median + diff) |>
  mutate(strain_type = factor(strain_type, levels = c("WT", "TS_mutant")))

p3 <- ggplot(data_plot, aes(strain_type, exo84_diameter, fill = strain_type)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = coluse) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

p4 <- ggplot(data_plot, aes(strain_type, cell_vol_fl, fill = strain_type)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = coluse) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "grey10"),
    plot.title = element_text(color = "grey10", face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
# 0.019 p-value
wilcox.test(exo84_diameter ~ strain_type, data_plot)
data_n <- data_plot |>
  group_by(strain_type) |>
  summarise(n = n())


ggsave(file.path(dir_save, "TS_ring_diameter.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "TS_ring_diameter_vol.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "TS_exo84_diameter.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "TS_exo84_diameter_vol.svg"), p4, width = BASE_WIDTH, height = BASE_HEIGHT)