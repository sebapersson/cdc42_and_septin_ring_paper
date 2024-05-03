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
