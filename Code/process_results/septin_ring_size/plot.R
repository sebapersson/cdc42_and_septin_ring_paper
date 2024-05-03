library(tidyverse)
library(stringr)
library(ggthemes)
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
# Old SBE model (Axl2 v1 model)
# -----------------------------------------------------------------------------------------------------------
tag = "Different_size_exo_default_main"
dir_data1 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v1", "Different_size_exo_default_main")
dir_data2 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v1", "Different_size_exo_default_main_wide")
case_list <- c("_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k13208P333k140P0k172P5k190P833k200P5k2132P5k226P0k235P0Kt0P001Dm0P072Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding/", 
               "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k13250P0k140P0k172P5k190P833k200P5k2132P5k226P0k235P0Kt0P001Dm0P072Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding/")
process_ring_size(tag, dir_data1, dir_data2, case_list, index_list = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

# -----------------------------------------------------------------------------------------------------------
# Try NF strain (Cdc24_38A understand decoupling)
# -----------------------------------------------------------------------------------------------------------
pole_ratio_larger <- c(1.2485451274997355, 1.2062215638556766, 1.2062215638556766, 1.2062215638556766, 1.2062215638556766, 
                       1.2273833456777061, 1.2168024547666914, 1.2485451274997355, 1.1850597820336473, 1.2273833456777061) 
pole_ratio_smaller <- c(0.8887948365252354, 0.8887948365252354, 0.8570521637921913, 0.8676330547032061, 0.8993757274362502, 
                        0.8993757274362502, 0.8993757274362502, 0.8993757274362502, 0.8676330547032061, 0.8676330547032061)
hline_val <- mean(sqrt(pole_ratio_larger) / sqrt(pole_ratio_smaller))
mean_val <- mean(sqrt(pole_ratio_smaller))
data_pole <- bind_rows(tibble(tag = "Cdc24_38A", 
                              condition = "Cdc42 diameter", 
                              norm_value = sqrt(pole_ratio_larger) / mean_val), 
                       tibble(tag = "WT", 
                              condition = "Cdc42 diameter", 
                              norm_value = sqrt(pole_ratio_smaller) / mean_val)) 

# Test k20 cond1 
dir_data1 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Test_k20", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data2 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Test_k20_nf", "_k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
data1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=300) |>
  mutate(tag = "k20=0.2 WT")
data2 <- process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=300) |>
  mutate(tag = "k20=0.2 mutant")
# Condition 4 - ring formation fails to frequently
data_plot1 <- bind_rows(data1, data2) |> 
  group_by(tag, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(t = "tmax", 
         k20 = "0.2")
mean_val_ <- mean((data_plot1 |> filter(tag == "k20=0.2 WT"))$median_mean)
data_norm1 <- data_plot1 |> 
  mutate(mean_val = mean_val_)

dir_data1 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Test_k20", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P5k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data2 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Test_k20_nf", "_k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P5k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
data1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=300) |>
  mutate(tag = "k20=0.5 WT")
data2 <- process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=300) |>
  mutate(tag = "k20=0.5 mutant")
# Condition 4 - ring formation fails to frequently
data_plot2 <- bind_rows(data1, data2) |> 
  group_by(tag, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(t = "tmax", 
         k20 = "0.5")
mean_val_ <- mean((data_plot2 |> filter(tag == "k20=0.5 WT"))$median_mean)
data_norm2 <- data_plot2 |> 
  mutate(mean_val = mean_val_)

# Ring : "Mean val = 1.033 +/- 0.017"
# Cluster : "Mean val = 1.213 +/- 0.139"
mean_pole <- median((data_pole |> filter(tag == "Cdc24_38A"))$norm_value)
data_norm_tot <- bind_rows(data_norm1, data_norm2)

data_plot <- data_norm_tot |> 
  mutate(norm_value = median_mean / mean_val) |> 
  rename("condition" = "k20") |> 
  select(tag, condition, norm_value) |> 
  bind_rows(data_pole) |> 
  mutate(condition = factor(condition, levels = rev(c("Cdc42 diameter", "0.2", "0.5"))))

p1 <- ggplot(data_plot, aes(condition, norm_value)) + 
  geom_violin(aes(group = tag), position = position_dodge(width=0.0), 
              draw_quantiles = 0.5) + 
  geom_jitter(width = 0.1) + 
  geom_hline(yintercept = 1.213) + 
  geom_hline(yintercept = 1.033) + 
  coord_flip() + 
  theme_bw(base_size = 14) +
  labs(x = "", y = "Diameter normalised with wild-type simulation") + 
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))

dir_save <- file.path("..", "..", "..", "Results", "Septin_ring_v2", "Ring_size")
if( !dir.exists(dir_save) ) dir.create(dir_save)
ggsave(file.path(dir_save, "Diameter_plot_nf.svg"), p1, width = BASE_WIDTH*1.2, height = BASE_HEIGHT)

data_mean <- data_plot |> 
  group_by(tag) |> 
  summarise(mean = mean(norm_value))

# -----------------------------------------------------------------------------------------------------------
# Try Cdc24_kk strain
# -----------------------------------------------------------------------------------------------------------
dir_data1 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_speed", "Test_cdc24kk", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P2k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data2 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_speed", "Test_cdc24kk", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")

col_use <- c("#9e9ac8", "#3f007d")
data1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=350) |>
  mutate(k19 = "0.8")
data2 <- process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99", filter_t=350) |>
  mutate(k19 = "1.0")

data_plot1 <- bind_rows(data1, data2) |> 
  group_by(k19, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(t_tag = "max")

p1 <- ggplot(data_plot1, aes(k19, median_mean, fill = k19)) + 
  geom_boxplot() + 
  labs(y = "Ring diameter", 
       x = "SPR value") +
  scale_fill_manual(values = col_use, name = "Hormone conc.") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
data_mean <- data_plot1 |> 
  group_by(k19) |> 
  summarise(mean = mean(median_mean))

# Time-lapse data 
data1 <- read_csv(file.path(dir_data1, "Distribution_data.csv"), col_types = cols()) |> 
  mutate(k19 = "0.8")
data2 <- read_csv(file.path(dir_data2, "Distribution_data.csv"), col_types = cols()) |> 
  mutate(k19 = "1.0")
data_plot <- bind_rows(data1, data2) |> 
  group_by(t, k19, index) |> 
  summarise(max_Cdc42 = max(Cdc42T), 
            max_P = max(P)) |> 
  filter(max_Cdc42 > 2)
p2 <- ggplot(data_plot, aes(t, max_Cdc42, color=k19, fill = k19)) + 
  geom_line(aes(group = index), linewidth=0.1) + 
  geom_smooth() + 
  geom_hline(yintercept = 184.33) + 
  scale_color_manual(values = col_use, name = "Hormone conc.") + 
  scale_fill_manual(values = col_use, name = "Hormone conc.") + 
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold"), 
        plot.subtitle = element_text(color="grey30"))
  
dir_save <- file.path("..", "..", "..", "Results", "Septin_ring_v2", "Ring_speed")
if( !dir.exists(dir_save) ) dir.create(dir_save, recursive = T)
ggsave(file.path(dir_save, "Area_cdc24kk.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Cdc42_max_cdc24kk.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# -----------------------------------------------------------------------------------------------------------
# Normal size 
# -----------------------------------------------------------------------------------------------------------
# Function for retreiving just the size data
dir_data1 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide1standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data2 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide2standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data3 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide3standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data4 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide4standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data5 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide5standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data6 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide6standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data7 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide7standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")
dir_data8 <- file.path("..", "..", "..", "Intermediate", "Simulations", "Ring_size", "Axl2_v4", "No_cable_wide8standard", "_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P072DmFOO0P016Dms0P288Dmss0P16Dc0P16DmGdp0P072Gd10P0Gk0P017k17_alt_noCrowding")

res1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.3") |> 
  bind_rows(process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.4")) |> 
  bind_rows(process_no_cable(dir_data3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.5")) |> 
  bind_rows(process_no_cable(dir_data4, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.6")) |> 
  bind_rows(process_no_cable(dir_data5, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.7")) |> 
  bind_rows(process_no_cable(dir_data6, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.8")) |> 
  bind_rows(process_no_cable(dir_data7, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.9")) |> 
  bind_rows(process_no_cable(dir_data8, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99")) 

data_plot <- res1 |> 
  group_by(tag, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(nodes_hit = case_when(tag == "0.3" ~ 147, 
                               tag == "0.4" ~ 120, 
                               tag == "0.5" ~ 95, 
                               tag == "0.6" ~ 70, 
                               tag == "0.7" ~ 50, 
                               tag == "0.8" ~ 33, 
                               tag == "0.9" ~ 17, 
                               tag == "0.99" ~ 3), 
         total_nodes = 9451) |> 
  mutate(nodes_hit_n = nodes_hit / total_nodes * 100)

data_mean <- data_plot |> 
  group_by(tag) |> 
  summarise(mean = mean(median_inner), 
            mean_out = mean(median_outer), 
            n_nodes =mean(nodes_hit_n) )

p1 <- ggplot(data_plot, aes(nodes_hit_n, median_mean)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (mean)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggplot(data_plot, aes(nodes_hit_n, median_inner)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (inner)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p3 <- ggplot(data_plot, aes(nodes_hit_n, median_outer)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (outer)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

dir_save <- file.path("..", "..", "..", "Results", "Septin_ring_v2", "Ring_size")
if(!dir.exists(dir_save)) dir.create(dir_save)
ggsave(file.path(dir_save, "Mean_dist_set1.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Mean_inner_set1.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Mean_outer_set1.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)

data_tot <- tibble()
data_tot <- bind_rows(data_tot, mutate(data_plot, size = "Normal"))

# TODO: Revise below when simulations are done

# -----------------------------------------------------------------------------------------------------------
# Bigger size 
# -----------------------------------------------------------------------------------------------------------
dir_data1 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide1medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data2 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide2medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data3 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide3medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data4 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide4medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data5 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide5medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data6 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide6medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data7 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide7medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data8 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide8medium/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P05DmFOO0P011Dms0P2Dmss0P111Dc0P111DmGdp0P05r1P2Gd10P0Gk0P017k17_alt_noCrowding/"

res1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.3") |> 
  bind_rows(process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.4")) |> 
  bind_rows(process_no_cable(dir_data3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.5")) |> 
  bind_rows(process_no_cable(dir_data4, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.6")) |> 
  bind_rows(process_no_cable(dir_data5, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.7")) |> 
  bind_rows(process_no_cable(dir_data6, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.8")) |> 
  bind_rows(process_no_cable(dir_data7, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.9")) |> 
  bind_rows(process_no_cable(dir_data8, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99")) 

data_plot <- res1 |> 
  group_by(tag, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(nodes_hit = case_when(tag == "0.3" ~ 137, 
                               tag == "0.4" ~ 114, 
                               tag == "0.5" ~ 86, 
                               tag == "0.6" ~ 65, 
                               tag == "0.7" ~ 49, 
                               tag == "0.8" ~ 33, 
                               tag == "0.9" ~ 16, 
                               tag == "0.99" ~ 1), 
         total_nodes = 9451) |> 
  mutate(nodes_hit_n = nodes_hit / total_nodes * 100)

data_mean <- data_plot |> 
  group_by(tag) |> 
  summarise(mean = mean(median_mean), 
            mean_outer = mean(median_outer), 
            mean_inner = mean(median_inner))

p1 <- ggplot(data_plot, aes(nodes_hit_n, median_mean*1.2)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (mean)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggplot(data_plot, aes(nodes_hit_n, median_inner)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (inner)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p3 <- ggplot(data_plot, aes(nodes_hit_n, median_outer)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (outer)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

data_tot <- bind_rows(data_tot, mutate(data_plot, size = "Big"))

# -----------------------------------------------------------------------------------------------------------
# Small size 
# -----------------------------------------------------------------------------------------------------------
dir_data1 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide1small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data2 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide2small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data3 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide3small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data4 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide4small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data5 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide5small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data6 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide6small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data7 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide7small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
dir_data8 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide8small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P2k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"

#dir_data1 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide1small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data2 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide2small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data3 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide3small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data4 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide4small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data5 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide5small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data6 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide6small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data7 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide7small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"
#dir_data8 <- "../../Intermediate/Simulations/Ring_size/Axl2_v4/No_cable_wide8small/_k2b1P8k5a4P0k5b32P0k12a200P0k12b1000P0k1362P5k140P0k172P5k191P5k200P3k2132P5k2221P0k235P2k240P02k255P5Kt0P001Dm0P112DmFOO0P025Dms0P45Dmss0P25Dc0P25DmGdp0P112r0P8Gd10P0Gk0P017k17_alt_noCrowding/"

res1 <- process_no_cable(dir_data1, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.3") |> 
  bind_rows(process_no_cable(dir_data2, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.4")) |> 
  bind_rows(process_no_cable(dir_data3, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.5")) |> 
  bind_rows(process_no_cable(dir_data4, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.6")) |> 
  bind_rows(process_no_cable(dir_data5, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.7")) |> 
  bind_rows(process_no_cable(dir_data6, c(1, 2, 3, 4, 5, 6, 7, 8, 9), "0.8")) |> 
  bind_rows(process_no_cable(dir_data7, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.9")) |> 
  bind_rows(process_no_cable(dir_data8, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), "0.99")) 

data_plot <- res1 |> 
  group_by(tag, index) |> 
  summarise(median_inner = mean(d_min), 
            median_outer = mean(d_max), 
            median_mean = mean((d_min+d_max)*0.5)) |> 
  mutate(nodes_hit = case_when(tag == "0.3" ~ 179, 
                               tag == "0.4" ~ 139, 
                               tag == "0.5" ~ 110, 
                               tag == "0.6" ~ 84, 
                               tag == "0.7" ~ 62, 
                               tag == "0.8" ~ 40, 
                               tag == "0.9" ~ 20, 
                               tag == "0.99" ~ 3), 
         total_nodes = 9451) |> 
  mutate(nodes_hit_n = nodes_hit / total_nodes * 100)

data_mean <- data_plot |> 
  group_by(tag) |> 
  summarise(mean = mean(median_mean), 
            mean_outer = mean(median_outer), 
            mean_inner = mean(median_inner))

p1 <- ggplot(data_plot, aes(nodes_hit_n, median_mean*0.8)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (mean)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggplot(data_plot, aes(nodes_hit_n, median_inner)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (inner)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p3 <- ggplot(data_plot, aes(nodes_hit_n, median_outer)) + 
  geom_violin(aes(group = tag), fill="grey80", draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter (outer)", 
       title = "Ring diameter increases with wide-spread exocytosis", 
       subtitle = "At around x=1.0 the trend saturates") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

data_tot <- bind_rows(data_tot, mutate(data_plot, size = "Small"))


write_csv(data_tot, "./Ring_tot.csv")
data_tot <- read_csv("./Ring_tot.csv")

# We look closer for three of the values
data_closer <- data_tot |> 
  filter(tag %in% c("0.99", "0.9", "0.6")) |> 
  mutate(median_mean = case_when(size == "Small" ~ median_mean * 0.8, 
                                 size == "Big" ~ median_mean * 1.2, 
                                 T ~ median_mean)) |> 
  mutate(size = factor(size, levels = c("Small", "Normal", "Big")))
data_norm <- data_closer |> 
  filter(size == "Small") |> 
  group_by(tag) |> 
  summarise(mean_val = mean(median_mean))
data_plot <- inner_join(data_closer, data_norm, by = c("tag")) |> 
  mutate(tag = factor(tag, levels = c("0.99", "0.9", "0.6")))

data_plot_mean <- data_plot |> 
  group_by(tag, size) |> 
  summarise(mean = mean(median_mean / mean_val))

p1 <- ggplot(data_plot, aes(tag, median_mean / mean_val, fill = size)) + 
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter", 
       title = "Ring diameter increases with wide-spread exocytosis") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))
p2 <- ggplot(data_plot, aes(tag, median_mean, fill = size)) + 
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter", 
       title = "Ring diameter increases with wide-spread exocytosis") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

ratio_normal <- tibble(ratio = c(0.8887948365252354, 0.8887948365252354, 0.8570521637921913, 0.8676330547032061, 0.8993757274362502, 
                                 0.8993757274362502, 0.8993757274362502, 0.8993757274362502, 0.8676330547032061, 0.8676330547032061), 
                       size = "Normal", 
                       R = 2.5,
                       cluster_area = R^2 * 4 * pi * ratio, 
                       cluster_d = sqrt(cluster_area / pi))
ratio_big <- tibble(ratio = c(0.8147286001481324, 0.835890381970161, 0.8041477092371178, 0.7935668183261031, 0.8147286001481324, 
                              0.7935668183261031, 0.782985927415088, 0.8041477092371178, 0.8147286001481324, 0.8253094910591472),
                    size = "Big", 
                    R = 2.5*1.2,
                    cluster_area = R^2 * 4 * pi * ratio, 
                    cluster_d = sqrt(cluster_area / pi))
ratio_small <- tibble(ratio = c(0.9416992910803089, 0.9099566183472648, 1.0580890911014706, 1.0369273092794413, 0.9311184001692943, 
                                0.8887948365252354, 0.9946037456353825, 1.0051846365463972, 1.0580890911014706, 0.9946037456353820),
                      size = "Small", 
                      R = 2.5*0.8,
                      cluster_area = R^2 * 4 * pi * ratio, 
                      cluster_d = sqrt(cluster_area / pi))
ratio_plot <- bind_rows(ratio_small, ratio_normal, ratio_big) |> 
  mutate(index = rep(1:10, 3))
data_norm <- ratio_plot |> 
  filter(size == "Small") |> 
  group_by(size) |> 
  summarise(mean = mean(cluster_d), 
            mean_a = mean(cluster_area))
ratio_plot <- ratio_plot |> 
  mutate(size = factor(size, levels = c("Small", "Normal", "Big")))

ratio_plot_mean <- ratio_plot |> 
  group_by(size) |> 
  summarise(mean_d = mean(cluster_d / data_norm$mean))

p3 <- ggplot(ratio_plot, aes(size, cluster_d / data_norm$mean)) + 
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(width=0.005, size=3.0) + 
  labs(x = "% of 9451 nodes that can be hit by exocytosis [%]", y = "Ring diameter", 
       title = "Ring diameter increases with wide-spread exocytosis") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))

data_both <- inner_join(data_plot, ratio_plot, by = c("index", "size"))
coluse <- c("#f0f0f0", "#bdbdbd", "#636363")
p3 <- ggplot(data_both, aes(cluster_area, median_mean)) + 
  geom_smooth() +
  geom_point(aes(fill = size, shape = tag), size=3.0, color="black") + 
  scale_fill_manual(values = coluse) + 
  labs(u = "Ring diameter [µm]", x = "Cluster area [µm^2]", 
       title = "Ring diameter increases with wide-spread exocytosis") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"))


dir_save <- "../../Results/Septin_ring_v2/Sim_results/"
ggsave(str_c(dir_save, "Dist_size_norm.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(str_c(dir_save, "Dist_size.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(str_c(dir_save, "Cluster_ring.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)

