# Help functions for plotting results
source("cluster_area_functions.R")

# This file produces the plots related to Cdc42 cluster area

dir_save <- file.path("..", "..", "..", "Results", "Cluster_area", "Model")
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)

# -----------------------------------------------------------------------------------------------
# Constant conc simulations
# -----------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_pos_const_conc")
data_tot <- process_files(dir_experiment) |> 
  filter(cdc42_max > 100) |> 
  mutate(V = 4*pi/3 * r_num^3, 
         pole_diameter = 2*sqrt(pole_size / pi))

col_use <- "#dadaeb"
p_const <- ggplot(data_tot, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, color=col_use, fill=col_use) + 
  geom_point(size=3.0, color=col_use) + 
  labs(x = "Cell volume [fL]", y = "Cdc42-GTP pole-diameter [µm]", 
       title = "β = 0: With constant amount pole size does not increase with V") + 
  scale_x_continuous(limits = c(1.8, 2.45), breaks = c(1.8, 2.0, 2.2, 2.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "Pos_size_vol_const.svg"), p_const, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Positive feedback model 
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_pos_weak_inc_beta100")
data_tot <- process_files(dir_experiment)
data_tot <- data_tot |> 
  mutate(V = r_num^3 * 4*pi / 3) |> 
  mutate(pole_diameter = 2*sqrt(pole_size / pi))

col_use <- c("#4eb3d3", "#08589e")
p1 <- ggplot(data_tot, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, method="lm", color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  labs(x = "Cell radius [µm]", y = "Cdc42-GTP pole-area [µm²]") +
  scale_color_manual(values = col_use, name = "Factor Cdc42 reduced") + 
  scale_fill_manual(values = col_use, name = "Factor Cdc42 reduced") +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, -0.25, 0.0, 0.25, 0.5)) + 
  labs(y = "log10(Cdc42-GTP pole area) [µm]", x = "log10(Cell volume) [fL]", 
       title = "Positive feedback model: Pole area increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

mm <- lm(log10(pole_size) ~ log10(V), data=data_tot)
summary(mm) # Slope ~ 0.33 (p-value 2e-16)

data_dist <- process_files(dir_experiment, file_name = "/Pole_dist_end.csv") 
data_dist_ <- data_dist |> 
  group_by(scale_r, scale_k8max) |> 
  summarise(max_val = max(BemGEF42 + Cdc42T ), 
            r_num = median(r_num)) |> 
  mutate(V = r_num^3 * 4*pi/3)
data_plot_both <- inner_join(data_dist_, data_tot) 

p1_alt <- ggplot(data_plot_both, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  labs(x = "Cell radius [µm]", y = "Cdc42-GTP pole-area [µm²]") +
  scale_color_manual(values = col_use, name = "Factor Cdc42 reduced") + 
  scale_fill_manual(values = col_use, name = "Factor Cdc42 reduced") +
  scale_y_continuous(limits = c(-0.1, 0.86)) +
  labs(y = "log10(Cdc42-GTP pole area) [µm]", x = "log10(Cell volume) [fL]", 
       title = "Positive feedback model: Pole area increases with cell volume") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

p2 <- ggplot(data_plot_both, aes(pole_size, max_val)) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  labs(y = "log10(Max Cdc42-GTP concentration) [A.U]", x = "log10(Cdc42-GTP pole area) [µm²]", title = "Model pos") + 
  scale_color_manual(values = col_use, name = "Factor Cdc42 reduced") + 
  scale_fill_manual(values = col_use, name = "Factor Cdc42 reduced") +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "Pos_size_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pos_size_vol_alt.svg"), p1_alt, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pos_size_nax.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Negative feedback model
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_nf_amount_weak_beta100")
data_part <- process_files(dir_experiment) |> 
  mutate(V = r_num^3 * 4*pi / 3) |> 
  mutate(pole_diameter = 2*sqrt(pole_size / pi))

p1 <- ggplot(data_part, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, method="lm", color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  scale_fill_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_color_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_y_continuous(limits = c(-0.1, 0.86)) +
  labs(x = "log10(Max Cdc42-GTP concentration) [A.U]", x = "log10(Cdc42-GTP pole area) [µm²]", title = "Model : k8h") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

mm <- lm(log10(pole_size) ~ log10(V), data=data_part)
summary(mm) # 0.35 (p-value <2e-16)

data_dist <- process_files(dir_experiment, file_name = "/Pole_dist_end.csv") 
data_dist_ <- data_dist |> 
  group_by(scale_r, scale_k8max) |> 
  summarise(max_val = max(BemGEF42 + Cdc42T + BemGEF42_star), 
            r_num = median(r_num), 
            star_max = max(BemGEF42_star)) |> 
  mutate(V = r_num^3 * 4*pi/3)
data_plot_both <- inner_join(data_dist_, data_part) 

p2 <- ggplot(data_plot_both, aes(pole_size, max_val)) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  scale_fill_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_color_manual(values = cbPalette[-1], name = "Strength k8h") + 
  labs(y = "log10(Max Cdc42-GTP concentration) [A.U]", x = "log10(Cdc42-GTP pole area) [µm²]", title = "Model : k8h") + 
  scale_y_continuous(breaks = c(80, 90, 100, 110)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

p3 <- ggplot(data_plot_both, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  scale_fill_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_color_manual(values = cbPalette[-1], name = "Strength k8h") + 
  labs(x = "log10(Cell volume) [fL]", y = "log10(Cdc42-GTP pole area) [µm²]", title = "Model : k8h") + 
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, -0.25, 0.0, 0.25, 0.5)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "NF_normal_size_max.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "NF_normal_vol_size.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Neg_size_vol.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Experimenting with feedback and GAP activity for negative feedback model
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_nf_test_GAP")
data_part <- process_files(dir_experiment) |> 
  mutate(V = r_num^3 * 4*pi / 3) |> 
  mutate(pole_diameter = 2*sqrt(pole_size / pi)) 

col_use <- c("#9ebcda", "#8c6bb1")
p2 <- ggplot(data_part, aes(V, pole_diameter, color=scale_k2b, fill=scale_k2b)) + 
  geom_smooth(linewidth=2.0, method="lm") + 
  geom_point(size=3.0) + 
  scale_fill_manual(values = col_use, name = "GAP activity") + 
  scale_color_manual(values = col_use, name = "GAP activity") + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]", 
       title = "Negative feedback : Pole diameter increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]", 
       title = "Negative feedback : Pole diameter increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "nothing", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

data_dist <- process_files(dir_experiment, file_name = "/Pole_dist_end.csv") 
data_dist_ <- data_dist |> 
  group_by(scale_r, scale_k2b) |> 
  summarise(max_val = max(BemGEF42 + Cdc42T + BemGEF42_star), 
            r_num = median(r_num), 
            star_max = max(BemGEF42_star)) |> 
  mutate(V = r_num^3 * 4*pi/3)
data_plot_both <- inner_join(data_dist_, data_part) 

p1 <- ggplot(data_plot_both, aes(V, max_val, color=scale_k2b, fill=scale_k2b)) + 
  geom_smooth(linewidth=2.0) + 
  geom_point(size=3.0) + 
  scale_fill_manual(values = col_use, name = "GAP activity") + 
  scale_color_manual(values = col_use, name = "GAP activity") + 
  labs(y = "Max Cdc42-GTP concentration [A.U]", x = "Cell volume [fl]", title = "GAP activity") + 
  scale_y_log10() + 
  scale_x_log10() + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

mm1 <- lm(log10(pole_diameter) ~ log10(V), data=filter(data_part, scale_k2b==0.457))
summary(mm1) # 0.27
mm2 <- lm(log10(pole_diameter) ~ log10(V), data=filter(data_part, scale_k2b==1.0))
summary(mm2) # 0.27
ggsave(file.path(dir_save, "Neg_size_vol_GAP.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Neg_conc_vol_GAP.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Experimenting with GAP activity positive feedback model
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_pos_test_GAP")
data_part <- process_files(dir_experiment) |> 
  mutate(V = r_num^3 * 4*pi / 3) |> 
  mutate(pole_diameter = 2*sqrt(pole_size / pi)) |> 
  filter(scale_k2b != 0.571)

col_use <- c("#6baed6", "#2171b5")
p1 <- ggplot(data_part, aes(log10(V), log10(pole_size), color=scale_k2b, fill=scale_k2b)) + 
  geom_smooth(linewidth=2.0) + 
  geom_point(size=3.0) + 
  scale_fill_manual(values = col_use, name = "GAP activity") + 
  scale_color_manual(values = col_use, name = "GAP activity") + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]") + 
  labs(y = "Maximum Cdc42-GTP cluster area [µm²]", x = "Cell volume [fL]", 
       title = "Negative feedback : Pole diameter increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

datamax <- get_cdc42max(dir_experiment) |> 
  filter(scale_k2b != 0.571)
p2 <- ggplot(datamax, aes(log10(r_num^3*4*pi/3), log10(cdc42max), color=scale_k2b, fill=scale_k2b)) + 
  geom_smooth() + 
  geom_point(size=3.0) + 
  scale_fill_manual(values = col_use, name = "GAP activity") + 
  scale_color_manual(values = col_use, name = "GAP activity") + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]", 
       title = "Negative feedback : Pole diameter increases with cell volume", 
       subtitle = "log10 transformed x- and y-axis") + 
  labs(y = "Maximum Cdc42-GTP diameter [µm]", x = "Cell volume [fL]") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "Pos_size_vol_GAP.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Pos_conc_vol_GAP.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Experimenting with smaller scaling pos (protein dilution positive feedback model)
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_pos_weak_inc_beta50")
data1 <- process_files(dir_experiment) |> mutate(beta = "0.50", V = r_num^3 * 4*pi/3)
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_pos_weak_inc_beta75")
data2 <- process_files(dir_experiment) |> mutate(beta = "0.75", V = r_num^3 * 4*pi/3)
col_use <- c("#99d8c9", "#238b45", "#238b48")

data1_pos <- data1
data2_pos <- data2

# ----------------------------------------------------------------------------------------------------
# Experimenting with smaller scaling (dilution) nf model
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_nf_amount_weak_beta50")
data1 <- process_files(dir_experiment) |> mutate(beta = "0.50")
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_nf_amount_weak_beta75")
data2 <- process_files(dir_experiment) |> mutate(beta = "0.75")
data_tot <- bind_rows(data1, data2) |> 
  filter(pole_ratio < 90) |> 
  mutate(V = r_num^3 * 4*pi/3, 
         pole_diameter = 2 * sqrt(pole_size / pi))
data50 <- data_tot |> filter(beta == "0.50")
data75 <- data_tot |> filter(beta == "0.75")

data1_nf <- data50
data2_nf <- data75

col_use <- "#9e9ac8"
data1_plot <- bind_rows(mutate(data1_pos, model = "Pos"), mutate(data1_nf, model = "Neg"))
p1 <- ggplot(data1_plot, aes(log10(V), log10(pole_size), color=model)) + 
  geom_point(aes(shape=model), size=3.0, color=col_use) + 
  geom_smooth(aes(linetype=model), linewidth=2.0, method="lm", color=col_use, fill=col_use) + 
  labs(y = "log10(Cdc42-GTP pole area) [µm²]", x = "log10(Cell volume) [fL]", 
       title = "β = 0.5") + 
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  scale_x_continuous(limits = c(1.8, 2.45), breaks = c(1.8, 2.0, 2.2, 2.4)) +
  theme_bw(base_size = 14) + 
  theme(axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"), 
        legend.position = "nothing")

col_use <- "#6a51a3"
data2_plot <- bind_rows(mutate(data2_pos, model = "Pos"), mutate(data2_nf, model = "Neg"))
p2 <- ggplot(data2_plot, aes(log10(V), log10(pole_size), color=model)) + 
  geom_point(aes(shape=model), size=3.0, color=col_use) + 
  geom_smooth(aes(linetype=model), linewidth=2.0, method="lm", color=col_use, fill=col_use) + 
  labs(y = "log10(Cdc42-GTP pole area) [µm²]", x = "log10(Cell volume) [fL]", 
       title = "β = 0.75") + 
  scale_x_continuous(limits = c(1.8, 2.45), breaks = c(1.8, 2.0, 2.2, 2.4)) +
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  theme_bw(base_size = 14) + 
  theme(axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"), 
        legend.position = "nothing")

ggsave(file.path(dir_save, "Both_beta50.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Both_beta75.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)


# Visual to explore the scaling 
r_ref <- 2.5
n_ref <- 1.0 * 4*pi/3*r_ref^3
r_list <- seq(from=2.5, to=4, length.out=100)
data_scaling1 <- tibble(r = r_list, 
                        beta=1.0,
                        V = pi*4/3*r^3, 
                        ratio_ = ((r-r_ref)/r_ref*beta+1)^3,
                        c = n_ref * ratio_ / V)
data_scaling2 <- tibble(r = r_list, 
                        beta=0.75,
                        V = pi*4/3*r^3, 
                        ratio_ = ((r-r_ref)/r_ref*beta+1)^3,
                        c = n_ref * ratio_ / V)
data_scaling3 <- tibble(r = r_list, 
                        beta=0.5,
                        V = pi*4/3*r^3, 
                        ratio_ = ((r-r_ref)/r_ref*beta+1)^3,
                        c = n_ref * ratio_ / V)
data_scaling4 <- tibble(r = r_list, 
                        beta=0.0,
                        V = pi*4/3*r^3, 
                        ratio_ = ((r-r_ref)/r_ref*beta+1)^3,
                        c = n_ref * ratio_ / V)
data_scaling <- bind_rows(data_scaling1, data_scaling2, data_scaling3, data_scaling4) |> 
  mutate(beta = as.factor(beta))

col_use <- rev(c("#3f007d", "#6a51a3", "#9e9ac8", "#dadaeb"))
p2 <- ggplot(data_scaling, aes(log10(V), log10(c), color = beta)) + 
  geom_line(linewidth=2.0) + 
  scale_x_continuous(limits = c(1.8, 2.45), breaks = c(1.8, 2.0, 2.2, 2.4)) +
  labs(x = "Cell volume [fL]", y = "Cdc42-GTP totalt concentration [µM]", 
       title = "β = 0 amount constant, and β = 1 concentration constant") + 
  scale_color_manual(values = col_use) + 
  theme_bw(base_size = 14) + 
  theme(axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"), 
        legend.position = "nothing")

ggsave(file.path(dir_save, "NF_weaker_50.svg"), p1_50, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "NF_weaker_75.svg"), p1_75, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "Illustration_beta.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# ----------------------------------------------------------------------------------------------------
# Stronger feedback flattening intensity for negative feedback model
# ----------------------------------------------------------------------------------------------------
dir_experiment <- file.path("..", "..", "..", "Intermediate", "Simulations", "Pole_size", "2d_model_nf_beta100_strong_feedback")
data_part <- process_files(dir_experiment) |> 
  mutate(V = r_num^3 * 4*pi/3, 
         pole_diameter = 2 * sqrt(pole_size / pi))

data_dist <- process_files(dir_experiment, file_name = "/Pole_dist_end.csv") 
data_dist_ <- data_dist |> 
  group_by(scale_r, scale_k8max) |> 
  summarise(max_val = max(BemGEF42 + Cdc42T + BemGEF42_star), 
            r_num = median(r_num), 
            star_max = max(BemGEF42_star)) |> 
  mutate(V = r_num^3 * 4*pi/3)

data_plot_both <- inner_join(data_dist_, data_part) |> 
  filter(scale_k8max == 2)

# Will export and only plot for the stronger feedback 
# value 
p1 <- ggplot(data_plot_both, aes(pole_size, max_val)) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  scale_fill_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_color_manual(values = cbPalette[-1], name = "Strength k8h") + 
  labs(y = "log10(Max Cdc42-GTP concentration) [A.U]", x = "log10(Cdc42-GTP pole area) [µm²]", title = "Model : k8h") + 
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

p2 <- ggplot(data_plot_both, aes(log10(V), log10(pole_size))) + 
  geom_smooth(linewidth=2.0, color="grey50") + 
  geom_point(size=3.0, color="grey10") + 
  scale_fill_manual(values = cbPalette[-1], name = "Strength k8h") + 
  scale_color_manual(values = cbPalette[-1], name = "Strength k8h") + 
  labs(x = "log10(Cell volume) [fL]", y = "log10(Cdc42-GTP pole area) [µm²]", title = "Model : k8h") + 
  scale_y_continuous(limits = c(-0.1, 0.86)) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=12), 
        plot.subtitle = element_text(color="grey30"))

ggsave(file.path(dir_save, "NF_strong_size_max.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(file.path(dir_save, "NF_strong_vol_size.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
