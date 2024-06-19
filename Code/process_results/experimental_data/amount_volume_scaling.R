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
# Amount volume experimental data from the Skottheim lab (Fig. 2)
# -----------------------------------------------------------------------------------------------------------------
path_data <- file.path("..", "..", "..", "Data", "Table_S3_Yeast.xlsx")
gene_list <- c("CDC24", "BUD3", "BEM2", "BEM3", "RGA1", "RGA2", "CLA4", "STE20", "BEM1", "GIC1", "GIC2", "AXL2", "RSR1", "CDC42", "CDC11", "CDC12", "CDC3", "CDC10")
gene_class <- tibble(Gene_name = gene_list, 
                     classification = c("GEF", "GEF", "GAP", "GAP", "GAP", "GAP", "Effector", "Effector", "Effector", "Effector", "Effector", "Effector", "Effector", "Cdc42", "Septin", "Septin", "Septin", "Septin"))

data_G1_ethanol_ <- read_excel(path_data, sheet = "G1_arrest_time_Ethanol")
data_G1_ethanol <- data_G1_ethanol_ |> 
  rename("Gene_name" = "Gene name") |> 
  filter(Gene_name %in% gene_class$Gene_name) |> 
  inner_join(gene_class, by = c("Gene_name")) |> 
  rename("Protein_slope_Exp1" = "Protein_slope_Exp#1", 
         "Protein_slope_Exp2" = "Protein_slope_Exp#2", 
         "Protein_slope_Exp3" = "Protein_slope_Exp#3") |> 
  mutate(tag = "G1_arrest_ethanol")

data_G1_glucose_ <- read_excel(path_data, sheet = "G1_arrest_time_Glucose")
data_G1_glucose <- data_G1_glucose_ |> 
  rename("Gene_name" = "Gene name") |> 
  filter(Gene_name %in% gene_class$Gene_name) |> 
  inner_join(gene_class, by = c("Gene_name")) |> 
  rename("Protein_slope_Exp1" = "Protein_slope_Exp#1", 
         "Protein_slope_Exp2" = "Protein_slope_Exp#2", 
         "Protein_slope_Exp3" = "Protein_slope_Exp#3") |> 
  mutate(tag = "G1_arrest_glucose")

data_plot <- bind_rows(data_G1_ethanol, data_G1_glucose)
data_plot <- data_plot |> mutate(Gene_name = factor(Gene_name, levels = gene_class$Gene_name))

p <- ggplot(data_plot, aes(color = classification)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(1.0, -1.0), linetype=2) +
  geom_point(aes(Gene_name, Protein_slope_Exp1), size=3.0) +
  geom_point(aes(Gene_name, Protein_slope_Exp2), size=3.0) +
  geom_point(aes(Gene_name, Protein_slope_Exp3), size=3.0) +
  geom_crossbar(aes(x=Gene_name, y=Mean_protein_slope, ymin = Mean_protein_slope, ymax = Mean_protein_slope), linewidth=1.0) +
  facet_wrap(~tag) +
  labs(x = "Gene name", y = "Protein slope", 
       title = "Most polarisation related genes have close to constant concentration (slope = 0)", 
       subtitle = "The bar is the mean from fitting three replication experiments") + 
  scale_y_continuous(breaks = seq(from = -1.25, to = 1.25, by = 0.25)) +
  scale_color_manual(values = cbPalette[-c(1, 5)],  name = "Type") + 
  theme_bw(base_size = 14) +
  theme(panel.grid.minor.y = element_blank(), 
        legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"), 
        axis.text.x = element_text(color="grey10", angle=45, size=17, hjust=1, vjust=1))

data_mean <- data_plot |> 
  group_by(Protein, tag) |> 
  summarise(mean_mean = mean(Mean_protein_slope)) |> 
  filter(Protein %in% c("CDC42_YEAST", "CDC24_YEAST", "BEM2_YEAST", "BEM1_YEAST"))

dir_save <- file.path("..", "..", "..", "Results", "Pole_size", "Amount_volume")
if( !dir.exists(dir_save) ) dir.create(dir_save, recursive = T)
ggsave(file.path(dir_save, "Amount_vol_scaling.svg"), p, width = BASE_WIDTH*2.9, height = BASE_HEIGHT*1.5)  

# Last data on size mutants 
data_size_mutants_ <- read_excel(path_data, sheet = "Size_mutants")
data_size_mutants <- data_size_mutants_ |> 
  rename("Gene_name" = "Gene name") |> 
  filter(Gene_name %in% gene_class$Gene_name) |> 
  inner_join(gene_class, by = c("Gene_name")) |> 
  rename("Protein_slope_Exp1" = "Protein_slope_Exp#1", 
         "Protein_slope_Exp2" = "Protein_slope_Exp#2", 
         "Protein_slope_Exp3" = "Protein_slope_Exp#3") |> 
  mutate(tag = "Size_mutants")

data_size_mutants <- data_size_mutants |> mutate(Gene_name = factor(Gene_name, levels = gene_class$Gene_name))

p2 <- ggplot(data_size_mutants, aes(color = classification)) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(1.0, -1.0), linetype=2) +
  geom_point(aes(Gene_name, Protein_slope_Exp1), size=3.0) +
  geom_point(aes(Gene_name, Protein_slope_Exp2), size=3.0) +
  geom_point(aes(Gene_name, Protein_slope_Exp3), size=3.0) +
  geom_crossbar(aes(x=Gene_name, y=Mean_protein_slope, ymin = Mean_protein_slope, ymax = Mean_protein_slope), linewidth=1.0) +
  facet_wrap(~tag) +
  labs(x = "Gene name", y = "Protein slope", 
       title = "Most polarisation related genes have close to constant concentration (slope = 0)", 
       subtitle = "The bar is the mean from fitting three replication experiments") + 
  scale_y_continuous(breaks = seq(from = -1.25, to = 1.25, by = 0.25), limits=c(-1.3, 1.3)) +
  scale_color_manual(values = cbPalette[-c(1, 5)],  name = "Type") + 
  theme_bw(base_size = 14) +
  theme(panel.grid.minor.y = element_blank(), 
        legend.position = "bottom", 
        axis.title = element_text(color="grey10"),
        plot.title = element_text(color="grey10", face ="bold", size=15), 
        plot.subtitle = element_text(color="grey30"), 
        axis.text.x = element_text(color="grey10", angle=45, size=17, hjust=1, vjust=1))

ggsave(file.path(dir_save, "Amount_vol_scaling_mutant.svg"), p2, width = BASE_WIDTH*1.3, height = BASE_HEIGHT*1.5, dpi=300)  
