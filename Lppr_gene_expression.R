## Single cell transcriptome dataset - Allen Brain Atlas
##
## downloaded from:
## SMART (2019): 
## https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq
## This data set includes single-cell transcriptomes from
## multiple cortical areas and the hippocampal formation, including 76,533 total
## cells. Samples were collected from dissections of brain regions from 8
## week-old male and female mice, primarily from pan-GABAergic,
## pan-glutamatergic, and pan-neuronal transgenic lines, with the addition of
## more specific transgenic lines and some retrogradely-labeled cells in VISp
## and ALM. 
##
##
## TenX (2020): https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
## This data set includes single-cell transcriptomes from multiple cortical
## areas and the hippocampal formation, including 1,093,785 total cells. Samples
## were collected from dissections of brain regions from ~8 week-old male and
## female mice, from pan-neuronal transgenic lines.
## 
##
## goal: visualize PRG-Expression over cell types



library(tidyverse)
library(ggbeeswarm)
library(scico)
library(glue)
library(cowplot)
source("find_gene.R")
source("plotting_setup.R")

SMART <- list(name = "SMART", 
              cellcount = "76,533")
TenX <- list(name = "10x",
             cellcount ="1,093,785")


# SMART seq
SMART$expression_medians <- read_csv(file.path("SMART", "medians.csv"))
SMART$expression_means <- read_csv(file.path("SMART", "trimmed_means.csv"))
SMART$metadata <- read_csv(file.path("SMART", "metadata.csv"))


# 10x Seq
TenX$expression_medians <- read_csv(file.path("10x", "medians.csv"))
TenX$expression_means <- read_csv(file.path("10x", "trimmed_means.csv"))
TenX$metadata <- read.csv(file.path("10x", "metadata.csv"))



# Decide for which gene (family) to analyse 

SMART <- find_gene(SMART, "Lppr")
TenX <- find_gene(TenX, "Plppr")



### Plots
SMART$Figure_A <- plot_class(SMART, "Lppr")
TenX$Figure_A <- plot_class(TenX, "Plppr")


# Each individual gene has to be named separately here
SMART$Figure_B1 <- plot_subclass(SMART, "Lppr1")
SMART$Figure_B2 <- plot_subclass(SMART, "Lppr2")
SMART$Figure_B3 <- plot_subclass(SMART, "Lppr3")
SMART$Figure_B4 <- plot_subclass(SMART, "Lppr4")
SMART$Figure_B5 <- plot_subclass(SMART, "Lppr5")

TenX$Figure_B1 <- plot_subclass(TenX, "Plppr1")
TenX$Figure_B2 <- plot_subclass(TenX, "Plppr2")
TenX$Figure_B3 <- plot_subclass(TenX, "Plppr3")
TenX$Figure_B4 <- plot_subclass(TenX, "Plppr4")
TenX$Figure_B5 <- plot_subclass(TenX, "Plppr5")


# Finish up SMART Figure
SMART$Figure_B <- plot_grid(SMART$Figure_B1, SMART$Figure_B2, SMART$Figure_B3, SMART$Figure_B4, SMART$Figure_B5,
  scale = 1,
  nrow = 1, align = "h", axis = "bt"
)


SMART$Figure <- plot_grid(SMART$Figure_A, SMART$Figure_B, nrow = 2, rel_heights = c(1, 3), labels = c("A", "B"))


ggsave("Lppr_means_SMART.png", SMART$Figure, device = "png", scale = 1.5, width = 210, height = 180, units = "mm")



# Finish up 10x Figure
TenX$Figure_B <- plot_grid(TenX$Figure_B1, TenX$Figure_B2, TenX$Figure_B3, TenX$Figure_B4, TenX$Figure_B5,
                     scale = 1,
                     nrow = 1, align = "h", axis = "bt"
)


TenX$Figure <- plot_grid(TenX$Figure_A, TenX$Figure_B, nrow = 2, rel_heights = c(1, 3), labels = c("A", "B"))


ggsave("Lppr_means_TenX.png", TenX$Figure, device = "png", scale = 1.5, width = 210, height = 180, units = "mm")
