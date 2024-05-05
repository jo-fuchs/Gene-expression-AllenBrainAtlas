## individual cells correlation

library(tidyverse)
library(scico)
library(ggbeeswarm)

## loading data ####

SMART <- list(name = "SMART", 
              cellcount = "76,533")

## load metadata

SMART$metadata_exons <- vroom::vroom(file.path("SMART", "2024",  "metadata.csv"), delim = ",") |>
  select(sample_name, cluster_label, class_label, subclass_label, neighborhood_label, donor_sex_label) |> 
  mutate(class_label = factor(class_label, levels = c("Non-Neuronal", "GABAergic", "Glutamatergic")))

# store gene names
SMART$genes <- unique(SMART$expression_means$feature)


# ## load exon-only sequencing data
# 
# # initial filtering, only done once per selection -> write to disk to save operation time later
# library(chunked)
# SMART$individual_exons <- read_csv_chunkwise(file.path("SMART", "expression_matrix_exons.csv"), chunk_size = 10000) %>%
#   select(sample_name, Lppr3, Lppr1, Lppr2, Lppr4, Lppr5, Pten, Rdx, Ntrk1, Ntrk2, Ntrk3, Insr, Egfr, Inpp5d, Inppl1, Akt1, Akt2, Akt3, Ppap2a, Ppap2b, Ppap2c, Ppapdc2, Ppapdc3, Ppapdc1a, Ppapdc1b, Cnr1, Lpar1, Lpar2, Lpar3, Lpar4, Lpar5, Lpar6)
# # sample name has to be in there, the other entries are gene names of interest
# # (make sure they are spelled exactly as in dataset using Shiny app)
# 
# write_chunkwise(SMART$individual_exons, file = file.path("SMART", "Lppr_friends_exon_matrix.csv"))


SMART$individual_exons <- vroom::vroom(file.path("SMART", "2024", "Lppr_friends_exon_matrix.csv"), delim = ",") |> 
  mutate(         
    # zero-expression levels conflict with log calculations -> add 1 to every expression level
    across(where(is.numeric), ~.x+1))

SMART$individual_cleaned <- SMART$individual_exons |> 
  left_join(SMART$metadata_exons, by = "sample_name") |> 
  filter(subclass_label != "V3d" & !is.na(subclass_label)) |> # only one cell
  group_by(subclass_label) |> 
  mutate(number_cells = n(),
         subclass_label = factor(subclass_label)) |> 
  ungroup()


SMART$exon_expression <- SMART$individual_cleaned |> 
  pivot_longer(cols =  !c("sample_name", "cluster_label", "class_label", "subclass_label", "neighborhood_label", "donor_sex_label"), names_to = "feature", values_to = "expression") |> 
  group_by(cluster_label, class_label, subclass_label, neighborhood_label, feature) |> 
  summarise(expression = (mean(log2(expression), trim = 0.25)),
            number = n() ) |> 
  ungroup()


## Plotting setup

sci_pal <- "batlow"

Branchtheme <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),  
    strip.background = element_blank(), 
    plot.title.position = "plot", plot.caption.position = "plot",
    legend.justification = "top",legend.position = "none",
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
    text = element_text(size = 11)
    )


LPPR_geneset = c ("Lppr1", "Lppr2", "Lppr3", "Lppr4", "Lppr5")


# Figure_3D exon_only
(Figure_3D_exon <- SMART$exon_expression |> 
  filter(feature %in% LPPR_geneset) |> 
  ggplot(aes(y = class_label, x = expression, color = class_label, size = number, alpha = log(number))) +
  geom_quasirandom(groupOnX = F, shape = 16) +
  facet_wrap(~feature, nrow = 1) +
  Branchtheme +
  scale_x_continuous(limits = c(0, 10)) +
  scale_size_continuous(range = c(1,3)) +
  labs(
    x = "mRNA abundance [log2(CPM(exons+1))]",
    y = ""  ) +
  scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0))

  ggsave("Figure_3D_exon.png", Figure_3D_exon, device = "png", width = 37, height = 7, units = "cm", dpi = 300, bg = "white")



