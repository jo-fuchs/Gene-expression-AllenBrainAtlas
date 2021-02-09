
## Plotting setup ----
Branchtheme <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    # plot.title.position = "plot", plot.caption.position = "plot",
    legend.position = "none", legend.justification = "top",
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
    # panel.grid.major.x = element_blank(),
    # text = element_text(family = "Source Sans Pro")
  )




## Colors from scico-package --
# Crameri, F. (2018). Scientific colour maps. Zenodo.
# http://doi.org/10.5281/zenodo.1243862

# for Genotype & Location
sci_pal <- "batlow"

# For Branchtypes
sci_pal2 <- "roma"



## plot cell class for selected genes



plot_class <- function(input, gene) {
  input$final %>%
    # mutate(class_label = fct_reorder(class_label, expression)) %>%
    
    ggplot(aes(y = class_label, x = expression, color = class_label, size = number, alpha = log(number))) +
    geom_quasirandom(groupOnX = F, shape = 16) +
    facet_wrap(~feature, nrow = 1) +
    Branchtheme +
    scale_x_continuous(limits = c(0, 10)) +
    scale_size_continuous(range = c(1,4)) +
    labs(
      x = "mRNA abundance [log2(CPM(exons+introns))]",
      y = "",
      title = glue("mRNA levels of {gene}-family in various P40 mouse brain cells"),
      subtitle = glue("single-cell transcriptomes ({input$name}) from multiple cortical areas and the hippocampal formation, including {input$cellcount} total cells")
    ) +
    scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
}



# plot cell types sorted individually
plot_subclass <- function(input, gene) {
  input$final %>%
    filter(feature == gene) %>%
    mutate(subclass_label = fct_reorder(subclass_label, expression)) %>%
    ggplot(aes(y = subclass_label, x = expression)) +
    geom_boxplot(color = "grey") +
    geom_point(shape = 16, aes(color = class_label, size = number, alpha = log(number))) +
    facet_wrap(~feature) +
    Branchtheme +
    labs(
      x = "",
      y = ""
    ) +
    scale_x_continuous(limits = c(0, 10)) +
    scale_size_continuous(range = c(1,4)) +
    scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
}


# # sorting in one function (parallelizable?)
# library(tidytext)
# SMART$final %>%  
#   #mutate(subclass_label = fct_reorder(subclass_label, expression)) %>%
# 
#   ggplot(aes(y = reorder_within(subclass_label, expression, feature), x = expression, color = class_label)) +
#   geom_boxplot() +
#   facet_wrap(~ feature, nrow = 1, scales = "free_y") +
#   Branchtheme +
#   labs(x = "log2(CPM(exons+introns))",
#        y = "") +
#   scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
# 

