# Plotting everything for only one gene
plot_individual_gene <- function(input, gene, max_val){
input <- find_gene(input, gene)
input$Figure_A <- plot_class(input, gene, max_val)
input$Figure_B <- plot_neighborhood(input, gene, max_val)
input$Figure_C <- plot_subclass(input, gene, max_val)


input$Header <- ggplot() + Branchtheme + 
  labs(title = glue("mRNA levels of {gene} in various P50-P60 mouse brain cells"),
       subtitle = glue("single-cell transcriptomes ({input$name}) from multiple cortical areas and the hippocampal formation, including {input$cellcount} total cells"))



input$Figure_left <- plot_grid(input$Figure_A, input$Figure_B, nrow = 2, rel_heights = c(1, 1), labels = c("A", "B"))

input$Figure_data <- plot_grid(input$Figure_left, input$Figure_C, ncol = 2, labels = c("", "C"))

input$Figure <- plot_grid(input$Header, input$Figure_data, nrow = 2, rel_heights = c(1, 10), labels = c("", ""))


# ggsave(glue("{gene}_expression_{input$name}.png"), input$Figure, device = "png", scale = 1.5, width = 210, height = 150, units = "mm")

# ggsave(glue("{gene}_expression_{input$name}.pdf"), input$Figure, device = "pdf", scale = 1.5, width = 210, height = 150, units = "mm")
return(input$Figure)
}

## Colors from scico-package --
# Crameri, F. (2018). Scientific colour maps. Zenodo.
# http://doi.org/10.5281/zenodo.1243862

# for Genotype & Location
sci_pal <- "batlow"



## Plotting functions individual stuff
Branchtheme <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),  
    strip.background = element_blank(), strip.text.x = element_blank(),
    # plot.title.position = "plot", plot.caption.position = "plot",
    legend.position = "none", legend.justification = "top",
    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()
    # panel.grid.major.x = element_blank(),
    # text = element_text(family = "Source Sans Pro")
  )



plot_class <- function(input, gene, max_val) {
  input$final %>%
    # mutate(class_label = fct_reorder(class_label, expression)) %>%
    
    ggplot(aes(y = class_label, x = expression, color = class_label, size = number, alpha = log(number))) +
    geom_quasirandom(groupOnX = F, shape = 16) +
   # facet_wrap(~feature, nrow = 1) +
    Branchtheme +
    coord_cartesian(xlim = c(0,max_val)) +
    #scale_x_continuous(limits = c(0, 10)) +
    scale_size_continuous(range = c(1,3)) +
    labs(
      x = "mRNA abundance [log2(CPM(exons+introns))]",
      y = "",
      title = glue("{gene}-expression in general brain cell types")
    ) +
    scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
}



# plot neighbourhoods sorted individually
plot_neighborhood <- function(input, gene, max_val) {
  input$final %>%
    filter(feature == gene) %>%
    mutate(neighborhood_label = fct_reorder(neighborhood_label, expression)) %>%
    ggplot(aes(y = neighborhood_label, x = expression)) +
    geom_boxplot(color = "grey") +
    geom_point(shape = 16, aes(color = class_label, size = number, alpha = log(number))) +
   # facet_wrap(~feature) +
    Branchtheme +
    labs(
      x = "mRNA abundance [log2(CPM(exons+introns))]",
      y = "",
      title = glue("{gene}-expression in neighborhoods")
    ) +
    coord_cartesian(xlim = c(0,max_val)) +
#    scale_x_continuous(limits = c(0, 10)) +
    scale_size_continuous(range = c(1,3)) +
    scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
}

# plot cell types sorted individually
plot_subclass <- function(input, gene, max_val) {
  input$final %>%
    filter(feature == gene) %>%
    mutate(subclass_label = fct_reorder(subclass_label, expression)) %>%
    ggplot(aes(y = subclass_label, x = expression)) +
    geom_boxplot(color = "grey") +
    geom_point(shape = 16, aes(color = class_label, size = number, alpha = log(number))) +
   # facet_wrap(~feature) +
    Branchtheme +
    labs(
      x = "mRNA abundance [log2(CPM(exons+introns))]",
      y = "",
      title = glue("{gene}-expression in individual cell types")
    ) +
    coord_cartesian(xlim = c(0,max_val)) +
#    scale_x_continuous(limits = c(0, 10)) +
    scale_size_continuous(range = c(1,4)) +
    scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0)
}

