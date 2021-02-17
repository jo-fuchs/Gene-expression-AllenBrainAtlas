
find_gene <- function(input, gene) {
  # # replace spaces with . in metadata (only necessary with because base R "corrects" column names)
  # input$metadata <- input$metadata %>% mutate(cluster_label = str_replace_all(cluster_label," ", "."),
  #                                 cluster_label = str_replace_all(cluster_label,"-", "."),
  #                                 cluster_label = str_replace_all(cluster_label,"/", "."))
  
  
  # only keep important input columns
  # input$clusters <- input$metadata %>%
  #   select(c("cluster_label", "class_label", "subclass_label")) %>%
  #   distinct()
  
  # lookup table for more informative neighborhood names
  neighborhoods <- tribble(~original_name, ~recoded_name,
                           "PT", "Pyramidal tract",
                           "CGE", "Caudal ganglionic eminence",
                           "L2_3_IT", "Layer 2/3 intra-telencephalic",
                           "MGE", "Medial ganglionic eminence",
                           "Other", "Other",
                           "NP_L6CT_L6b", "Near-projecting & Layer 6 & 6b",
                           "L4_5_6_IT_Car3", "Layer 4/5/6 intra-telencephalic",
                           "DG_SUB_CA" ,"Hippocampus & adjacent regions"
  )
  
  # select Gene of interest from expression levels
  input$gene_expression <- input$expression_means %>% filter(str_detect(feature, gene))
  
  
  
  # count individual cell types
  input$clusters <- input$metadata %>%
    group_by(cluster_label, class_label, subclass_label, neighborhood_label) %>%
    summarise(number = n())
  
  
  input$gene_expression_l <- input$gene_expression %>%
    pivot_longer(!feature,
                 names_to = "cluster_label",
                 # names_pattern = "X(.+)",
                 values_to = "expression"
    )
  
  
  # join with clusters
  input$final <- input$gene_expression_l %>%
    full_join(input$clusters, by = "cluster_label") %>%
    left_join(neighborhoods, by = c("neighborhood_label" = "original_name")) %>%  
    filter(
      class_label != is.na(class_label),
      feature != is.na(feature)
    ) %>%
    mutate(class_label = factor(class_label, levels = c("Non-Neuronal", "GABAergic", "Glutamatergic")))
  
  return(input)
}



