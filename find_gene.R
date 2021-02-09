

find_gene <- function(input, gene) {
  # # replace spaces with . in metadata (only necessary with because base R "corrects" column names)
  # input$metadata <- input$metadata %>% mutate(cluster_label = str_replace_all(cluster_label," ", "."),
  #                                 cluster_label = str_replace_all(cluster_label,"-", "."),
  #                                 cluster_label = str_replace_all(cluster_label,"/", "."))
  
  
  # only keep important input columns
  # input$clusters <- input$metadata %>%
  #   select(c("cluster_label", "class_label", "subclass_label")) %>%
  #   distinct()
  
  # select Gene of interest from expression levels
  input$gene_expression <- input$expression_means %>% filter(str_detect(feature, gene))
  
  
  
  # count individual cell types
  input$clusters <- input$metadata %>%
    group_by(cluster_label, class_label, subclass_label) %>%
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
    filter(
      class_label != is.na(class_label),
      feature != is.na(feature)
    ) %>%
    mutate(class_label = factor(class_label, levels = c("Non-Neuronal", "GABAergic", "Glutamatergic")))
  
  return(input)
}

