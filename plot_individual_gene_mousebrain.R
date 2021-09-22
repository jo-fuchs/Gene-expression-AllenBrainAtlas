
## Expression 
plot_individual_gene_mousebrain <- function(dataset, gene, max_val, y_axis , color) {
  
# Colors
  ## blue-glut: 09215f
  ## green-gaba: 4e734d
  ## orange-rest: f3a57c
  cols <- c("Glutamate" = "#09215f", "GABA" = "#4e734d", "Unclear" = "#f3a57c",
            "Acetylcholine" = "#3b8db7", "Serotonin" = "#9d9548", "Dopamine" = "#752d16", 
            "Nitric oxide" = "#929292", "Noradrenaline" = "#fe9ac5")
  
plot_expression_mousebrain <- function(dataset, gene, max_val, y_axis , color) {
  #print(dataset  )
  if (gene %in% dataset$genes){
    
    expression1 <- (dataset$expression_means[["matrix"]][,dataset$expression_means$row.attrs$Gene[] == gene])*1000000 + 1
   
    # merge with metadata
    metadata <- extract_metadata_mousebrain(dataset)
    
    total <- tibble(expression1, metadata)
    
    
    if(paste(deparse(substitute(y_axis))) %in% names(total)) {
      
      # if(paste(deparse(substitute(color))) == "Transmitter") {
      #   total %>% mutate(Y_axis = fct_reorder({{y_axis}}, expression1, median)) %>% 
      #     ggplot(aes(x = log2(expression1), y = Y_axis, col = {{color}}, size = NCells)) + 
      #     geom_boxplot(fill = NA, color = "grey", size = 0.5) +
      #     geom_quasirandom(alpha = 0.5, groupOnX = F) + 
      #     scale_x_continuous(limits = c(0, max_val)) +
      #     labs(title = paste("Expression of ", gene, "in adolescent mouse brain"),
      #          x = "log10 expression value",
      #          y = "",
      #          subtitle = paste(deparse(substitute(y_axis)))) + 
      #     scale_color_brewer(palette = "Set1", type = "div") +
      #     Branchtheme
      # }
      # else  {
        total %>% mutate(Y_axis = fct_reorder({{y_axis}}, log2(expression1), median)) %>% 
          ggplot(aes(x = log2(expression1), y = Y_axis, col = {{color}}, size = NCells, alpha = log(NCells))) + 
          geom_boxplot(fill = NA, color = "grey", size = 0.5) +
          geom_quasirandom(groupOnX = F, shape = 16) + 
          scale_x_continuous(limits = c(0, max_val)) +
          labs(x = "log2 expression value per million",
               y = "",
               subtitle = paste(deparse(substitute(y_axis)))) + 
          scale_size_continuous(range = c(1,3)) +
          scale_colour_manual(values = cols) +
          #scale_color_scico_d(palette = sci_pal, begin = 0.75, end = 0) +
          Branchtheme
      #}
    }
    else print(paste("Column not found as metadata"))
    
  }
  
  else print(paste("Gene not found in dataset. Did you mean one of the following: ", paste(genes[grep(gene, genes)], collapse = ", "),", or ",  paste(head(genes[agrep(gene, genes)], 10), collapse = ", " ), "...", sep = ""))
}

input <- list()

if(dataset$name == Adolescent$name){
  input$Figure_A <- plot_expression_mousebrain(dataset, gene, max_val, Transmitter, Transmitter)
  input$Figure_B <- plot_expression_mousebrain(dataset, gene, max_val, TaxonomyRank2, Transmitter)
  input$Figure_C <- plot_expression_mousebrain(dataset, gene, max_val, TaxonomyRank4, Transmitter)
  
  
  input$Header <- ggplot() + Branchtheme + 
    labs(title = glue("mRNA levels of {gene} in various P20 mouse brain cells"),
         subtitle = glue("single-cell transcriptomes ({dataset$name}) from multiple brain regions, including {dataset$cellcount} total cells"))
  
} else {
  input$Figure_A <- plot_expression_mousebrain(dataset, gene, max_val, Transmitter, Transmitter )
  input$Figure_B <- plot_expression_mousebrain(dataset, gene, max_val, Class, Transmitter)
  input$Figure_C <- plot_expression_mousebrain(dataset, gene, max_val, Subclass, Transmitter) + theme(axis.text.y = element_text(size = 7))
  
  
  input$Header <- ggplot() + Branchtheme + 
    labs(title = glue("mRNA levels of {gene} in various E10 mouse brain cells"),
         subtitle = glue("single-cell transcriptomes ({dataset$name}) from multiple brain regions, including {dataset$cellcount} total cells"))
  
} 
  
  input$Figure_left <- plot_grid(input$Figure_A, input$Figure_B, nrow = 2, align = "v", rel_heights = c(1, 2), labels = c("A", "B"))
  
  input$Figure_data <- plot_grid(input$Figure_left, input$Figure_C, ncol = 2, labels = c("", "C"))
  
  input$Figure <- plot_grid(input$Header, input$Figure_data, nrow = 2, rel_heights = c(1, 10), labels = c("", ""))
  
  
  # ggsave(glue("{gene}_expression_{input$name}.png"), input$Figure, device = "png", scale = 1.5, width = 210, height = 150, units = "mm")
  
  # ggsave(glue("{gene}_expression_{input$name}.pdf"), input$Figure, device = "pdf", scale = 1.5, width = 210, height = 150, units = "mm")
  return(input$Figure)

}

