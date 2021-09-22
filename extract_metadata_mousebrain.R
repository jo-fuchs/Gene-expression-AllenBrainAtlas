extract_metadata_mousebrain <- function(dataset) {
  
  working_data = dataset$expression_means
  
  # get column names
  metadata_headers <- names(working_data$col.attrs[])
  metadata_headers <- metadata_headers[str_detect(metadata_headers, "SampleID|Age|Tissue|Sex|Bucket|Comment|Marker", negate = T)] # remove unnecessary cols
  
  
  # create metadata table
  metadata = tibble(index = 1:length(working_data$col.attrs$Class[]))
  for (col in metadata_headers){
    metadata[, col] <-working_data$col.attrs[[col]][]
  }
  
  
  # find neurotransmitter info
  extract_neurotransmitter <- function(input, searchcolumn){
    input %>% mutate(Transmitter = case_when(
      str_detect({{searchcolumn}}, "Serotonin") ~ "Serotonin",
      str_detect({{searchcolumn}}, "Acetylcholine") ~ "Acetylcholine",
      str_detect({{searchcolumn}}, "Dopamine") ~ "Dopamine",
      str_detect({{searchcolumn}}, "Nitric") ~ "Nitric oxide",
      str_detect({{searchcolumn}}, "Noradren")~ "Noradrenaline",
      str_detect({{searchcolumn}}, "lutamat")~ "Glutamate",
      str_detect({{searchcolumn}}, "GABA") ~ "GABA",
      
      TRUE ~ "Unclear" ),
      Transmitter <- factor(Transmitter, 
                            levels = c("Unclear", "GABA", "Glutamate", "ACh", "Dopamine", "Serotonin", "Noradrenaline", "Nitric oxide"))
      
    ) # the rest
  }
  
  if(dataset$name == "Adolescent (P20) - mousebrain.org") {
    metadata <- extract_neurotransmitter(metadata, Neurotransmitter)
    
  } else {
    metadata <- extract_neurotransmitter(metadata, Subclass)
  }
  return(metadata)
  
  
}
