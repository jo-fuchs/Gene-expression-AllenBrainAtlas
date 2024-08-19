#############################################################################
##
## Shiny app to browse Allen Brain atlas Single-cell RNA Seq mouse datasets 
##
## v0.2 (21.09.21)
##
#############################################################################

library(shiny)
library(tidyverse)
library(ggbeeswarm)
library(scico)
library(glue)
library(cowplot)
source("find_gene.R")



# SMART seq
SMART <- list(name = "SMART", 
              cellcount = "76,533")

# SMART$expression_medians <- read_csv(file.path("SMART", "medians.csv"))
SMART$expression_means <- vroom::vroom(file.path("SMART", "trimmed_means.csv"), delim = ",")
SMART$metadata <- vroom::vroom(file.path("SMART", "metadata.csv"), delim = ",")
SMART$genes <- unique(SMART$expression_means$feature)


# 10x Seq
TenX <- list(name = "10x",
             cellcount ="1,093,785")

# TenX$expression_medians <- read_csv(file.path("10x", "medians.csv"))
TenX$expression_means <- vroom::vroom(file.path("10x", "trimmed_means.csv"), delim = ",")
TenX$metadata <- vroom::vroom(file.path("10x", "metadata.csv"), delim = ",")
TenX$genes <- unique(TenX$expression_means$feature)


# Adolescent dataset
Adolescent <- list(name = "Adolescent (P20) - mousebrain.org", 
              cellcount = "509,876")

Adolescent$expression_means <- loomR::connect(filename = file.path("Adolescent", "l5_all.agg.loom"), mode = "r+", skip.validate = T)
Adolescent$genes <- Adolescent$expression_means$row.attrs$Gene[]


# Adolescent dataset
Developing <- list(name = "Developing (E10-15) - mousebrain.org", 
                   cellcount = "292,495")

Developing$expression_means <- loomR::connect(filename = file.path("Developing", "dev_all.agg.loom"), mode = "r+", skip.validate = T)
Developing$genes <- Developing$expression_means$row.attrs$Gene[]




## Individual gene
source("plot_individual_gene.R")
source("plot_individual_gene_mousebrain.R")
source("extract_metadata_mousebrain.R")


# Shiny app
ui = fluidPage(
    column(4, selectInput("dataset", "Dataset:",
                 c("Developing mouse brain (10x - mousebrain.org)" = "Developing",
                   "Adolescent mouse brain (10x - mousebrain.org)" = "Adolescent",
                   "Adult mouse brain (10x - Allen brain atlas)" = "TenX",
                   "Adult mouse brain (SMART - Allen brain atlas)" = "SMART"))),
    column(4,uiOutput("variable")),
    column(4,  sliderInput("max",
                           "Maximal expression value to plot:",
                           min = 2,
                           max = 30,
                           value = 25)),
    
  plotOutput("plot", height = 900))


server = function(input, output) {
  output$variable <- renderUI({
    selectInput("variable", "Gene:",
              get(input$dataset)$genes)})
  output$plot <-  renderPlot({
    print(input$dataset)
    if(input$dataset %in% c("SMART", "TenX")) {
      plot_individual_gene(get(input$dataset), input$variable, input$max)
    } else {
      plot_individual_gene_mousebrain(get(input$dataset), input$variable, input$max, Class, TaxonomyRank1)
    }
    
  })
}

 #shinyApp(ui, server)
 runGadget(ui, server, viewer = dialogViewer(width = 1000, height = 900, dialogName = "view"))
 
 
 
 # SMART-specific-App
 
 # ui_SMART = fluidPage(
 #   selectInput("variable", "Gene:",
 #               SMART$genes),
 #   plotOutput("plot"),
 #   sliderInput("max",
 #               "Max value of plots:",
 #               min = 2,
 #               max = 20,
 #               value = 10)
 # )
 # 
 # server_SMART = function(input, output) {
 #   output$plot <-  renderPlot({
 #     plot_individual_gene(SMART, input$variable, input$max)
 #   })
 # }
 # 
 # shinyApp(ui_SMART, server_SMART)
 # 
 # 

# run it as 900 x 500 window
#runGadget(ui, server, viewer = dialogViewer("Gene expression browser", width = 900, height = 500))

