#############################################################################
##
## Shiny app to browse Allen Brain atlas Single-cell RNA Seq mouse datasets 
##
## v0.1 (15.02.21)
##
#############################################################################

library(shiny)
library(tidyverse)
library(ggbeeswarm)
library(scico)
library(glue)
library(cowplot)
source("find_gene.R")


SMART <- list(name = "SMART", 
              cellcount = "76,533")
TenX <- list(name = "10x",
             cellcount ="1,093,785")


# SMART seq
# SMART$expression_medians <- read_csv(file.path("SMART", "medians.csv"))
SMART$expression_means <- read_csv(file.path("SMART", "trimmed_means.csv"))
SMART$metadata <- read_csv(file.path("SMART", "metadata.csv"))
SMART$genes <- unique(SMART$expression_means$feature)


# 10x Seq
# TenX$expression_medians <- read_csv(file.path("10x", "medians.csv"))
TenX$expression_means <- read_csv(file.path("10x", "trimmed_means.csv"))
TenX$metadata <- read.csv(file.path("10x", "metadata.csv"))
TenX$genes <- unique(TenX$expression_means$feature)



## Individual gene
source("plot_individual_gene.R")



# Shiny app
ui = fluidPage(
  selectInput("dataset", "Dataset:",
               c("10x" = "TenX",
                 "SMART" = "SMART")),
  uiOutput("variable"),
  plotOutput("plot"),
  sliderInput("max",
              "Maximal expression value to plot:",
              min = 2,
              max = 20,
              value = 15)
)

server = function(input, output) {
  output$variable <- renderUI({
    selectInput("variable", "Gene:",
              get(input$dataset)$genes)})
  output$plot <-  renderPlot({
    plot_individual_gene(get(input$dataset), input$variable, input$max)
  })
}

 shinyApp(ui, server)
 
 
 
 
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

