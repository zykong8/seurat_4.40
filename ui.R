library(shiny)
library(dplyr)
library(patchwork)
library(Seurat)
library(tidydr)
library(scCustomize)
library(ggplot2)


ui <- fluidPage(
  h1("Single Cell Analysis"),
  fluidRow(
    column(
      width = 2,
      h2("Re-clustering"),
      sliderInput("pcas", label = "Determine the ‘dimensionality’ of the dataset",
                  min = 10, max = 80, value = 15, step = 5),
      sliderInput("resolution", label = "Determine the resolution of clustering cells", 
                  min = 0.4, max = 1.2, value = 0.6, step = 0.2)
    ),
    column(
      width = 5,
      h2("UMAP"),
      plotOutput("umap")
    ),
    column(
      width = 5,
      h2("tSNE"),
      plotOutput("tsne")
    )
  ),
  
  h1("Feature Plots"),
  fluidRow(
    column(
      width = 2,
      h2("Select Genes"),
      textInput("genelist", label = "Please enter genes:", value = NULL)
    ),
    column(
      width = 5,
      h2("UMAP"),
      plotOutput("featureplotUMAP")
    ),
    column(
      width = 5,
      h2("t-SNE"),
      plotOutput("featureplotTSNE")
    ),
    
  ),
  
  h1("Clustered Dotplot"),
  fluidRow(
    column(
      width = 2,
      h2("Select Genes"),
      textAreaInput("dotgenelist", label = "Please enter gene-list:", value = NULL,
                    width = 300, height = 200)
    ),
    column(
      width = 5,
      h2("Dotplot"),
      plotOutput("dotplot")
    ),
    column(
      width = 5,
      h2("Heatmap"),
      plotOutput("htmap")
    )
  )
  
  
  
)

# Start UI
shinyUI(ui)