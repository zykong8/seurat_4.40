library(shiny)

server <- function(input, output, session){
  obj <- reactive(
    {
      readRDS("/Users/xiaofei/Desktop/SeuratApp/app/data/pbmc_tutorial.rds")
    }
  )

  pcas <- reactive({
    input$pcas
  })
  
  res <- reactive({
    input$resolution
  })
  
  newObj <- reactive({
    UpdateObj(pcas(), res(), obj())
  })
  
  output$umap <- renderPlot({
    DimPlot(newObj(), reduction = "umap", label = TRUE) + theme_dr()
  })
    
  
  output$tsne <- renderPlot({
    DimPlot(newObj(), reduction = "tsne", label = TRUE) + theme_dr()
  })
  
  # All genes
  allGenes <- reactive({
    rownames(newObj())
  })
  
  # Feature plot
  output$featureplotUMAP <- renderPlot({
    req(input$genelist)
    if (input$genelist %in% allGenes()){
      FeaturePlot(newObj(), features = input$genelist, reduction = "umap") + theme_dr()
    }
  })
  
  output$featureplotTSNE <- renderPlot({
    req(input$genelist)
    if (input$genelist %in% allGenes()){
      FeaturePlot(newObj(), features = input$genelist, reduction = "tsne") + theme_dr()
    }
  })
  
  # Dotplot
  dotgenelist <- reactive({
    req(input$dotgenelist)
    strsplit(input$dotgenelist, split = ",|\\s+", fixed = FALSE)[[1]]
  })
  
  output$dotplot <- renderPlot({
    req(input$dotgenelist)
    if (length(intersect(unique(dotgenelist()), unique(allGenes()))) >0){
      my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                     '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
                     '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
                     '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                     '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
                     '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
      
      Clustered_DotPlot(seurat_object = newObj(), 
                        features = intersect(unique(dotgenelist()), unique(allGenes())), 
                        colors_use_exp = c('#330066','#336699','#66CC66','#FFCC33'),
                        colors_use_idents = my36colors,
                        flip = TRUE, cluster_ident = FALSE, plot_km_elbow = FALSE,
                        x_lab_rotate = 90, row_label_size = 16)
    }
  })
  
  # Heatmap
  output$htmap <- renderPlot({
    req(input$genelist)
    DoHeatmap(object = newObj(), features = intersect(unique(dotgenelist()), unique(allGenes()))) +
      scale_fill_gradientn(colors = c("navy","white","firebrick3"))
  })

  
  
}


# Start Server
shinyServer(server)