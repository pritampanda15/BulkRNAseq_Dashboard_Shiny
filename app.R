#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(DT)

# Set maximum upload size (e.g., 1 GB)
options(shiny.maxRequestSize = 1024 * 1024 * 1024)  # 1 GB

ui <- dashboardPage(
  dashboardHeader(title = "RNA-seq Analysis Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Differential Expression", tabName = "de", icon = icon("chart-line")),
      menuItem("Plots", tabName = "plots", icon = icon("image"))
    )
  ),
  dashboardBody(
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
              fluidRow(
                fileInput("counts_file", "Upload Counts Data (.csv, .tsv, .txt)", 
                          accept = c(".csv", ".tsv", ".txt")),
                fileInput("metadata_file", "Upload Metadata (.csv, .tsv, .txt)", 
                          accept = c(".csv", ".tsv", ".txt")),
                actionButton("analyze", "Run Analysis")
              )
      ),
      # Differential Expression Tab
      tabItem(tabName = "de",
              DTOutput("de_table")
      ),
      # Plots Tab
      tabItem(tabName = "plots",
              fluidRow(
                plotOutput("volcano_plot"),
                plotOutput("heatmap_plot")
              ),
              fluidRow(
                plotOutput("pca_plot"),
                plotOutput("sample_distance_plot")
              ),
              fluidRow(
                plotOutput("ma_plot"),
                plotOutput("dispersion_plot")
              )
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive values to store the uploaded data
  data <- reactiveValues(counts = NULL, metadata = NULL, dds = NULL, results = NULL)
  
  # Helper function to read different file types
  read_data <- function(filepath) {
    ext <- tools::file_ext(filepath)
    if (ext == "csv") {
      return(read.csv(filepath, row.names = 1))
    } else if (ext %in% c("tsv", "txt")) {
      return(read.delim(filepath, row.names = 1))
    } else {
      stop("Unsupported file type: ", ext)
    }
  }
  
  # Load and validate data with progress feedback
  observeEvent(input$analyze, {
    req(input$counts_file, input$metadata_file) # Ensure files are uploaded
    
    withProgress(message = "Running RNA-seq Analysis", value = 0, {
      # Step 1: Read files
      incProgress(0.2, detail = "Reading input files...")
      counts <- read_data(input$counts_file$datapath)
      metadata <- read_data(input$metadata_file$datapath)
      
      # Step 2: Validate data
      incProgress(0.4, detail = "Validating data...")
      if (!all(rownames(metadata) %in% colnames(counts))) {
        showNotification("Metadata row names must match Counts column names", type = "error")
        return()
      }
      
      counts <- counts[ , !(colnames(counts) %in% "Gene Name")]
      
      # Store in reactive values
      data$counts <- counts
      data$metadata <- metadata
      
      # Step 3: Create DESeq2 object
      incProgress(0.6, detail = "Creating DESeq2 dataset...")
      dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = ~ condition  # Use 'sample' as the design factor
      )
      
      # Step 4: Run DESeq2 Analysis
      incProgress(0.8, detail = "Running DESeq2 analysis...")
      dds <- DESeq(dds, fitType = "mean")
      results <- results(dds)
      
      # Step 5: Store Results
      incProgress(1, detail = "Finalizing results...")
      data$dds <- dds
      data$results <- results
    })
    
    showNotification("Analysis complete!", type = "message")
  })
  
  # Render Differential Expression Table
  output$de_table <- DT::renderDT({
    req(data$results)
    results_df <- as.data.frame(data$results)
    DT::datatable(results_df, options = list(pageLength = 10))
  })
  
  # Volcano Plot
  output$volcano_plot <- renderPlot({
    req(data$results)
    res <- as.data.frame(data$results)
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "pvalue",
                    title = "Volcano Plot",
                    pCutoff = 0.05,
                    FCcutoff = 1
    )
  })
  
  # Heatmap Plot
  output$heatmap_plot <- renderPlot({
    req(data$dds)
    vsd <- vst(data$dds, blind = FALSE)
    mat <- assay(vsd)
    mat <- mat[order(rowMeans(mat), decreasing = TRUE), ][1:50, ]  # Top 50 genes
    heatmap(mat, scale = "row", main = "Heatmap of Top Expressed Genes")
  })
  
  # PCA Plot
  output$pca_plot <- renderPlot({
    req(data$dds)
    vsd <- vst(data$dds, blind = FALSE)
    pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_minimal() +
      ggtitle("PCA Plot")
  })
  
  # Sample Distance Heatmap
  output$sample_distance_plot <- renderPlot({
    req(data$dds)
    vsd <- vst(data$dds, blind = FALSE)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(vsd)
    colnames(sampleDistMatrix) <- colnames(vsd)
    
    pheatmap::pheatmap(sampleDistMatrix,
                       clustering_distance_rows = sampleDists,
                       clustering_distance_cols = sampleDists,
                       main = "Sample Distance Heatmap")
  })
  
  # MA Plot
  output$ma_plot <- renderPlot({
    req(data$results)
    plotMA(data$results, main = "MA Plot")
  })
  
  # Dispersion Plot
  output$dispersion_plot <- renderPlot({
    req(data$dds)
    plotDispEsts(data$dds, main = "Dispersion Estimates")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
