library(shiny)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(EnhancedVolcano)
library(DT)

options(shiny.maxRequestSize = 100 * 1024^2)

ui <- fluidPage(
  titlePanel("Differential Expression Analysis with DESeq2"),
  sidebarLayout(
    sidebarPanel(
      fileInput("count_file", "Upload Count Data", accept = c(".csv", ".txt")),
      fileInput("condition_file", "Upload Condition Data", accept = c(".csv", ".txt")),
      actionButton("run_analysis", "Run DESeq2 Analysis"),
      tags$hr(),
      h4("Processing Log"),
      verbatimTextOutput("log_output")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview",
                 h4("Count Data Preview"),
                 DTOutput("count_table"),
                 h4("Condition Data Preview"),
                 DTOutput("condition_table")
        ),
        tabPanel("MA Plot", 
                 plotOutput("maplot"),
                 downloadButton("download_ma", "Download MA Plot")
        ),
        tabPanel("Boxplot", 
                 plotOutput("boxplot"),
                 downloadButton("download_box", "Download Boxplot")
        ),
        tabPanel("Scatter Plot", 
                 plotOutput("scatterplot"),
                 downloadButton("download_scatter", "Download Scatter Plot")
        ),
        tabPanel("Volcano Plot", 
                 plotOutput("volcanoplot"),
                 downloadButton("download_volcano", "Download Volcano Plot")
        ),
        tabPanel("Heatmap", 
                 plotOutput("heatmap"),
                 downloadButton("download_heatmap", "Download Heatmap")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  count_data <- reactive({
    req(input$count_file)
    tryCatch({
      if (grepl("\\.csv$", input$count_file$name)) {
        read.csv(input$count_file$datapath, row.names = 1)
      } else {
        read.delim(input$count_file$datapath, row.names = 1)
      }
    }, error = function(e) {
      showNotification("Error reading count file!", type = "error")
      return(NULL)
    })
  })
  
  condition_data <- reactive({
    req(input$condition_file)
    tryCatch({
      if (grepl("\\.csv$", input$condition_file$name)) {
        read.csv(input$condition_file$datapath)
      } else {
        read.delim(input$condition_file$datapath)
      }
    }, error = function(e) {
      showNotification("Error reading condition file!", type = "error")
      return(NULL)
    })
  })
  
  output$count_table <- renderDT({ datatable(count_data(), options = list(scrollX = TRUE)) })
  output$condition_table <- renderDT({ datatable(condition_data(), options = list(scrollX = TRUE)) })
  
  analysis_results <- eventReactive(input$run_analysis, {
    req(count_data(), condition_data())
    
    withProgress(message = 'Processing', value = 0, {
      incProgress(0.1, detail = "Preparing data...")
      
      counts <- count_data()
      cond_data <- condition_data()
      condition <- as.character(cond_data$Condition)
      coldata <- data.frame(row.names = cond_data$Sample, condition)
      
      incProgress(0.3, detail = "Creating DESeq2 object...")
      dds <- DESeqDataSetFromMatrix(
        countData = counts[rowSums(counts) > 0, ],
        colData = coldata,
        design = ~condition
      )
      
      incProgress(0.5, detail = "Running DESeq2...")
      ddsDE <- DESeq(dds)
      
      incProgress(0.8, detail = "Extracting results...")
      res <- results(ddsDE, independentFiltering = TRUE)
      normcounts <- counts(ddsDE, normalized = TRUE)
      
      incProgress(1.0, detail = "Complete!")
      
      return(list(ddsDE = ddsDE, results = res, normcounts = normcounts, coldata = coldata))
    })
  })
  
  output$log_output <- renderText({
    analysis_results()
    "Analysis completed successfully!"
  })
  
  output$maplot <- renderPlot({
    res <- analysis_results()$results
    if (!inherits(res, "DESeqResults")) return(NULL)
    DESeq2::plotMA(res, main = "MA Plot", ylim = c(-5, 5))
  })
  
  output$boxplot <- renderPlot({
    dds <- analysis_results()$ddsDE
    if (is.null(dds)) return(NULL)
    
    boxplot(log2(counts(dds, normalized = TRUE) + 1),
            col = rainbow(ncol(dds)), 
            outline = FALSE,
            main = "Boxplot of Normalized Counts",
            xlab = "Samples", 
            ylab = "log2(normalized counts + 1)")
  })
  
  output$scatterplot <- renderPlot({
    res <- as.data.frame(analysis_results()$results)
    normcounts <- analysis_results()$normcounts
    coldata <- analysis_results()$coldata
    
    if (length(unique(coldata$condition)) < 2) {
      showNotification("Error: Not enough conditions for scatter plot!", type = "error")
      return(NULL)
    }
    
    condition_levels <- levels(factor(coldata$condition))
    res$baseMeanA <- rowMeans(normcounts[, coldata$condition == condition_levels[1], drop = FALSE])
    res$baseMeanB <- rowMeans(normcounts[, coldata$condition == condition_levels[2], drop = FALSE])
    
    ggplot(res, aes(x = log10(baseMeanA + 1), y = log10(baseMeanB + 1))) + 
      geom_point(alpha = 0.5) + 
      geom_abline(slope = 1, intercept = 0, color = "red") + 
      ggtitle("Expression Scatter Plot") + 
      theme_minimal()
  })
  
  output$volcanoplot <- renderPlot({
    res <- analysis_results()$results
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    title = 'Volcano Plot',
                    pCutoff = 0.05,
                    FCcutoff = 1)
  })
  
  output$heatmap <- renderPlot({
    normcounts <- analysis_results()$normcounts
    top_genes <- order(rowVars(normcounts), decreasing = TRUE)[1:50]
    mat <- normcounts[top_genes, ]
    
    pheatmap(mat, scale = "row", clustering_distance_rows = "correlation",
             show_rownames = FALSE, annotation_col = analysis_results()$coldata)
  })
  
  # === Download Handlers ===
  output$download_ma <- downloadHandler(
    filename = function() { "MA_plot.png" },
    content = function(file) {
      png(file)
      DESeq2::plotMA(analysis_results()$results, main = "MA Plot")
      dev.off()
    }
  )
  
  output$download_box <- downloadHandler(
    filename = function() { "Boxplot.png" },
    content = function(file) {
      png(file)
      boxplot(log2(counts(analysis_results()$ddsDE, normalized = TRUE) + 1),
              col = rainbow(ncol(analysis_results()$ddsDE)), 
              outline = FALSE,
              main = "Boxplot of Normalized Counts")
      dev.off()
    }
  )
  
  output$download_scatter <- downloadHandler(
    filename = function() { "Scatter_Plot.png" },
    content = function(file) {
      res <- as.data.frame(analysis_results()$results)
      normcounts <- analysis_results()$normcounts
      coldata <- analysis_results()$coldata
      condition_levels <- levels(factor(coldata$condition))
      res$baseMeanA <- rowMeans(normcounts[, coldata$condition == condition_levels[1], drop = FALSE])
      res$baseMeanB <- rowMeans(normcounts[, coldata$condition == condition_levels[2], drop = FALSE])
      png(file)
      print(
        ggplot(res, aes(x = log10(baseMeanA + 1), y = log10(baseMeanB + 1))) + 
          geom_point(alpha = 0.5) + 
          geom_abline(slope = 1, intercept = 0, color = "red") + 
          ggtitle("Expression Scatter Plot") + 
          theme_minimal()
      )
      dev.off()
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() { "Volcano_Plot.png" },
    content = function(file) {
      res <- analysis_results()$results
      png(file)
      print(
        EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = 'Volcano Plot',
                        pCutoff = 0.05,
                        FCcutoff = 1)
      )
      dev.off()
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() { "Heatmap.png" },
    content = function(file) {
      normcounts <- analysis_results()$normcounts
      top_genes <- order(rowVars(normcounts), decreasing = TRUE)[1:50]
      mat <- normcounts[top_genes, ]
      png(file)
      pheatmap(mat, scale = "row", clustering_distance_rows = "correlation",
               show_rownames = FALSE, annotation_col = analysis_results()$coldata)
      dev.off()
    }
  )
}

shinyApp(ui, server)
