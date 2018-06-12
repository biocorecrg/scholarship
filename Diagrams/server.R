#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(gplots)
library(VennDiagram)
# Define server logic
shinyServer(function(input, output) {
   
  #Put all the inputs to null again
  observeEvent(input$reset1, {
    output$distPlot <- renderPlot({NULL})
    output$newOne <- renderPlot({NULL})
    output$tabledgenes <- renderDataTable({NULL})
    reset('matrix')
    reset('genlist')
  })
  
  #Execute graphs
  observeEvent(input$heat, {
    ## Matrix counts heatmap
    output$distPlot <- renderPlot({
      if(is.null(input$heatmap_matrix)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        thread <- read.delim(file = path$datapath, sep = "\t", header = TRUE)
        rownames(thread) <- thread$Gene
        thread <- thread[-1]
        thread <- data.matrix(thread, rownames.force = NA)
        if(is.null(thread)){
          return(NULL)
        }else{
          if(input$hide_label_matrix){
            label <- ""
          }else{
            label <- NULL
          }
          heatmap.2(main = input$titlematrix,
                    density.info = "density",
                    thread,
                    col=bluered(75),
                    trace = "none",
                    scale = input$scalefull,
                    labRow = label
          )
        }
      }
    }) #Closing matrix counts heatmap
    
    #Table of selected genes
    output$tabledgenes <- renderDataTable({
      if(is.null(input$heatmap_matrix)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        if(is.null(input$heat_genlist)){
          return (NULL)
        }else{
          path2 <- input$heat_genlist
          genlist <- as.vector(input$tocheck)
          fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE)
          genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE)
          selection <- as.vector(t(genlist))
          final <- rbind(subset(fullgens, Gene == selection[1]))
          for (row in 2:length(selection)) {
            final <- rbind(final, subset(fullgens, Gene == selection[row]))
          }
          final
        }
      }
    })#closing selected genes table
    ## Heatmap with selected genes
    output$newOne <- renderPlot({
      if(is.null(input$heatmap_matrix) | is.null(input$heat_genlist)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        path2 <- input$heat_genlist
        fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE)
        genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE)
        selection <- as.vector(t(genlist))
        thread2 <- rbind(subset(fullgens, Gene == selection[1]))
        for (row in 2:length(selection)) {
          thread2 <- rbind(thread2, subset(fullgens, Gene == selection[row]))
        }
        rownames(thread2) <- thread2$Gene
        thread2 <- thread2[-1]
        thread2 <- data.matrix(thread2, rownames.force = NA)
        if(is.matrix(thread2)){
          if(input$hide_label_selected){
            label <- ""
          }else{
            label <- NULL
          }
          heatmap.2(main = input$title_new_one,
                    density.info = "density",
                    thread2,
                    col=bluered(75),
                    trace = "none",
                    labRow = label,
                    scale = input$scale_selected_gens)
        }
      }
    })#Closing genes heatmap
    
    output$download1 <- downloadHandler("Heatmaps", function(theFile) {
      #open pdf
      pdf(theFile)
      
      #plot 1
      if(is.null(input$heatmap_matrix)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        thread <- read.delim(file = path$datapath, sep = "\t", header = TRUE)
        rownames(thread) <- thread$Gene
        thread <- thread[-1]
        thread <- data.matrix(thread, rownames.force = NA)
        if(is.null(thread)){
          return(NULL)
        }else{
          if(input$hide_label_matrix){
            label <- ""
          }else{
            label <- NULL
          }
          heatmap.2(main = input$titlematrix,
                    density.info = "density",
                    thread,
                    col=bluered(75),
                    trace = "none",
                    scale = input$scalefull,
                    labRow = label
          )
        }
      }
      #plot 2
      if(is.null(input$heatmap_matrix) | is.null(input$heat_genlist)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        path2 <- input$heat_genlist
        fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE)
        genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE)
        selection <- as.vector(t(genlist))
        thread2 <- rbind(subset(fullgens, Gene == selection[1]))
        for (row in 2:length(selection)) {
          thread2 <- rbind(thread2, subset(fullgens, Gene == selection[row]))
        }
        rownames(thread2) <- thread2$Gene
        thread2 <- thread2[-1]
        thread2 <- data.matrix(thread2, rownames.force = NA)
        if(is.matrix(thread2)){
          if(input$hide_label_selected){
            label <- ""
          }else{
            label <- NULL
          }
          heatmap.2(main = input$title_new_one,
                    density.info = "density",
                    thread2,
                    col=bluered(75),
                    trace = "none",
                    labRow = label,
                    scale = input$scale_selected_gens)
        }
      }
      
      #Close pdf
      dev.off()

      # file.copy(from = "heatmaps.pdf", to = theFile)
    })  
  })
  #Closing event
  
  observeEvent(input$venn, {
    output$vennDiagram <-renderPlot({
      if(is.null(input$venn_matrix) | is.null(input$venn_genlist1) | is.null(input$venn_genlist2)){
        return (NULL)
      }else{
        path_matrix <- input$venn_matrix
        path_genlist1 <-input$venn_genlist1
        path_genlist2 <- input$venn_genlist2
        fullgens <- read.table(path_matrix$datapath, header = TRUE, sep = "\t", as.is = TRUE)
        genlist1 <- read.table(path_genlist1$datapath, header = TRUE, sep = "\t", as.is = TRUE)
        genlist2 <- read.table(path_genlist2$datapath, header = TRUE, sep = "\t", as.is = TRUE)
        
        selection1 <- as.vector(t(genlist1))
        selection2 <- as.vector(t(genlist2))
        
        final1 <- rbind(subset(fullgens, Gene == selection1[1]))
        for (row in 2:length(selection1)) {
          final1 <- rbind(final1, subset(fullgens, Gene == selection1[row]))
        }
        
        final2 <- rbind(subset(fullgens, Gene == selection2[1]))
        for (row in 2:length(selection2)) {
          final2 <- rbind(final2, subset(fullgens, Gene == selection2[row]))
        }
        
        overlapp <- subset(final1, Gene %in% final2$Gene)
        
        grid.newpage()
        draw.pairwise.venn(area1 = nrow(final1),
                           area2 = nrow(final2),
                           cross.area = nrow(overlapp),
                           category = c("Selection 1", "Selection 2"),
                           fill = c("lightblue", "pink"),
                           scaled = TRUE)
      }
    })
    
  })#Close event
  
})#Close server
