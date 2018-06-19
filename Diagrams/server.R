#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

#import libraries to use
library(shiny)
library(shinyjs)
library(gplots)
library(VennDiagram)


# Define server logic
shinyServer(function(input, output, session) {
  
  #Put all the inputs to null again
  observeEvent(input$reset1, {
    output$distPlot <- renderPlot({NULL})
    output$newOne <- renderPlot({NULL})
    output$tabledgenes <- renderDataTable({NULL})
    reset('heatmap_matrix')
    reset('heat_genlist')
  })
  
  
  observeEvent(input$add2, {
    if(input$add2 %% 2 == 1){
      shinyjs::hide(id = "input2")
    }else{
      shinyjs::show(id = "input2")
    }
  })
  observeEvent(input$add3, {
    if(input$add3 %% 2 == 1){
      shinyjs::hide(id = "input3")
    }else{
      shinyjs::show(id = "input3")
    }
  })
  observeEvent(input$add4, {
    if(input$add4 %% 2 == 1){
      shinyjs::hide(id = "input4")
    }else{
      shinyjs::show(id = "input4")
    }
  })
  
  observe({
    x <- 0
    y <- 0
    if(!is.null(input$heatmap_matrix)){
      path <- input$heatmap_matrix
      df <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
      x <- as.vector(colnames(df)[sapply(df, class) == "numeric"])
      updateCheckboxGroupInput(session, "samples",
                         choices = x, selected = x, inline = TRUE)
      y <- as.vector(colnames(df)[sapply(df, class) == "character"])
      updateRadioButtons(session, "id_column",
                         choices = y, inline = TRUE)
      updateCheckboxGroupInput(session, "clustering", choices = c("Rowv", "Colv"), inline = TRUE)
      updateCheckboxGroupInput(session, "dendogram", choices = c("Rowv", "Colv"), inline = TRUE)
    }
    #Updating textarea with selected genes uploaded from file.
    if(!is.null(input$heat_genlist)){
      path <- input$heat_genlist
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist", value = list)
    }
    
    if(!is.null(input$venn_genlist1)){
      path <- input$venn_genlist1
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist1", value = list)
    }
    
    if(!is.null(input$venn_genlist2)){
      path <- input$venn_genlist2
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist2", value = list)
    }
    
    if(!is.null(input$venn_genlist3)){
      path <- input$venn_genlist3
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist3", value = list)
    }
    
    if(!is.null(input$venn_genlist4)){
      path <- input$venn_genlist4
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist4", value = list)
    }
  })
  
  #Execute graphs
    dataInput <- reactive({
      if(is.null(input$heatmap_matrix)){
        return (NULL)
      }else{
        path <- input$heatmap_matrix
        df <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
        rownames(df) <- df$input$id_column
        df <- df[, sapply(df, class) == 'numeric']
        df <- data.matrix(df, rownames.force = NA)
        return (df)
      }
    })
    
    
    finalInput <- reactive({
      if(is.null(input$heat_genlist)) return (dataInput()) # if are not genes informed, then make the first treatment
      # otherwise ... 
      column_name <- as.character(input$id_column) #id column is the radio button to select type of annotation input.
      if(column_name == "ensembl_id"){
        path <- input$heatmap_matrix
        # path2 <- as.vector(input$genlist)
        samples <- input$samples # Selected columns to show. 
        fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
        # genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE, as.is = TRUE)
        selection <- strsplit(input$genlist, ",")[[1]]
        df <- rbind(subset(fullgens, ensembl_id == selection[1]))
        for (row in 2:length(selection)) {
          df <- rbind(df, subset(fullgens, ensembl_id == selection[row]))
        }
        rownames(df) <- df$ensembl_id
        df <- df[, colnames(df) %in% samples]
        df <- data.matrix(df, rownames.force = NA)
      }else if(column_name == "gene_name"){
        path <- input$heatmap_matrix
        # path2 <- input$heat_genlist
        samples <- input$samples
        fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
        # genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE, as.is = TRUE)
        selection <- strsplit(input$genlist, ",")[[1]]  # convert the character vector and extracting elements.
        df <- rbind(subset(fullgens, gene_name == selection[1]))
        for (row in 2:length(selection)) {
          df <- rbind(df, subset(fullgens, gene_name == selection[row]))
        }
        rownames(df) <- df$gene_name
        df <- df[, colnames(df) %in% samples] # Df with Selected columns to show. 
        df <- data.matrix(df, rownames.force = NA) # Convert df to matrix, require to make a heatmap.
        return (df)
      }
    })
    
    #Execute graphs
    observeEvent(input$heat, {
      output$distPlot <- renderPlot({
        thread <- finalInput()
        
        if(is.null(thread)){
          return(NULL)
        }else{
          #Control hide/show row label
          if(input$hide_label_matrix){
            label <- ""
          }else{
            label <- NULL
          }
          #end control hide/show row label
          
          #Control clustering
          boolCol <- NULL
          boolRow <- NULL
          if(length(input$clustering) == 1){ 
            if(input$clustering == "Colv")   boolCol <- TRUE
              
            if(input$clustering == "Rowv")   boolRow <- TRUE
          }
          if(length(input$clustering) == 2){
            boolCol <- TRUE
            boolRow <- TRUE
          }
          if(is.null(input$clustering)){
            boolRow <- FALSE
            boolCol <- NULL
          }
          #end Control clustering
          
          #Control show/hide dendogram
          if(length(input$dendogram) == 1){
            if(input$dendogram == "Rowv")   show_dendo <- "row"
              
            if(input$dendogram == "Colv")   show_dendo <- "column"
          }
          if(length(input$dendogram) == 2) show_dendo <- "both"
          if(is.null(input$dendogram)) show_dendo <- "none"
          #end Control show/hide dendogram
          
          
          return (heatmap.2(main = input$titlematrix,
                    density.info = "density",
                    thread,
                    keysize = 2,
                    Rowv = boolRow,
                    Colv = boolCol,
                    dendrogram = show_dendo,
                    col= input$heat_col,
                    trace = "none",
                    scale = input$scalefull,
                    labRow = label,
                    symm = FALSE))
        }
      })
      
      output$downloadHeat <- downloadHandler("Heatmaps", function(theFile) {
        #open file
        w <- input$width
        h <- input$height
        if(input$format == "pdf")   pdf(theFile, width = w, height = h)
        if(input$format == "jpeg")  jpeg(theFile, width = w, height = h)
        if(input$format == "png")  png(theFile, width = w, height = h)
        if(input$format == "tiff")  tiff(theFile, width = w, height = h)
        if(input$format == "bmp")  bmp(theFile, width = w, height = h)
        
        #plot
        thread <- finalInput()
        if(is.null(thread)){
          return(NULL)
        }else{
          if(input$hide_label_matrix){
            label <- ""
          }else{
            label <- NULL
          }
          #Control clustering
          boolCol <- NULL
          boolRow <- NULL
          if(length(input$clustering) == 1){ 
            if(input$clustering == "Colv")   boolCol <- TRUE
            
            if(input$clustering == "Rowv")   boolRow <- TRUE
          }
          if(length(input$clustering) == 2){
            boolCol <- TRUE
            boolRow <- TRUE
          }
          if(is.null(input$clustering)){
            boolRow <- FALSE
            boolCol <- NULL
          }
          #end Control clustering
          
          #Control show/hide dendogram
          if(length(input$dendogram) == 1){
            if(input$dendogram == "Rowv")   show_dendo <- "row"
            
            if(input$dendogram == "Colv")   show_dendo <- "column"
          }
          if(length(input$dendogram) == 2) show_dendo <- "both"
          if(is.null(input$dendogram)) show_dendo <- "none"
          #end Control show/hide dendogram
          heatmap.2(main = input$titlematrix,
                    density.info = "density",
                    thread,
                    keysize = 2,
                    Rowv = boolRow,
                    Colv = boolCol,
                    dendrogram = show_dendo,
                    col= input$heat_col,
                    trace = "none",
                    scale = input$scalefull,
                    labRow = label,
                    symm = FALSE
          )
        }
        #Close pdf
        dev.off()
      })
      #Control error message
      output$error_content <- renderText({
        if(!is.null(input$heat_genlist)){
          df <- finalInput() #load DATA.FRAME
          path <- input$heat_genlist #load Gene list
          genlist <- strsplit(input$genlist, ",")[[1]]  
          if(any(genlist %in% rownames(df))){
            return (NULL)
          }else{
           return (createAlert(session, "error_message",
                               content = "This annotation file doesn't match with given matrix",
                               style = "warning"))
          }
        }else{
          return (NULL)
        }
      })
      #Closing event
    })
    
    vennInputList <- reactive({
      list <- list()
      if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != "" & input$genlist4 != ""){
        list <- list( A = strsplit(input$genlist1, ",")[[1]],
                      B = strsplit(input$genlist2, ",")[[1]],
                      C = strsplit(input$genlist3, ",")[[1]],
                      D = strsplit(input$genlist4, ",")[[1]])
      }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
        list <- list( A = strsplit(input$genlist1, ",")[[1]],
                      B = strsplit(input$genlist2, ",")[[1]],
                      C = strsplit(input$genlist3, ",")[[1]])
      }else if(input$genlist1 != "" & input$genlist2 != ""){
        list <- list( A = strsplit(input$genlist1, ",")[[1]],
                      B = strsplit(input$genlist2, ",")[[1]])
      }else if(input$genlist1 != ""){
        list <- list( A = strsplit(input$genlist1, ",")[[1]])
      }
      
      return (list)
    })
    
    observeEvent(input$venn, {
      output$vennDiagram <-renderPlot({
        listVenn <- vennInputList()
        listNames <- c(input$title_genlist1, input$title_genlist2, input$title_genlist3, input$title_genlist4)
        listCol <- c(input$col1, input$col2, input$col3, input$col4)
        grid.newpage()
        return (grid.draw(venn.diagram(x = listVenn,
                              filename = NULL,
                              category.names = listNames[1:length(listVenn)],
                              fill = listCol[1:length(listVenn)],
                              cex = 2,
                              cat.cex = 2
                              )))
        
      })
      
    })#Close event
  
    output$downloadVenn <- downloadHandler("VennDiagram", function(theFile) {
      jpeg(theFile);
      listVenn <- vennInputList()
      listNames <- c(input$title_genlist1, input$title_genlist2, input$title_genlist3, input$title_genlist4)
      listCol <- c(input$col1, input$col2, input$col3, input$col4)
      grid.newpage()
      grid.draw(venn.diagram(x = listVenn,
                                     filename = NULL,
                                     category.names = listNames[1:length(listVenn)],
                                     fill = listCol[1:length(listVenn)],
                                     cex = 2,
                                     cat.cex = 2
                                ))
      
      dev.off();
    })#close downloadfile
  
    ##DEBUGGING
    
    # output$hola <- renderText({
    #   # paste(class(as.vector(input$genlist)),input$genlist)
    #   selection <- input$genlist
    #   paste(class(selection), selection)
    # })
    
})#Close server
