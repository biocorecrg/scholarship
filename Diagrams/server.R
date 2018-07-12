
# Define server logic
shinyServer(function(input, output, session) {
  
  #Show panels
  observeEvent(input$add2, {
    shinyjs::toggle("input2")
    # updateCollapse(session, "inputs", open = "input2")
  })
  observeEvent(input$add3, {
    shinyjs::toggle("input3")
  })
  observeEvent(input$add4, {
    shinyjs::toggle("input4")
  })
  #End show panels
  
  #Methods for interactive UI
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
    }
    #Updating textarea with selected genes uploaded from file.
    if(!is.null(input$heat_genlist)){
      path <- input$heat_genlist
      list <- scan(path$datapath, what = "character", sep = "\n")
      updateTextAreaInput(session, "genlist", value = paste(list, collapse = "\n"))
    }
    
    if(!is.null(input$venn_genlist1)){
      path <- input$venn_genlist1
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist1", value = paste(list, collapse = "\n"))
      
    }
    
    if(!is.null(input$venn_genlist2)){
      path <- input$venn_genlist2
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist2", value = paste(list, collapse = "\n"))
    }
    
    if(!is.null(input$venn_genlist3)){
      path <- input$venn_genlist3
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist3", value = paste(list, collapse = "\n"))
    }
    
    if(!is.null(input$venn_genlist4)){
      path <- input$venn_genlist4
      list <- scan(path$datapath, what = "character")
      updateTextAreaInput(session, "genlist4", value = paste(list, collapse = "\n"))
    }
  })
  
  #Create matrix from fileinput [HEATMAP]
  dataInput <- reactive({
    if(is.null(input$heatmap_matrix)){
      return (NULL)
    }else{
      path <- input$heatmap_matrix
      samples <- input$samples
      df <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
      column_name <- as.character(input$id_column) #id column is the radio button to select type of annotation input.
      if (column_name == "ensembl_id") rownames(df) <- df$ensembl_id
      df <- df[, colnames(df) %in% samples]
      df <- data.matrix(df, rownames.force = NA)
      return (df)
    }
  })
  
  #"Treatment of [HEATMAP] panel textArea"
  finalInput <- reactive({
    if(input$genlist == "") return (dataInput()) # if are not genes informed, then make the first treatment
    # otherwise ... 
      path <- input$heatmap_matrix
      # path2 <- as.vector(input$genlist)
      samples <- input$samples # Selected columns to show. 
      fullgens <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
      # genlist <- read.delim(file = path2$datapath, sep = "\t", header = TRUE, as.is = TRUE)
      selection <- strsplit(input$genlist, "\n")[[1]]
      df <- rbind(subset(fullgens, fullgens[,input$id_column] %in% selection[1]))
      for (row in 2:length(selection)) {
        df <- rbind(df, subset(fullgens, fullgens[,input$id_column] %in% selection[row]))
      }
      rownames(df) <- df[,input$id_column]
      df <- df[, colnames(df) %in% samples]
      df <- data.matrix(df, rownames.force = NA)
  })
  
  #Check size of matrix
  observe({
    df <- finalInput()
    if(!is.null(df)){
      if(nrow(df) > 10000){
        shinyjs::disable("heat")
      }else{
        shinyjs::enable("heat")
      }
      if(nrow(df) > 5000){
        shinyjs::disable("clustering")
        shinyjs::disable("dendogram")
      }else{
        shinyjs::enable("clustering")
        shinyjs::enable("dendogram")
      }
      updateSliderInput(session, "minInput", min = round(min(df)-median(df)), max = min(df), value = min(df))
      updateSliderInput(session, "maxInput", min = round(max(df)), max = round(max(df)+median(df)) , value = round(max(df)) )
    }
    
  })
  
  output$time <- renderText({
    df <- finalInput()
    if(!is.null(df)){
      if(nrow(df) > 2000){
        message <- paste("This action can take until 5 minutes, be patient.")
      }
    }
    print (message)
  })
  
  #Treatment of list for showing according to radioButton selection in venn Diagram panel
  vennIntersectList <- function(){
    #Treatment of intersections between list of genes
    A <- strsplit(input$genlist1, "\n")[[1]]
    B <- strsplit(input$genlist2, "\n")[[1]]
    C <- strsplit(input$genlist3, "\n")[[1]]
    D <- strsplit(input$genlist4, "\n")[[1]]
    AB <- intersect(A,B)
    AC <- intersect(A,C)
    AD <- intersect(A,D)
    BC <- intersect(B,C)
    BD <- intersect(B,D)
    CD <- intersect(C,D)
    ABC <- intersect(AB,C)
    ABD <- intersect(AB,D)
    ACD <- intersect(AC,D)
    BCD <- intersect(BC,D)
    ABCD <-intersect(ABC,D)
    
    A_NOT_B <- A[!A %in% B]
    A_NOT_BC <-A_NOT_B[!A_NOT_B %in% C]
    A_NOT_BCD <- A_NOT_BC[!A_NOT_BC %in% D]
    AB_NOT_C <- AB[!AB %in% C]
    AC_NOT_B <- AC[!AC %in% B]
    AB_NOT_CD <- AB_NOT_C[!AB_NOT_C %in% D]
    AC_NOT_BD <- AC_NOT_B[!AC_NOT_B %in% D]
    AD_NOT_B <- AD[!AD %in% B]
    AD_NOT_BC <- AD_NOT_B[!AD_NOT_B %in% C]
    ABD_NOT_C <-ABD[!ABD %in% C]
    ABC_NOT_D <- ABC[!ABC %in% D]
    ACD_NOT_B <- ACD[!ACD %in% B]
    B_NOT_A <- B[!B %in% A]
    B_NOT_AC <- B_NOT_A[!B_NOT_A %in% C]
    B_NOT_ACD <- B_NOT_AC[!B_NOT_AC %in% D]
    BC_NOT_A <- BC[!BC %in% A]
    BD_NOT_A <- BD[!BD %in% A ]
    BD_NOT_AC <- BD_NOT_A[!BD_NOT_A %in% C]
    BC_NOT_AD <- BC_NOT_A[!BC_NOT_A %in% D]
    BCD_NOT_A <- BCD[!BCD %in% A]
    C_NOT_A <- C[!C %in% A]
    C_NOT_AB <- C_NOT_A[!C_NOT_A %in% B]
    C_NOT_ABD <- C_NOT_AB[!C_NOT_AB %in% D]
    CD_NOT_A <- CD[!CD %in% A]
    CD_NOT_AB <- CD_NOT_A[!CD_NOT_A %in% B]
    D_NOT_A <-D[!D %in% A]
    D_NOT_AB <- D_NOT_A[!D_NOT_A %in% B]
    D_NOT_ABC <- D_NOT_AB[!D_NOT_AB %in% C]
    #End of Treatment 
    #Control the radioButton selection for returning an specific list 
    if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != "" & input$genlist4 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "(", length(ABCD), ")", sep = "") ) return (ABCD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist4, "(", length(ABC_NOT_D), ")", sep = "") ) return (ABC_NOT_D)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist3, "(", length(ABD_NOT_C), ")", sep = "") ) return (ABD_NOT_C)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist2, "(", length(ACD_NOT_B), ")", sep = "") ) return (ACD_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "(", length(BCD_NOT_A), ")", sep = "") ) return (BCD_NOT_A)
      if( as.character(input$select_list) == paste(input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(CD_NOT_AB), ")", sep = "") ) return (CD_NOT_AB)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist4, "(", length(BC_NOT_AD), ")", sep = "") ) return (BC_NOT_AD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(BD_NOT_AC), ")", sep = "") ) return (BD_NOT_AC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist4, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(AD_NOT_BC), ")", sep = "") ) return (AD_NOT_BC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(AC_NOT_BD), ")", sep = "") ) return (AC_NOT_BD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(AB_NOT_CD), ")", sep = "") ) return (AB_NOT_CD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(A_NOT_BCD), ")", sep = "") ) return (A_NOT_BCD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(B_NOT_ACD), ")", sep = "") ) return (B_NOT_ACD)
      if( as.character(input$select_list) == paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(C_NOT_ABD), ")", sep = "") ) return (C_NOT_ABD)
      if( as.character(input$select_list) == paste(input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(D_NOT_ABC), ")", sep = "") ) return (D_NOT_ABC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB_NOT_CD), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ABC_NOT_D) , ")", sep = "") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "(", length(AC_NOT_BD), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(ABC_NOT_D) , ")", sep = "") ) return (AC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist4, "(", length(AD_NOT_BC), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ACD_NOT_B) , ")", sep = "") ) return (AD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "(", length(BC_NOT_AD), "+", length(ABC_NOT_D), "+", length(ABCD), "+", length(BCD_NOT_A) , ")", sep = "") ) return (BC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist4, "(", length(BD_NOT_AC), "+", length(BCD_NOT_A), "+", length(ABCD), "+", length(ABD_NOT_C) , ")", sep = "") ) return (BD)
      if( as.character(input$select_list) == paste(input$title_genlist3, "+", input$title_genlist4, "(", length(CD_NOT_AB), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(BCD_NOT_A) , ")", sep = "") ) return (CD)
    }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "(",length(ABC),")", sep = "") ) return (ABC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "(", length(AB_NOT_C),")", sep = "") ) return (AB_NOT_C)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "(", length(AC_NOT_B),")", sep = "") ) return (AC_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "(", length(BC_NOT_A),")", sep = "") ) return (BC_NOT_A)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(A_NOT_BC),")", sep = "") ) return (A_NOT_BC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(B_NOT_AC),")", sep = "") ) return (B_NOT_AC)
      if( as.character(input$select_list) == paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(C_NOT_AB),")", sep = "") ) return (C_NOT_AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(",length(AB_NOT_C), "+", length(ABC), ")", sep = "") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "(",length(AC_NOT_B), "+", length(ABC), ")", sep = "") ) return (AC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "(",length(AB_NOT_C), "+", length(ABC), ")", sep = "") ) return (BC)
    }else if(input$genlist1 != "" & input$genlist2 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB),")", sep = "") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "(", length(A_NOT_B),")", sep = "") ) return (A_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "(", length(B_NOT_A),")", sep = "") ) return (B_NOT_A)
    }else if(input$genlist1 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "(", length(A), ")", sep = "") ) return (A)
    }
  }
  #End of Treatment
  
  uiState <- reactiveValues()
  uiState$readyFlag <- 0
  uiState$readyCheck <- 0
  
  #Observer for contrast button in Heatmap Panel
  observeEvent(input$heat, {
    #Draw heatmap in ui
    output$distPlot <- renderPlotly({
      uiState$readyFlag
      
      isolate({
        if (uiState$readyFlag == uiState$readyCheck) {
          uiState$readyFlag <- uiState$readyFlag+1
          return(NULL)
        }
      })
      
      thread <- finalInput() #Get matrix from "Treatment of heatmap panel textArea"
      #CONTROL HIDE / SHOW ROW LABEL
      if(input$rowlabel){
        label <- NULL
      }else{
        label <- rownames(thread)
      }
      if(!is.null(thread)){
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
        
        color.palette  <- colorRampPalette(c(input$heat_col1, input$heat_col2, input$heat_col3))(256)
        
        #end Control show/hide dendogram
        #Execute graphs
        plot <- heatmaply(thread,
                          main = input$titlematrix,
                          keysize = 2,
                          Rowv = boolRow,
                          Colv = boolCol,
                          dendrogram = show_dendo,
                          limits = c(input$minInput,input$maxInput),
                          col = color.palette,
                          trace = "none",
                          scale = input$scalefull,
                          row_dend_left = F,
                          labRow = label,
                          labCol = colnames(thread),
                          symm = FALSE,
                          plot_method = "plotly",
                          margins = c(50,160,50,100)
                          ) %>% layout( width = input$width, height = input$height)
                 
        plot$elementId <- NULL
        uiState$readyCheck <- uiState$readyFlag
        return (plot)
        }
    })
    
   #Control error message for heatmap Panel
   output$error_content <- renderText({
     if(!is.null(input$heat_genlist)){
       df <- finalInput() #load DATA.FRAME
       genlist <- strsplit(input$genlist, "\n")[[1]]  
       if(any(genlist %in% rownames(df))){
         return (NULL)
       }else{
        return (createAlert(session, "error_message",
                            content = "This annotation file doesn't match with given matrix",
                            style = "warning"))
       }
     }
     return (NULL)
   })
   
   #Closing event
  })
  
  output$downloadHeat <- downloadHandler(
    filename = function(){
      paste("heatmap", input$download_type_heat, sep = ".")
      },
    content = function(file){
      thread <- finalInput() #Get matrix from "Treatment of heatmap panel textArea"
      #CONTROL HIDE / SHOW ROW LABEL
      if(input$rowlabel){
        label <- NULL
      }else{
        label <- rownames(thread)
      }
      if(!is.null(thread)){
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
        
        quantile.range <- quantile(thread, probs = seq(0, 1, 0.01))
        palette.breaks <- seq(quantile.range[input$id], quantile.range[input$id2], 0.01)
        color.palette  <- colorRampPalette(c(input$heat_col1, input$heat_col2, input$heat_col3))
        #end Control show/hide dendogram
        #Execute graphs
        heatmaply(thread,
                  main = input$titlematrix,
                  keysize = 2,
                  Rowv = boolRow,
                  Colv = boolCol,
                  dendrogram = show_dendo,
                  breaks = palette.breaks,
                  colors = color.palette,
                  trace = "none",
                  scale = input$scalefull,
                  row_dend_left = F,
                  height = "1200",
                  labRow = label,
                  labCol = colnames(thread),
                  symm = FALSE,
                  plot_method = "ggplot",
                  margins = c(50,160,50,100),
                  file = file
          ) %>% layout(width = input$width, height = input$height)
      }
  })
  output$hola <- renderText({
    print(paste(input$id22[1],input$id22[2]))
  })
  
  #Treatment of input textArea
  vennInputList <- reactive({
    list <- list()
    if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != "" & input$genlist4 != ""){
      list <- list( A = strsplit(input$genlist1, "\n")[[1]],
                    B = strsplit(input$genlist2, "\n")[[1]],
                    C = strsplit(input$genlist3, "\n")[[1]],
                    D = strsplit(input$genlist4, "\n")[[1]])
    }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
      list <- list( A = strsplit(input$genlist1, "\n")[[1]],
                    B = strsplit(input$genlist2, "\n")[[1]],
                    C = strsplit(input$genlist3, "\n")[[1]])
    }else if(input$genlist1 != "" & input$genlist2 != ""){
      list <- list( A = strsplit(input$genlist1, "\n")[[1]],
                    B = strsplit(input$genlist2, "\n")[[1]])
    }else if(input$genlist1 != ""){
      list <- list( A = strsplit(input$genlist1, "\n")[[1]])
    }
    return (list)
  })
  
  #Observer for contrast in VennDiagram Panel  
  observeEvent(input$venn, {
    #Draw venn Diagram
    output$vennDiagram <-renderPlot({
      #Get list from "Treatment of input textArea"
      listVenn <- vennInputList()
      #Set vector with names of lists
      listNames <- c(input$title_genlist1, input$title_genlist2, input$title_genlist3, input$title_genlist4)
      #Set vector of colors of lists
      listCol <- c(input$col1, input$col2, input$col3, input$col4)
      grid.newpage()
      return (grid.draw(venn.diagram(x = listVenn,
                                     filename = NULL,
                                     category.names = listNames[1:length(listVenn)],
                                     fill = listCol[1:length(listVenn)],
                                     cex = 2,
                                     cat.cex = 1.5)))
      
    })
    
    #VerbatimTextOutput for list of selection according venn Diagram
    output$info <- renderPrint({
      lista <- unlist(vennIntersectList())
      print(lista)
    })
    
    output$title <- renderPrint({
      cat(input$select_list)
    })
    
    observe({
      #Treatment of intersections between list of genes
      A <- strsplit(input$genlist1, "\n")[[1]]
      B <- strsplit(input$genlist2, "\n")[[1]]
      C <- strsplit(input$genlist3, "\n")[[1]]
      D <- strsplit(input$genlist4, "\n")[[1]]
      AB <- intersect(A,B)
      AC <- intersect(A,C)
      AD <- intersect(A,D)
      BC <- intersect(B,C)
      BD <- intersect(B,D)
      CD <- intersect(C,D)
      ABC <- intersect(AB,C)
      ABD <- intersect(AB,D)
      ACD <- intersect(AC,D)
      BCD <- intersect(BC,D)
      ABCD <-intersect(ABC,D)
      
      A_NOT_B <- A[!A %in% B]
      A_NOT_BC <-A_NOT_B[!A_NOT_B %in% C]
      A_NOT_BCD <- A_NOT_BC[!A_NOT_BC %in% D]
      AB_NOT_C <- AB[!AB %in% C]
      AC_NOT_B <- AC[!AC %in% B]
      AB_NOT_CD <- AB_NOT_C[!AB_NOT_C %in% D]
      AC_NOT_BD <- AC_NOT_B[!AC_NOT_B %in% D]
      AD_NOT_B <- AD[!AD %in% B]
      AD_NOT_BC <- AD_NOT_B[!AD_NOT_B %in% C]
      ABD_NOT_C <-ABD[!ABD %in% C]
      ABC_NOT_D <- ABC[!ABC %in% D]
      ACD_NOT_B <- ACD[!ACD %in% B]
      B_NOT_A <- B[!B %in% A]
      B_NOT_AC <- B_NOT_A[!B_NOT_A %in% C]
      B_NOT_ACD <- B_NOT_AC[!B_NOT_AC %in% D]
      BC_NOT_A <- BC[!BC %in% A]
      BD_NOT_A <- BD[!BD %in% A ]
      BD_NOT_AC <- BD_NOT_A[!BD_NOT_A %in% C]
      BC_NOT_AD <- BC_NOT_A[!BC_NOT_A %in% D]
      BCD_NOT_A <- BCD[!BCD %in% A]
      C_NOT_A <- C[!C %in% A]
      C_NOT_AB <- C_NOT_A[!C_NOT_A %in% B]
      C_NOT_ABD <- C_NOT_AB[!C_NOT_AB %in% D]
      CD_NOT_A <- CD[!CD %in% A]
      CD_NOT_AB <- CD_NOT_A[!CD_NOT_A %in% B]
      D_NOT_A <-D[!D %in% A]
      D_NOT_AB <- D_NOT_A[!D_NOT_A %in% B]
      D_NOT_ABC <- D_NOT_AB[!D_NOT_AB %in% C]
      #END of Treatment of intersections between list of genes 
      #Interactive update of radiobuttons depending of number of lists inputted for user.
      if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != "" & input$genlist4 != ""){
        updateRadioButtons(
                          session, "select_list",
                          choices = c(paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "(", length(ABCD), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist4, "(", length(ABC_NOT_D), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist3, "(", length(ABD_NOT_C), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist2, "(", length(ACD_NOT_B), ")", sep = ""),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "(", length(BCD_NOT_A), ")", sep = ""),
                                      paste(input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(CD_NOT_AB), ")", sep = ""),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist4, "(", length(BC_NOT_AD), ")", sep = ""),
                                      paste(input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(BD_NOT_AC), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist4, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(AD_NOT_BC), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(AC_NOT_BD), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(AB_NOT_CD), ")", sep = ""),
                                      paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(A_NOT_BCD), ")", sep = ""),
                                      paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(B_NOT_ACD), ")", sep = ""),
                                      paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(C_NOT_ABD), ")", sep = ""),
                                      paste(input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(D_NOT_ABC), ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB_NOT_CD), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ABC_NOT_D) , ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "(", length(AC_NOT_BD), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(ABC_NOT_D) , ")", sep = ""),
                                      paste(input$title_genlist1, "+", input$title_genlist4, "(", length(AD_NOT_BC), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ACD_NOT_B) , ")", sep = ""),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "(", length(BC_NOT_AD), "+", length(ABC_NOT_D), "+", length(ABCD), "+", length(BCD_NOT_A) , ")", sep = ""),
                                      paste(input$title_genlist2, "+", input$title_genlist4, "(", length(BD_NOT_AC), "+", length(BCD_NOT_A), "+", length(ABCD), "+", length(ABD_NOT_C) , ")", sep = ""),
                                      paste(input$title_genlist3, "+", input$title_genlist4, "(", length(CD_NOT_AB), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(BCD_NOT_A) , ")", sep = ""))
                        )
      }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
        updateRadioButtons( session, "select_list",
                            choices = c(paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "(",length(ABC),")", sep = ""),
                                        paste(input$title_genlist1, "+", input$title_genlist2, "(",length(AB_NOT_C), "+", length(ABC), ")", sep = ""),
                                        paste(input$title_genlist1, "+", input$title_genlist3, "(",length(AC_NOT_B), "+", length(ABC), ")", sep = ""),
                                        paste(input$title_genlist2, "+", input$title_genlist3, "(",length(AB_NOT_C), "+", length(ABC), ")", sep = ""),
                                        paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "(", length(AB_NOT_C),")", sep = ""),
                                        paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "(", length(AC_NOT_B),")", sep = ""),
                                        paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "(", length(BC_NOT_A),")", sep = ""),
                                        paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(A_NOT_BC),")", sep = ""),
                                        paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(B_NOT_AC),")", sep = ""),
                                        paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(C_NOT_AB),")", sep = "")))
      }else if(input$genlist1 != "" & input$genlist2 != ""){
        updateRadioButtons(session, "select_list", choices = c( paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB),")", sep = ""),
                                                                paste(input$title_genlist1, "-", input$title_genlist2, "(", length(A_NOT_B),")", sep = ""),
                                                                paste(input$title_genlist2, "-", input$title_genlist1, "(", length(B_NOT_A),")", sep = "")) )
      }else if(input$genlist1 != ""){
        updateRadioButtons(session, "select_list", choices = c( paste(input$title_genlist1, "(", length(A), ")", sep = "")))
      }
    })
  })#Close event
  
  #Download list of selected genes  
  output$download_list <- downloadHandler(
    filename = function(){
      #Get selection to give name of file
      name <- input$select_list
      paste(name, "_List", ".txt", sep = "")
    },
    content = function(file) {
      #Get list to write file
      lista <- vennIntersectList()
      writeLines(paste(lista, collapse = "\n"), file)
      # write.table(paste(text,collapse=", "), file,col.names=FALSE)
    }
  )
  
  #Download vennDiagram  
  output$downloadVenn <- downloadHandler(
    filename = function(){
      paste("vennDiagram_",Sys.time())
      },
    content = function(theFile) {
      if(input$download_type_venn == "pdf") pdf(theFile, width = input$width_venn, height = input$height_venn)
      if(input$download_type_venn == "png") png(theFile, width = input$width_venn, height = input$height_venn)
      if(input$download_type_venn == "jpeg") jpeg(theFile, width = input$width_venn, height = input$height_venn)
      if(input$download_type_venn == "tiff") tiff(theFile, width = input$width_venn, height = input$height_venn)
      listVenn <- vennInputList()
      listNames <- c(input$title_genlist1, input$title_genlist2, input$title_genlist3, input$title_genlist4)
      listCol <- c(input$col1, input$col2, input$col3, input$col4)
      grid.newpage()
      grid.draw(venn.diagram(x = listVenn,
                             filename = NULL,
                             category.names = listNames[1:length(listVenn)],
                             fill = listCol[1:length(listVenn)],
                             cex = 2,
                             cat.cex = 1
                              ))
  
      dev.off()
  })#close downloadfile
  
  #Call a method from scatter.R file that return a dataframe from the inputted file in df_scatter input.
  df <- callModule(dfFile, "df_scatter", stringAsFactor = FALSE)
  
  #Observer to change options in inputs, depending from file input.
  observe({
    data <- df()
    if(!is.null(data)){
      x <- as.vector(colnames(data)[sapply(data, class) == "numeric"])
      updateSelectInput(session, "y_expression_scatter", choices = x)
      updateSelectInput(session, "x_expression_scatter", choices = x)
    }
    if(!is.null(data)){
      x <- cbind(" ", as.vector(colnames(data)))
      updateSelectInput(session, "colorby", choices = x, selected = NULL)
    }
    if(!is.null(data)){
      x <- cbind(" ", as.vector(colnames(data)[sapply(data, class) == "character"]))
      updateSelectInput(session, "typeby", choices = x, selected = NULL)
    }
    if(!is.null(data)){
      if(nrow(data) > 5000){
        shinyjs::disable("doScatter")
      }else{
        shinyjs::enable("doScatter")
      }
    }
  })
  
  #Get parameters to create scatter plot
  #Return a scatter plot
  drawplot <- reactive({
    data <- df()
    rownames(data) <- data$ensembl_id
    x <- input$x_expression_scatter
    y <- input$y_expression_scatter
    if(!is.null(data)){
      if( input$typeby != " " & input$colorby != " "){
        p <- ggplot( data, aes_string(x = x,
                                      y = y,
                                      color = input$colorby,
                                      shape = input$typeby) ) +
          geom_point( size = input$size) +
          ggtitle(input$scatter_title)+
          xlab(input$scatter_xlab) +
          ylab(input$scatter_ylab)
      }else if(input$colorby != " "){
        p <-ggplot( data, aes_string(x = x,
                                     y = y,
                                     color = input$colorby) ) +
          geom_point( shape = input$type,
                      size = input$size) +
          ggtitle(input$scatter_title)+
          xlab(input$scatter_xlab) +
          ylab(input$scatter_ylab)
      }else if(input$typeby != " "){
        p <- ggplot( data, aes_string(x = x,
                                      y = y,
                                      shape = input$typeby) ) +
          geom_point( color = input$color,
                      size = input$size)+
          ggtitle(input$scatter_title)+
          xlab(input$scatter_xlab) +
          ylab(input$scatter_ylab)
      }else{
        p <- ggplot( data, aes_string(x = x,
                                      y = y)) +
          geom_point( color = input$color,
                      shape = input$type,
                      size = input$size) +
          # xlim(c( min(data$x), max(data$x)*2)) +
          # ylim(c( min(data$y), max(data$y)*2)) +
          ggtitle(input$scatter_title)+
          xlab(input$scatter_xlab) +
          ylab(input$scatter_ylab)
        p$elementId <- NULL
        ggplotly(p)
      }
      
    }
  })
  
  #When diplay button in scatter plot clicked, will show the plot
  observeEvent(input$doScatter, {
    output$scatter <- renderPlotly({
      drawplot()
    })
  })
})#Close server
