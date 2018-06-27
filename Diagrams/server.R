#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# install.packages("gplots")
# install.packages("VennDiagram")
#import libraries to use
library(shiny)
library(shinyjs)
library(gplots)
library(VennDiagram)
library(affycoretools)
library(VennDetail)
# Define server logic
shinyServer(function(input, output, session) {
  
  #Show panels
  observeEvent(input$add2, {
    shinyjs::toggle("input2")
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
  
  #Treatment of list for showing according to radioButton selection in venn Diagram panel
  vennIntersectList <- reactive({
    #Treatment of intersections between list of genes
    A <- strsplit(input$genlist1, ",")[[1]]
    B <- strsplit(input$genlist2, ",")[[1]]
    C <- strsplit(input$genlist3, ",")[[1]]
    D <- strsplit(input$genlist4, ",")[[1]]
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
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "(", length(ABCD), ")") ) return (ABCD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist4, "(", length(ABC_NOT_D), ")") ) return (ABC_NOT_D)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist3, "(", length(ABD_NOT_C), ")") ) return (ABD_NOT_C)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist2, "(", length(ACD_NOT_B), ")") ) return (ACD_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "(", length(BCD_NOT_A), ")") ) return (BCD_NOT_A)
      if( as.character(input$select_list) == paste(input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(CD_NOT_AB), ")") ) return (CD_NOT_AB)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist4, "(", length(BC_NOT_AD), ")") ) return (BC_NOT_AD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(BD_NOT_AC), ")") ) return (BD_NOT_AC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist4, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(AD_NOT_BC), ")") ) return (AD_NOT_BC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(AC_NOT_BD), ")") ) return (AC_NOT_BD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(AB_NOT_CD), ")") ) return (AB_NOT_CD)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(A_NOT_BCD), ")") ) return (A_NOT_BCD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(B_NOT_ACD), ")") ) return (B_NOT_ACD)
      if( as.character(input$select_list) == paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(C_NOT_ABD), ")") ) return (C_NOT_ABD)
      if( as.character(input$select_list) == paste(input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(D_NOT_ABC), ")") ) return (D_NOT_ABC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB_NOT_CD), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ABC_NOT_D) , ")") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "(", length(AC_NOT_BD), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(ABC_NOT_D) , ")") ) return (AC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist4, "(", length(AD_NOT_BC), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ACD_NOT_B) , ")") ) return (AD)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "(", length(BC_NOT_AD), "+", length(ABC_NOT_D), "+", length(ABCD), "+", length(BCD_NOT_A) , ")") ) return (BC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist4, "(", length(BD_NOT_AC), "+", length(BCD_NOT_A), "+", length(ABCD), "+", length(ABD_NOT_C) , ")") ) return (BD)
      if( as.character(input$select_list) == paste(input$title_genlist3, "+", input$title_genlist4, "(", length(CD_NOT_AB), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(BCD_NOT_A) , ")") ) return (CD)
    }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "(",length(ABC),")") ) return (ABC)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "(", length(AB_NOT_C),")") ) return (AB_NOT_C)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "(", length(AC_NOT_B),")") ) return (AC_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "(", length(BC_NOT_A),")") ) return (BC_NOT_A)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(A_NOT_BC),")") ) return (A_NOT_BC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(B_NOT_AC),")") ) return (B_NOT_AC)
      if( as.character(input$select_list) == paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(C_NOT_AB),")") ) return (C_NOT_AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(",length(AB_NOT_C), "+", length(ABC), ")") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist3, "(",length(AC_NOT_B), "+", length(ABC), ")") ) return (AC)
      if( as.character(input$select_list) == paste(input$title_genlist2, "+", input$title_genlist3, "(",length(AB_NOT_C), "+", length(ABC), ")") ) return (BC)
    }else if(input$genlist1 != "" & input$genlist2 != ""){
      if( as.character(input$select_list) == paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB),")") ) return (AB)
      if( as.character(input$select_list) == paste(input$title_genlist1, "-", input$title_genlist2, "(", length(A_NOT_B),")") ) return (A_NOT_B)
      if( as.character(input$select_list) == paste(input$title_genlist2, "-", input$title_genlist1, "(", length(B_NOT_A),")") ) return (B_NOT_A)
    }
  })
  #End of Treatment
  #Create matrix from fileinput [HEATMAP]
  dataInput <- reactive({
    if(is.null(input$heatmap_matrix)){
      return (NULL)
    }else{
      path <- input$heatmap_matrix
      samples <- input$samples
      df <- read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
      rownames(df) <- df$input$id_column
      df <- df[, colnames(df) %in% samples]
      df <- data.matrix(df, rownames.force = NA)
      return (df)
    }
  })
  #"Treatment of [HEATMAP] panel textArea"
  finalInput <- reactive({
    if(input$genlist == "") return (dataInput()) # if are not genes informed, then make the first treatment
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
  
  #Observer for contrast button in Heatmap Panel
  observeEvent(input$heat, {
    #Draw heatmap in ui
    output$distPlot <- renderPlot({
      
      thread <- finalInput() #Get matrix from "Treatment of heatmap panel textArea"
      
      if(!is.null(thread)){
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
        
        #Type of density info
        if(as.character(input$density) == "histogram") density <- "histogram"
        if(as.character(input$density) == "density") density <- "density"
        if(as.character(input$density) == "none") density <- "none"
        #Control show/hide dendogram
        if(length(input$dendogram) == 1){
          if(input$dendogram == "Rowv")   show_dendo <- "row"
            
          if(input$dendogram == "Colv")   show_dendo <- "column"
        }
        if(length(input$dendogram) == 2) show_dendo <- "both"
        if(is.null(input$dendogram)) show_dendo <- "none"
        #end Control show/hide dendogram
        #Execute graphs
        return (heatmap.2(main = input$titlematrix,
                          density.info = density,
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
    
    #Create a file with heatmap plot
    output$downloadHeat <- downloadHandler("Heatmaps", function(theFile) {
     
       #set width and height of file to create
       w <- input$width
       h <- input$height
       #Select type of document to create
       if(input$format == "pdf")   pdf(theFile, width = w, height = h)
       if(input$format == "jpeg")  jpeg(theFile, width = w, height = h)
       if(input$format == "png")  png(theFile, width = w, height = h)
       if(input$format == "tiff")  tiff(theFile, width = w, height = h)
       if(input$format == "bmp")  bmp(theFile, width = w, height = h)
       
       thread <- finalInput() #Get matrix from "Treatment of heatmap panel textArea"
       if(!is.null(thread)){
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
         
         #plot
         heatmap.2( main = input$titlematrix,
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
   #Control error message for heatmap Panel
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
     }
     return (NULL)
   })
   #Closing event
  })
  
  #Treatment of input textArea
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
    
    observe({
      #Treatment of intersections between list of genes
      A <- strsplit(input$genlist1, ",")[[1]]
      B <- strsplit(input$genlist2, ",")[[1]]
      C <- strsplit(input$genlist3, ",")[[1]]
      D <- strsplit(input$genlist4, ",")[[1]]
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
                          choices = c(paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "(", length(ABCD), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist4, "(", length(ABC_NOT_D), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist3, "(", length(ABD_NOT_C), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist2, "(", length(ACD_NOT_B), ")"),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "(", length(BCD_NOT_A), ")"),
                                      paste(input$title_genlist3, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(CD_NOT_AB), ")"),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist4, "(", length(BC_NOT_AD), ")"),
                                      paste(input$title_genlist2, "+", input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(BD_NOT_AC), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist4, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(AD_NOT_BC), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(AC_NOT_BD), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(AB_NOT_CD), ")"),
                                      paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(A_NOT_BCD), ")"),
                                      paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "-", input$title_genlist4, "(", length(B_NOT_ACD), ")"),
                                      paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist4, "(", length(C_NOT_ABD), ")"),
                                      paste(input$title_genlist4, "-", input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(D_NOT_ABC), ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB_NOT_CD), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ABC_NOT_D) , ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist3, "(", length(AC_NOT_BD), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(ABC_NOT_D) , ")"),
                                      paste(input$title_genlist1, "+", input$title_genlist4, "(", length(AD_NOT_BC), "+", length(ABD_NOT_C), "+", length(ABCD), "+", length(ACD_NOT_B) , ")"),
                                      paste(input$title_genlist2, "+", input$title_genlist3, "(", length(BC_NOT_AD), "+", length(ABC_NOT_D), "+", length(ABCD), "+", length(BCD_NOT_A) , ")"),
                                      paste(input$title_genlist2, "+", input$title_genlist4, "(", length(BD_NOT_AC), "+", length(BCD_NOT_A), "+", length(ABCD), "+", length(ABD_NOT_C) , ")"),
                                      paste(input$title_genlist3, "+", input$title_genlist4, "(", length(CD_NOT_AB), "+", length(ACD_NOT_B), "+", length(ABCD), "+", length(BCD_NOT_A) , ")"))
                        )
      }else if(input$genlist1 != "" & input$genlist2 != "" & input$genlist3 != ""){
        updateRadioButtons( session, "select_list",
                            choices = c(paste(input$title_genlist1, "+", input$title_genlist2, "+", input$title_genlist3, "(",length(ABC),")"),
                                        paste(input$title_genlist1, "+", input$title_genlist2, "(",length(AB_NOT_C), "+", length(ABC), ")"),
                                        paste(input$title_genlist1, "+", input$title_genlist3, "(",length(AC_NOT_B), "+", length(ABC), ")"),
                                        paste(input$title_genlist2, "+", input$title_genlist3, "(",length(AB_NOT_C), "+", length(ABC), ")"),
                                        paste(input$title_genlist1, "+", input$title_genlist2, "-", input$title_genlist3, "(", length(AB_NOT_C),")"),
                                        paste(input$title_genlist1, "+", input$title_genlist3, "-", input$title_genlist2, "(", length(AC_NOT_B),")"),
                                        paste(input$title_genlist2, "+", input$title_genlist3, "-", input$title_genlist1, "(", length(BC_NOT_A),")"),
                                        paste(input$title_genlist1, "-", input$title_genlist2, "-", input$title_genlist3, "(", length(A_NOT_BC),")"),
                                        paste(input$title_genlist2, "-", input$title_genlist1, "-", input$title_genlist3, "(", length(B_NOT_AC),")"),
                                        paste(input$title_genlist3, "-", input$title_genlist1, "-", input$title_genlist2, "(", length(C_NOT_AB),")")))
      }else if(input$genlist1 != "" & input$genlist2 != ""){
        updateRadioButtons(session, "select_list", choices = c( paste(input$title_genlist1, "+", input$title_genlist2, "(", length(AB),")"),
                                                                paste(input$title_genlist1, "-", input$title_genlist2, "(", length(A_NOT_B),")"),
                                                                paste(input$title_genlist2, "-", input$title_genlist1, "(", length(B_NOT_A),")")) )
      }
    })
  })#Close event
  
  #Download list of selected genes  
  output$download_list <- downloadHandler(
    filename = function(){
      #Get selection to give name of file
      name <- input$select_list
      paste("List_", name, ".txt", sep = "")
    },
    content = function(file) {
      #Get list to write file
      lista <- vennIntersectList()
      writeLines(paste(lista, collapse = "\n"), file)
      # write.table(paste(text,collapse=", "), file,col.names=FALSE)
    }
  )
  
  #Download vennDiagram  
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
                           cat.cex = 1
                            ))
    
    dev.off();
  })#close downloadfile
  
})#Close server
