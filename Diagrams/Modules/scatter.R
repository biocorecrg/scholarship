

dfFileInput <- function(id, label = "Upload file") {
  
  ns <- NS(id)
  
  tagList(
    fileInput(ns("file_df"), "Upload gene expression data frame")  )
}

dfFile <- function(input, output, session, stringAsFactors){
  
  userFile <- reactive({
    # If no file is selected, don't do anything
    validate(need(input$file_df, message = FALSE))
    input$file_df
  })
  
  dataframe <- reactive({
    ##read.delim(file = path$datapath, sep = "\t", header = TRUE, as.is = TRUE)
    read.delim(userFile()$datapath,
               sep = "\t",
               header = TRUE,
               as.is = TRUE,
               stringsAsFactors = stringsAsFactors)
  })
  
  return (dataframe)
}