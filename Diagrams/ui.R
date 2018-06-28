#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
# install.packages("shinyjs")
# install.packages("shinyBS")
# install.packages("shinydashboard")
# install.packages("shinythemes")
# install.packages("colourpicker")
# source("https://bioconductor.org/biocLite.R")
# biocLite("affycoretools")
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinythemes)
library(colourpicker)
library(ggplot2)
library(plotly)
library(heatmaply)
# Define UI for application
shinyUI(
  navbarPage( theme = shinytheme("lumen"), windowTitle = "Biocore",
              img( src = 'logo.png', height = "40" ),
              position = "fixed-top", inverse = TRUE, useShinyjs(), 
              ##Home view
    # tabPanel( title = "Home", icon = icon("home"),
    #   box(width = 2),
    #   box(width = 8, style = "border-bottom: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB;
    #                       border-radius: 1em; height: 88vh; margin-top: 50px;",
    #                 box(width = 2,
    #                   img( src = 'biocore-logo.png', height = "70")
    #                   )
    #                 ),
    #   box(width = 2)
    # ), ##Heatmap view
    tabPanel( title = "Heatmap", icon = icon("fas fa-braille", class = "fa-2x"),
             value = useShinyjs(),
             ## container box
             box( style = "color: #5DADE2;",
                  #Control error
                  box( textOutput("error_content"),bsAlert("error_message"), width = 12),
                  ###Inputs for matrix count
               box(  h1("Configure your plot", style = "color: #707B7C;"),
                 fileInput("heatmap_matrix", "Upload your matrix count", accept = ".txt"),
                 fileInput("heat_genlist", "Upload a list of genes", accept = ".txt"),
                 textAreaInput("genlist", "Paste genes", width = 150, height = 100),
                 textInput("titlematrix", "Title of heatmap", value = "heatmap"),
                 checkboxGroupInput("clustering", "Clustering", choices = c("Rowv", "Colv"),
                                    inline = TRUE),
                 checkboxGroupInput("dendogram", "show/hide dendograms", choices = c("Rowv", "Colv"),
                                    inline = TRUE),
                 selectInput("scalefull", "Scale", choices = c("none", "row", "column"), selected = "none"),
                 checkboxGroupInput("samples", "show/hide samples", choices= c("a","b","c"),
                              inline = TRUE),
                 radioButtons("id_column", "Select the ID column", choices = list("A","B","C"), selected = NULL,
                              inline = TRUE),
                 selectInput("heat_col", "Select colour: ", choices = c("viridis",
                                                                        "heat.colors",
                                                                        "terrain.colors",
                                                                        "BrBG",
                                                                        "cool_warm",
                                                                        "topo.colors",
                                                                        "cm.colors",
                                                                        "bluered",
                                                                        "redblue")),
                 numericInput("width", "choose width: ", value = 1270),
                 numericInput("height", "choose height: ", value = 900),
                 actionButton('heat', 'Contrast', icon("bar-chart-o")),
                 width = 12, style = "border-right: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB; border-bottom: 4px solid #D5DBDB;
                                      border-top-right-radius: 1em; border-bottom-right-radius: 1em;
                                      padding : 0.5em; margin-bottom: 1em; margin-top: 70px;"
               ),
               br(),
               ##Container of buttons
               box( verbatimTextOutput("hola")),
               width = 3, style = "border-radius: 1em;"
               ),
             box( style = "border: 4px solid #D5DBDB;
                          height: 88vh; border-radius: 1em; margin-top: 70px;",
               box( plotlyOutput("distPlot"), width = 11),
             width = 9)
             ),
    ##VennDiagram view
    tabPanel(title = "VennDiagram", icon = icon("adjust", class = "fa-2x")
             ##Container box
             ,box( style = "color: #5DADE2;", width = 6, style = "border-radius: 1em; margin-top: 50px;",
                   #Inputs for selected genes
                   box( h1("Configure your plot", style = "color: #707B7C;"), 
                        box(id = "input1", width = 6,  
                           textInput("title_genlist1", "Name the selection", "A"),
                           fileInput("venn_genlist1", "Upload a list of genes", accept = ".txt"),
                           textAreaInput("genlist1", "Paste genes", width = 150, height = 100),
                           colourInput("col1", "Select colour", "lightblue", allowTransparent = TRUE),
                           actionButton('add2', 'Add', icon("fas fa-plus"))
                          ),
                        box(  id = "input2",
                          width = 6,
                          textInput("title_genlist2", "Name the selection", "B"),
                          fileInput("venn_genlist2", "Upload a list of genes", accept = ".txt"),
                          textAreaInput("genlist2", "Paste genes", width = 150, height = 100),
                          colourInput("col2", "Select colour", "pink", allowTransparent = TRUE),
                          actionButton('add3', 'Add', icon("fas fa-plus"))
                        ),
                        box(  id = "input3",
                          width = 6,
                          textInput("title_genlist3", "Name the selection", "C"),
                          fileInput("venn_genlist3", "Upload a list of genes", accept = ".txt"),
                          textAreaInput("genlist3", "Paste genes", width = 150, height = 100),
                          colourInput("col3", "Select colour", "green", allowTransparent = TRUE),
                          actionButton('add4', 'Add', icon("fas fa-plus"))
                        ),
                        box( id = "input4", 
                          width = 6,
                          textInput("title_genlist4", "Name the selection", "D"),
                          fileInput("venn_genlist4", "Upload a list of genes", accept = ".txt"),
                          textAreaInput("genlist4", "Paste genes", width = 150, height = 100),
                          colourInput("col4", "Select colour", "yellow", allowTransparent = TRUE)
                        ),
                        # tags$div(id = 'input2'),
                        # tags$div(id = 'input3'),
                        # tags$div(id = 'input4'),
                        br(),
                        box(
                          actionButton('venn', 'Contrast', icon("bar-chart-o")),
                          downloadButton("downloadVenn", "Generate plot", icon("fas fa-download")),
                          style = "bottom: 0px;"
                        ),
                        width = 12 #style = "border-right: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB; 
                        #               border-top-right-radius: 1em; border-bottom-right-radius: 1em; 
                        #               padding: 0.5em; margin-bottom: 1em; border-bottom: 4px solid #D5DBDB;
                        #               height: 91vh;"
                      )
                  ),
             box( style = "border: 4px solid #D5DBDB;
                          height: 88vh; border-radius: 1em; margin-top: 70px;",
                  box( plotOutput("vennDiagram"), width = 7),
                  box( 
                    conditionalPanel(
                      condition = "input.genlist1 != '' & input.genlist2 != ''",
                      radioButtons("select_list", " ", choices = c("a", "b"))
                    ),
                  conditionalPanel(
                    condition = "input.genlist1 != '' & input.genlist2 != ''",
                    downloadButton("download_list", "Generate report list", icon("fas fa-download"))
                  ),
                  width = 4
               ),
               box( verbatimTextOutput("info"), width = 12 )
          )
  )
)
)
