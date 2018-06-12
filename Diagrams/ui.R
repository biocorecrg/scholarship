#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinythemes)
# Define UI for application
shinyUI(
  navbarPage( theme = shinytheme("lumen"), windowTitle = "Biocore",
              img( src = 'logo.png', height = "45", style = "position: relative; top: -10px;" ),
              position = "static-top", inverse = TRUE,
    tabPanel( title = "Home", icon = icon("home"),
      box(width = 2),
      box(width = 8, box( width = 12, style = "border-bottom: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB;
                          border-radius: 1em; height: 88vh;", img( src = 'biocore-logo.png', height = "70"))),
      box(width = 2)
    ),
    tabPanel( title = "Heatmap", icon = icon("bar-chart-o"),
             value = useShinyjs(),
             #Inputs for matrix count
             box( style = "color: #5DADE2;",
               box(  h1("Configure your plot", style = "color: #707B7C;"),
                 fileInput("heatmap_matrix", "Upload your matrix count", accept = ".txt"),
                 textInput("titlematrix", "Title of first heatmap"),
                 selectInput("scalefull", "Scale", choices = c("none", "row", "column"), selected = "none"),
                 checkboxInput("hide_label_matrix", "Hide label", value = FALSE), 
                 width = 12, style = "border-right: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB; 
                                      border-top-right-radius: 1em;padding : 0.5em; margin-bottom: 1em;"
               ),#Inputs for selected genes
               br(),
               box(
                 fileInput("heat_genlist", "Upload a list of genes", accept = ".txt"),
                 textInput("title_new_one", "Title of first heatmap"),
                 selectInput("scale_selected_gens", "Scale", choices = c("none", "row", "column"), selected = "none"),
                 checkboxInput("hide_label_selected", "Hide label", value = FALSE),
                 width = 12, style = "border-right: 4px solid #D5DBDB; padding : 0.5em; margin-bottom: 1em;"
               ),#Action buttons
               box( 
                 actionButton('heat', 'Contrast', icon("bar-chart-o")),br(),
                 actionButton('reset1', 'Reset', icon("refresh")),br(),br(),
                 actionButton("tabBut", "Show table", icon("table")),br(),
                 downloadButton("download1", "Generate report"),
                 textOutput("hola"), # this is just for testing.
                 width = 12, style = "border-right: 4px solid #D5DBDB; border-bottom: 4px solid #D5DBDB;
                                      border-bottom-right-radius: 1em; padding : 0.5em;"
               ),
               width = 3, style = "border-radius: 1em;"
               ),
             box( style = "border-top: 4px solid #D5DBDB; border-bottom: 4px solid #D5DBDB; border-left: 4px solid #D5DBDB;
                          height: 88vh; border-radius: 1em;",
               box( plotOutput("distPlot"), status = "primary", width = 11),br(),
               box( plotOutput("newOne"), status = "warning", width = 11, style = "margin-top : 50px"),
               bsModal("modalExample", "Data Table", "tabBut", size = "large", dataTableOutput("tabledgenes")),
             width = 9)
             ),
    tabPanel(title = "VennDiagram", icon = icon("adjust")
             ,box( style = "color: #5DADE2;", width = 3, style = "border-radius: 1em;",
                   #Inputs for selected genes
                   box(  h1("Configure your plot", style = "color: #707B7C;"),
                         fileInput("venn_matrix", "Upload your matrix count", accept = ".txt"),
                         fileInput("venn_genlist1", "Upload a list of genes", accept = ".txt"),
                         fileInput("venn_genlist2", "Upload a list of genes", accept = ".txt"),
                         width = 12, style = "border-right: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB; 
                                      border-top-right-radius: 1em;padding : 0.5em; margin-bottom: 1em;"
                      ),
                   br(),#Action buttons
                   box( 
                     actionButton('venn', 'Contrast', icon("bar-chart-o")),br(),
                     actionButton('reset2', 'Reset', icon("refresh")),br(),br(),
                     actionButton("tabBut2", "Show table", icon("table")),br(),
                     downloadButton("download2", "Generate report"),
                     width = 12, style = "border-right: 4px solid #D5DBDB; border-bottom: 4px solid #D5DBDB;
                                      border-bottom-right-radius: 1em; padding : 0.5em;"
                      )
                  ),
             box(
               box( plotOutput("vennDiagram") )
             )
             )
  )
)

