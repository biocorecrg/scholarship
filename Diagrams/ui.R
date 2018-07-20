# Define UI for application
shinyUI( ui <- function(request){
  navbarPage( 
    theme = shinytheme("lumen"), windowTitle = "Biocore tools",
    img( src = 'logo.png', height = "40" ), 
    position = "static-top", inverse = TRUE, useShinyjs(),
    tags$head(tags$style(type="text/css", "
                                    .busy {
                                       # position: fixed;
                                       # bottom: 0%;
                                        right: 0px;
                                       # width: 20%;
                                       # padding-top: 20px;
                                       # padding-bottom: 20px;
                                       text-align: center;
                                       font-weight: bold;
                                       font-size: 100%;
                                       color: #000000;
                                       background-color: red;
                                       z-index: 105;
                                    }
                                    p{
                                      margin: 0px;
                                    }
                                    h1{
                                      color: #707B7C;
                                      text-align: center;
                                    }
                                     ")
    ),
    header = 
      box(width = 12,
          box(width = 11,
            div(class = "busy",
                p("Please wait, the server is busy..."),
                textOutput("time")
            )
          ),box(width = 1,
            bookmarkButton("Save session")
          )
    ),
    ##Heatmap view
    tabPanel(title = "Heatmap", icon = icon("fas fa-braille", class = "fa-2x"),
             value = useShinyjs(),
             ## Message busy
             
             ## container box
             box( width = 3, style = "border-radius: 1em; color: #5DADE2;",
                  #Control error
               box(width = 12,
                   textOutput("error_content"),
                   bsAlert("error_message")
                  ),
               ###Inputs for matrix count
               box(width = 12, style = "padding : 0.5em; margin-bottom: 1em;",
                   br(),
                   br(),
                   fileInput("heatmap_matrix", "Upload your matrix count", accept = ".txt"),
                   radioButtons("id_column", "Select the ID column", choices = list("A","B","C"), inline = TRUE),
                   fileInput("heat_genlist", "Upload a list of genes", accept = ".txt"),
                   textAreaInput("genlist", "Paste genes"),
                   actionBttn("clear", "Clear selection", icon("fas fa-eraser"),style = "jelly", size = "sm"),
                   dropdownButton(
                     circle = TRUE, status = "secundary", icon = icon("gear"), width = "300px",
                     tooltip = tooltipOptions(title = "Click to see options !"), size = "sm",
                     h1("Configure your plot"),
                     textInput("titlematrix", "Title of heatmap", value = "Heatmap"),
                     checkboxGroupInput("clustering", "Clustering", choices = c("Rowv", "Colv"),
                                        inline = TRUE),
                     checkboxGroupInput("dendogram", "Show/hide dendograms", choices = c("Rowv", "Colv"),
                                        inline = TRUE),
                     selectInput("scalefull", "Scale", choices = c("none", "row", "column"), selected = "none"),
                     checkboxGroupInput("samples", "Show/hide samples", choices= c("A","B","C"),
                                  inline = TRUE),
                     sliderInput("fontsizerow_hm", "Row fontsize: ", min = 8, max = 20, value = 12),
                     sliderInput("fontsizecol_hm", "Col fontsize: ", min = 8, max = 20, value = 12),
                     materialSwitch("rowlabel", "Hide row label", status = "primary"),
                     box(width = 12,
                         h1("Colour options", style = "text-align: left;"),
                         box(width = 3,
                             colourInput("heat_col1", "Select colour", "green", allowTransparent = TRUE)
                         ),
                         box(width = 3,
                             colourInput("heat_col2", "Select colour", "black", allowTransparent = TRUE)
                         ),
                         box(width = 3,
                             colourInput("heat_col3", "Select colour", "red", allowTransparent = TRUE)
                         ),
                         box(width = 12,
                             box(width = 6,
                               sliderInput("minInput", label =  "Range under", min = 0, max = 1, value = 0, dragRange = FALSE) 
                             ),
                             box(width = 6,
                                 sliderInput("maxInput", label =  "Range above", min = 0, max = 1, value = 0, dragRange = FALSE) 
                             )
                         )
                     )
                   ) 
               )
            ),
             box( width = 9, style = "height: 88vh; border-radius: 1em;",
                  box(width = 6,
                      h1("Instructions", style = "text-align: center;"),
                      hr(),
                      box(width = 12,
                          p("1 - The uploaded file should not be of the full genes."),
                          p("2 - Please try to upload file with a maximum of 10000 genes or display button will be disable."),
                          p("3 - You can add a filter file with selected genes or paste directly in the text area."),
                          p("4 - Selection of genes must match with id column, so you should choose type of selection."),
                          p("5 - Input in text area should be separated 1 per line."),
                          p("6 - Clustering and dendogram option buttons will be disable with files above 5000 genes.")
                      )
                  ),
                  box(width = 6,
                    box(width = 6,
                        h1("View options"),
                        hr(),
                        box(width = 6,
                            numericInput("width", "Choose plot width: ", value = 1270)
                        ),
                        box(width = 6,
                            numericInput("height", "Choose plot height: ", value = 750)
                        ),
                        actionBttn('heat', 'Display', icon("bar-chart-o"), style = "jelly", size = "sm")
                    ),
                    box(width = 6,
                        conditionalPanel(
                          condition = "input.heat == true",
                          h1("Download options"),
                          hr(),
                          box(width = 12,
                            radioButtons("download_type_heat", "Select type of file:", choices = c("pdf", "png", "jpeg"), inline = TRUE),
                            downloadBttn("downloadHeat", "Save Heatmap", size = "sm", style = "jelly")
                          )
                        )
                      )
                  ),
                  box(width = 11,
                    plotlyOutput("distPlot")
                  )
             )
    ),
    ##VennDiagram view
    tabPanel(title = "VennDiagram", icon = icon("adjust", class = "fa-2x")
             ##Container box
             ,box( style = "color: #5DADE2;", width = 6, style = "border-radius: 1em;",
                   #Inputs for selected genes
                   box( width = 12, h1("Configure your plot"), 
                        box(id = "input1", width = 6,  
                             textInput("title_genlist1", "Name the selection", "A"),
                             fileInput("venn_genlist1", "Upload a list of genes", accept = ".txt"),
                             textAreaInput("genlist1", "Paste genes"),
                             colourInput("col1", "Select colour", "lightblue", allowTransparent = TRUE),
                             actionBttn('add2', 'Add', icon("fas fa-plus"), style = "material-circle", size = "sm")
                          ),
                        conditionalPanel(
                          condition = "input.add2 > 0",
                          box(id = "input2",
                            width = 6,
                            textInput("title_genlist2", "Name the selection", "B"),
                            fileInput("venn_genlist2", "Upload a list of genes", accept = ".txt"),
                            textAreaInput("genlist2", "Paste genes"),
                            colourInput("col2", "Select colour", "pink", allowTransparent = TRUE),
                            actionBttn('add3', 'Add', icon("fas fa-plus"), style = "material-circle", size = "sm")
                          )
                        ),
                        conditionalPanel(
                          condition = "input.add3 > 0",
                          box(id = "input3",
                            width = 6,
                            textInput("title_genlist3", "Name the selection", "C"),
                            fileInput("venn_genlist3", "Upload a list of genes", accept = ".txt"),
                            textAreaInput("genlist3", "Paste genes"),
                            colourInput("col3", "Select colour", "green", allowTransparent = TRUE),
                            actionBttn('add4', 'Add', icon("fas fa-plus"), style = "material-circle", size = "sm")
                          )
                        ),
                        conditionalPanel(
                          condition = "input.add4 > 0",
                          box(id = "input4", 
                            width = 6,
                            textInput("title_genlist4", "Name the selection", "D"),
                            fileInput("venn_genlist4", "Upload a list of genes", accept = ".txt"),
                            textAreaInput("genlist4", "Paste genes"),
                            colourInput("col4", "Select colour", "yellow", allowTransparent = TRUE)
                          )
                        )
                        # tags$div(id = 'input2'),
                        # tags$div(id = 'input3'),
                        # tags$div(id = 'input4'),
                         #style = "border-right: 4px solid #D5DBDB; border-top: 4px solid #D5DBDB; 
                        #               border-top-right-radius: 1em; border-bottom-right-radius: 1em; 
                        #               padding: 0.5em; margin-bottom: 1em; border-bottom: 4px solid #D5DBDB;
                        #               height: 91vh;"
                      )
                  ),
             box( width = 6, style = "height: 88vh; border-radius: 1em;",
                conditionalPanel(
                  condition = "input.venn > 0",
                  box(width = 12,
                    h1("Download options"),
                    hr(),
                    box(width = 3,
                      numericInput("width_venn", "Choose image width: ", value = 600)
                    ),
                    box(width = 3,
                      numericInput("height_venn", "Choose image height: ", value = 600)
                    ),
                    box(width = 3,
                      selectInput("download_type_venn", "Select type of file:", choices = c("tiff","pdf", "png", "jpeg"))
                    ),
                    box(width = 3,
                      downloadBttn("downloadVenn", "Save venn", size = "sm", style = "jelly")
                    ),
                    box(width = 7, plotOutput("vennDiagram")),
                    box(width = 5,
                        radioButtons("select_list", " ", choices = c("a", "b"))
                    )
                  )
               ),
               box(width = 12,
                   conditionalPanel(
                     condition = "(input.genlist1 != '' | input.genlist2 != '' | input.genlist3 != '' | input.genlist4 != '')",
                     box(width = 6,
                       actionBttn('venn', 'Display', icon("bar-chart-o"), style = "jelly")
                    )
                   ),
                   box(width = 3,
                     verbatimTextOutput("title")
                   ),
                   box(width = 2,
                       conditionalPanel(
                         condition = "input.venn > 0",
                         downloadBttn("download_list", "Save list", size = "sm", style = "jelly")
                       )
                   ),
                   box(width = 12,
                       verbatimTextOutput("info")
                   )
                )
          )
    ),
    tabPanel(
      title = "Scatter plot", icon = icon("fas fa-industry", class = "fa-2x"),
      box(width = 3, style = "border-radius: 1em; color: #5DADE2;",
          dfFileInput("df_scatter", "Upload gene expresion data frame"),
          selectInput("y_expression_scatter", "Select sample for Y axis:", choices = c("a","b","c")),
          sliderInput("ylim_scatter", "Y limits:", min = 0, max = 1, value = c(0,1), step = 1),
          selectInput("x_expression_scatter", "Select sample for X axis:", choices = c("a","b","c")),
          sliderInput("xlim_scatter", "X limits:", min = 0, max = 1, value = c(0,1), step = 1),
          dropdownButton(
            circle = TRUE, status = "secundary", icon = icon("gear"), width = "300px",
            tooltip = tooltipOptions(title = "Click to see options !"), size = "sm",
            textInput("scatter_title", "Name this plot", "Graph"),
            sliderInput("scatter_title_size", "Title size: ", 8, 20, value = 12),
            textInput("scatter_ylab", "Name of Y axis", "Sample Y"),
            sliderInput("scatter_ylab_size", "Y label size: ", 8, 20, value = 12),
            textInput("scatter_xlab", "Name of X axis:", "Sample X"),
            sliderInput("scatter_xlab_size", "X label size: ", 8, 20, value = 12),
            colourInput("color_scatter", "Select points colour", "#48AB6F", allowTransparent = TRUE),
            selectInput("type_scatter", "Point type: ", choices = c("1", "2", "3", "4", "5", "9", "10", "13", "22")),
            sliderInput("size_scatter", "Point size: ", 1, 10, value = 2),
            selectInput("colorby_scatter", "Colour by: ", choices = c("a","b","c")),
            selectInput("typeby_scatter", "Type by: ", choices = c("a","b","c"))
          )
      ),
      box(width = 9,
          conditionalPanel(
            condition = "input.doScatter == false",
            box( actionBttn("doScatter", "Display", icon("bar-chart-o"), style = "jelly", size = "sm"))
          ),
          box (width = 12, plotlyOutput("scatter",height = "80vh"))
      )
    ),
    tags$script( HTML(
      " setInterval(function(){
        if( ($('html').attr('class')=='shiny-busy') ){
                setTimeout(function() {
                  if ($('html').attr('class')=='shiny-busy') {
                      if($('#textit').html()!='Waiting...' ){
                          $('div.busy').show()
                      }
                  }   
                },100) 
              } else {
                $('div.busy').hide()
              }
            },3000) "
    ))
  )
})
