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
                                    #help{
                                      position: auto;
                                      bottom: 0px;
                                    }
                                     ")
    ),
    header = 
      box(width = 12,
          box(width = 11,
            div(id = "help",
                p("If you want to try the app and you don't have data files, you can download it from "),
                a("here", href="https://github.com/ogonza46/scholarship/tree/master/Diagrams/www/files"), style = "float:left;"
            ),
            div(class = "busy",
                p("Please wait, the server is busy..."),
                textOutput("time")
            )
          ),box(width = 1,
            bookmarkButton("Save session", status = "success")
          )
    ),
    tabPanel( title = "Help", icon = icon("fas fa-info-circle", class = "fa-2x"),
              box(width = 4,
                  h1("Heatmap instructions"),
                  dropdownButton(
                    circle = TRUE, status = "secundary", icon = icon("fas fa-braille"), width = "31vw",
                    tooltip = tooltipOptions(title = "Click to see instructions !"), size = "sm",
                    box(width = 12,
                        p("· This tool allows you to generate and customize heatmaps."),br(),
                        p("· The data requested for generate heatmaps is a data frame of gene expresion, with at least one annotation's column."),br(),
                        p("· The data frame should not be longer than 10000 rows, otherwise, you can select an annotation's column to make a selection of genes to display."),br(),
                        p("· You can upload a filter file with selection of genes or paste directly in the text area."),br(),
                        p("· In order to display the heatmap, button will be disabled at least the data frame or selection it is 2000 or less rows."),br(),
                        p("· Once you display, you can change the data frame or the selection to show different data, but TAKE CARE. If you clear selection, and you already changed data frame and this one is larger than 10000 row, heatmap won't be available. "),br(),
                        icon("gear"),
                        p("· the clustering and dendogram buttons will be enable/disable to click depending on number of rows of your data, because a long data frame would make a long process to clustering"),br(),
                        p("· Variance input initialize to 0, but it's acording to the data. That way you can make selection of data to show by variance."),br(),
                        p("· For a better costumization of your heatmap, you can select colors, make scale by row or column, hide row label and change the size of the labels."),br(),
                        p("· Once you click display, just wait. Clicking again will hide download options."),br(),
                        p("· You can select the format file, or directly download png of the graph in the plot widget."),br(),
                        p("· For every tool you have the chance to save the state and restore it after, using the provided link from \"Save Session\" button.", style = "color: green;")
                    )
                  )
              ),
              box(width = 4,
                  h1("Venn Diagram instructions"),
                  dropdownButton(
                    circle = TRUE, status = "secundary", icon = icon("adjust"), width = "31vw",
                    tooltip = tooltipOptions(title = "Click to see instructions !"), size = "sm",
                    box(width = 12,
                        p("· This tool allows you to find differences between even four lists and generate a venn Diagram."),br(),
                        p("· You can upload a list file, or paste list in the text area, one per row."),br(),
                        p("· To add new list, you can click the + button."),br(),icon("fas fa-plus"),
                        p("· At least, one list should be introduced to show display button."),br(),
                        p("· Once you click display button, you have a graph and with options to download, select section, and show a list with data selected."),br(),
                        p("· You can change the name of lists and color of sections."),br(),
                        p("· You can select the format file to download the image, and size of this one."),br(),
                        p("· For every tool you have the chance to save the state and restore it after, using the provided link from \"Save Session\" button.", style = "color: green;")
                    )
                  )
              ),
              box(width = 4,
                  h1("Scatter plot instructions"),
                  dropdownButton(
                    circle = TRUE, status = "secundary", icon = icon("fas fa-industry"), width = "31vw",
                    tooltip = tooltipOptions(title = "Click to see instructions !"), size = "sm",
                    box(width = 12,
                        p("· This tool requires a dataframe of gene expression to calculate data and generate a scatter plot."),br(),
                        p("· The uploaded dataframe will be threated and you will be able to select the column to represent in x and y axis, by numerics columns inside the dataframe."),br(),
                        icon("gear"),p("· You can, set up the limits, change name and size of axis labels or title, choose colour or shape or size of the points, and finally group the expressed points for colour or shape."),br(),
                        p("· The button \"Display\" will be enable to click unless the dataframe it's bigger than 5000 genes, TAKE CARE."),br(),
                        p("· Once you click the button, if colour by or type by selects were filled, you may wait until plot calculate and it generate."), br(),
                        p("· You will be able to show the correlation, p-value inside the plot. Show a table with (y = a + bx) function data, and download it in \"pdf\", \"png\", \"jpeg\" format."),br(),
                        p("· For every tool you have the chance to save the state and restore it after, using the provided link from \"Save Session\" button.", style = "color: green;")
                    )
                  )
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
                     sliderInput("variance", "Variance", min = 0, max = 1, value = 0),
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
                  box(width = 9,
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
                            radioButtons("download_type_heat", "Select type of file:", choices = c("png", "pdf", "jpeg"), inline = TRUE),
                            downloadBttn("downloadHeat", "Save Heatmap", size = "sm", style = "jelly")
                          )
                        )
                      )
                  )
                  ,box(width = 11,
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
          conditionalPanel(
            condition = "input.doScatter == true",
            box (width = 12, 
                 box(width = 2,
                     materialSwitch("cor_scatter", "Show correlation", status = "primary")
                 ),
                 box(width = 2,
                     materialSwitch("p_value_scatter", "Show p-Value", status = "primary")
                 ),
                 box(width = 2,
                     materialSwitch("stats_table", "Show table", status = "primary")
                 ),
                 box(width = 4,
                       h1("Download options"),
                       hr(),
                       box(width = 6,
                           selectInput("download_type_scatter", "Select type of file:", choices = c("pdf", "png", "jpeg"))
                       ),
                       box(width = 6,
                           downloadBttn("downloadScatter", "Save plot", size = "sm", style = "jelly")
                       )
                 ),
                 dataTableOutput("summaryStats", width = 600),
                 plotlyOutput("scatter",height = "80vh")
                 )
          )
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
