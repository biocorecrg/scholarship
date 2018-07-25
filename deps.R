list.of.packages <- c("shiny", "ggplot2", "processx", "shinyjs", "DT", "shinyBS", "gridExtra", "shinyWidgets", "tidyverse", "cowplot", "shinycssloaders", "shinydashboard", "shinythemes", "colourpicker", "plotly", "heatmaply", "gplots", "VennDiagram")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

