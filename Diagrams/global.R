list.of.packages <- c("ggplot2", "shinyjs", "shinyBS", "shinyWidgets", "shinycssloaders", "shinydashboard", "shinythemes", "colourpicker", "plotly", "heatmaply", "gplots", "VennDiagram")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(shinythemes)
library(colourpicker)
library(ggplot2)
library(plotly)
library(heatmaply)
library(gplots)
library(VennDiagram)

source("Modules/scatter.R")