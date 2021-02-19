library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(DT)
library(tidyverse)
library(tidytext)
library(fgsea)
library(data.table)
library(RColorBrewer)
library(plotly)
library(htmltools)
library(waiter)
library(shinyalert)
library(shinydisconnect)
library(rintrojs)

options(repos = BiocManager::repositories())

bcols = c("#FFFDE7","#FFEBEE") # colors for parameter wellpanels

source("global/functions.R")
