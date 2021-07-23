library(shiny)
library(shinydashboard)
library(shinydashboardPlus, include.only=c("boxPad","descriptionBlock")) # "appButton","box","loadingState"
library(shinythemes)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(DT)
library(tidyverse)
library(tidytext)
library(fgsea)
library(survival) # to do the survival analysis
library(survminer) # to plot the survival analysis nicer
library(data.table)
library(grid)
library(RColorBrewer)
library(plotly)
library(htmltools)
library(waiter)
library(shinyalert)
library(shinydisconnect)
library(rintrojs)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(ids)
library(pacman)
library(colorspace)
library(ggsci)
library(plyr)
library(plotly)
library(parallel)
# library(arrangements) # permutations
options(repos = BiocManager::repositories())

source("global/functions.R")
source("global/variables.R")

# --------------- Initialize introjs -------------------
intropath <- paste0(getwd(), "/intro/")
filepaths <- list.files(intropath, full.names=T)
intros <- lapply(filepaths, function(x){
  df <- data.frame(read.csv(x, header=T, sep="\t"))
  df$element <- sub("^", "#", df$element)
  df[df=="#"] <- NA
  df
})
names(intros) <- tools::file_path_sans_ext(list.files(intropath, full.names=F))
rownames(intros) <- NULL
