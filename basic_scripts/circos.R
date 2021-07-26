library(circlize)
library(RColorBrewer)

# display.brewer.all()
# 
# categories <- c("DepMap","TARGET","TCGA")
# categories_cols = brewer.pal(3, "Set2")
# 
# projectss <- unname(unlist(projects[categories]))
# projectss
# circos.par("track.height" = 0.3)
# circos.initialize(categories, xlim = cbind(rep(0, 3), sapply(categories,function(x) length(projects[[x]]))))
# circos.trackHist(graph_df$bigproject, x = graph_df$survival_days,  bin.size = 0.1,
#                  col = "#999999", border = "#999999")
# circos.track(categories, ylim = c(0, 1), track.index = 1, track.height = 0.3, panel.fun = function(x, y) {
#   circos.axis(h = 1, major.tick = F, minor.ticks = F,
#               labels.cex = 0.1, col = categories_cols[CELL_META$sector.numeric.index], labels.col="#ffffff")
#   circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
#               facing = "bending.inside",col = categories_cols[CELL_META$sector.numeric.index])
# }, bg.border = NA)
# 
# circos.clear()

# par(new = TRUE) # <- magic
# circos.par("canvas.xlim" = c(-2, 1), "canvas.ylim" = c(-2, 1))
# circos.initialize(projectss, xlim = c(0,1))
# circos.track(projectss, ylim = c(-1, 1), track.index = 1, track.height = 0.1)
# circos.clear()






# 
# make_box_data <- function(graph_df,df_inner, bigproject){
#   result <- list()
#   for(index in 1:max(df_inner[df_inner$bigproject == bigproject,"subproject_mapping"])){
#     temp = graph_df$survival_days[graph_df$subproject_mapping == index & graph_df$bigproject == bigproject]
#     result[[index]] = temp 
#   }
#   return(result)
# }
# 
# TCGA <- make_box_data(graph_df,df_inner,"TCGA")
# TARGET <- make_box_data(graph_df,df_inner,"TARGET")
# 
# 
# 
# 
# saveRDS(DepMap,file = "DepMap.rds")
# saveRDS(TCGA,file = "TCGA.rds")
# saveRDS(TARGET,file = "TARGET.rds")
# saveRDS(df_inner,file = "df_inner.rds")
# saveRDS(graph_df,file = "graph_df.rds")


#START OF THE CODES FOR GRAPH----
setwd("/Applications/Codes/cSurvival_dev/basic_scripts")

#load RDS

DepMap <- readRDS("DepMap.rds")
TCGA <- readRDS("TCGA.rds")
TARGET <- readRDS("TARGET.rds")
df_inner <- readRDS("df_inner.rds")
graph_df <- readRDS("graph_df.rds")

categories <- c("DepMap","TARGET","TCGA")
categories_cols = brewer.pal(3, "Set2")

#This function calculates the survival rate based on year
calculate_survival_rate <- function(data,year){
  result <- c()
  for(index in 1:length(data)){
    result[[index]] <- length(which(data[[index]] <= (365.25 * year))) / length(data[[index]])
  }
  return(result)
}


#have to manually allocate more space in xlim so that the box plot can be aligned
xlim_matrix <- matrix(nrow = 3,ncol = 2)
xlim_matrix[1,] = c(0.5,3.5)
xlim_matrix[2,] = c(0.5,9.5)
xlim_matrix[3,] = c(0.5,33.5)
row.names(xlim_matrix) <- c("DepMap","TARGET","TCGA")

calculate_ylim <- function(DepMap,TARGET,TCGA,bigproject = c("DepMap","TARGET","TCGA"),survival_days = "survival_days"){
  box_ylim <- matrix(nrow = 3,ncol = 2)
  survival_days <- rep(survival_days,3)
  box_ylim[1,] = c(min(unlist(DepMap)),max(unlist(DepMap)))
  box_ylim[2,] = c(min(unlist(TARGET)),max(unlist(TARGET)))
  box_ylim[3,] = c(min(unlist(TCGA)),max(unlist(TCGA)))
  row.names(box_ylim) <- bigproject
  return(box_ylim)
}


box_ylim
box_ylim <- calculate_ylim(DepMap = DepMap,TARGET = TARGET, TCGA = TCGA)




#START OF GRAPH
circos.par(points.overflow.warning = FALSE)
circos.initialize(df_inner$bigproject, xlim = xlim_matrix)#df_inner$subproject_mapping)
#circos.trackPoints(df_inner$bigproject, x = df_inner$subproject_mapping, y = df_inner$count, pch = 16, cex = 0.5)
circos.track(ylim = c(0,11),bg.border = NA, panel.fun = function(x, y) {
  
})

circos.update(sector.index = "DepMap", track.index = 1,bg.border = NA)
circos.text(x = c(1,1.75:2.75), y = CELL_META$ycenter, labels = df_inner$subproject[df_inner$bigproject == "DepMap"], cex = 0.35, facing = "clockwise", niceFacing = TRUE)
circos.update(sector.index = "TARGET", track.index = 1,bg.border = NA)
circos.text(x = c(1,1.75:8.75), y = CELL_META$ycenter, labels = df_inner$subproject[df_inner$bigproject == "TARGET"], cex = 0.35, facing = "clockwise", niceFacing = TRUE)
circos.update(sector.index = "TCGA", track.index = 1,bg.border = NA)
circos.text(x = c(1,1.75:32.75), y = CELL_META$ycenter, labels = df_inner$subproject[df_inner$bigproject == "TCGA"], cex = 0.35, facing = "clockwise", niceFacing = TRUE)

circos.track(ylim = box_ylim, panel.fun = function(x, y) {
})


# This is for the layer of box plot and survival rate, I wrote a function to calculate survival rate above
circos.update(sector.index = "TCGA", track.index = 2)
circos.text(x = c(1,1.75:32.75), y = CELL_META$cell.ylim[2] * 0.9, labels = paste0(as.character(round(unlist(calculate_survival_rate(TCGA,5)),2) * 100),"%"), cex = 0.4, facing = "bending.outside", niceFacing = TRUE, sector.index	= "TCGA", track.index = 2)
circos.boxplot(TCGA, c(1,1.75:32.75), box_width = 0.5,cex = 0.25)
circos.update(sector.index = "DepMap", track.index = 2)
circos.boxplot(DepMap, c(1,1.75:2.75),box_width = 0.5, cex = 0.25)
circos.update(sector.index = "TARGET", track.index = 2)
circos.text(x = c(0.75,1.75:8.75), y = CELL_META$cell.ylim[2] * 0.9, labels = paste0(as.character(round(unlist(calculate_survival_rate(TARGET,5)),2) * 100),"%"), cex = 0.4, facing = "bending.outside", niceFacing = TRUE, sector.index	= "TARGET", track.index = 2)
circos.boxplot(TARGET, c(1,1.75:8.75),box_width = 0.5, cex = 0.25)

# This is for the layer of histogram and count, the reason I multiply 1.1 here is because I want to create extra room so that the count text can be written above of the histogram
circos.trackHist(ylim = c(min(df_inner$count),max(df_inner$count)*1.1),graph_df$bigproject,x = graph_df$subproject_mapping, col = "#999999", border = "#999999", bin.size = 0.5)
circos.text(x = c(1,1.75:32.75), y = CELL_META$cell.ylim[2] * 0.9, labels = as.character(df_inner$count[df_inner$bigproject == "TCGA"]), cex = 0.4, facing = "bending.outside", niceFacing = TRUE, sector.index	= "TCGA", track.index = 3)
circos.text(x = c(0.75,1.75:8.75), y = CELL_META$cell.ylim[2] * 0.9, labels = as.character(df_inner$count[df_inner$bigproject == "TARGET"]), cex = 0.4, facing = "bending.outside", niceFacing = TRUE, sector.index	= "TARGET", track.index = 3)
circos.text(x = c(1,1.75:2.75), y = CELL_META$cell.ylim[2] * 0.9, labels = as.character(df_inner$count[df_inner$bigproject == "DepMap"]), cex = 0.4, facing = "bending.outside", niceFacing = TRUE, sector.index	= "DepMap", track.index = 3)

#This is for the most inner track
circos.track(df_inner$bigproject,ylim = c(0,1), track.height = 0.2, panel.fun = function(x, y){
  circos.axis(h = 1, major.tick = F, minor.ticks = F,
              labels.cex = 0.1, col = categories_cols[CELL_META$sector.numeric.index], labels.col="#ffffff")
  circos.text(CELL_META$xcenter, y = 0, CELL_META$sector.index, cex = 0.75,
              facing = "clockwise",col = categories_cols[CELL_META$sector.numeric.index], niceFacing = TRUE)
}, bg.border = NA)

circos.clear()



#Calculate sum of tcga and target----
sum( df$columnA )

#sum of TCGA and TARGET, so big project does not equal to DepMap
sum(df_inner[df_inner$bigproject != "DepMap",]$count)

sum(df_inner[df_inner$bigproject == "TARGET",]$count)
sum(df_inner[df_inner$bigproject == "TCGA",]$count)