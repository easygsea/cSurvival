install.packages("plotly")
install.packages("tidyverse")
install.packages("htmlwidgets")
install.packages("EnvStats")
library(circlize)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(colorspace)
library(EnvStats)
library(RColorBrewer)

#SET UP----
#read data.rds
data <- readRDS("~/Downloads/rstudio-export/rna_data.rds")
df_o <- readRDS("~/Downloads/rstudio-export/df_o.rds")
rvdepmap <- readRDS("~/Downloads/rstudio-export/rv_depmap.rds")



#Helper function
# original survival df


# generate survival df
generate_surv_df <- function(df, patient_ids, exp, q){
  # generate the data from for model first
  gene_quantiles <- exp %>%
    sapply(function(x) ifelse(x > q, "High", "Low"))
  names(gene_quantiles) <- patient_ids
  
  # # generate survival analysis df
  # df$level <- gene_quantiles[match(df$patient_id,names(gene_quantiles))]
  df$level <- gene_quantiles
  lels <- unique(df$level) %>% sort(.,decreasing = T)
  df$level <- factor(df$level, levels = lels)
  
  df <- df %>% dplyr::filter(!is.na(patient_id))
  return(df)
}



surv_cox <- function(df, mode=1){
  if(rvdepmap){
    if(mode == 1){
      coxph(Surv(dependency) ~ level, data = df)
    }else if(mode == 2){
      coxph(Surv(dependency) ~ level.x * level.y, data = df)
    }
  }else{
    if(mode == 1){
      coxph(Surv(survival_days, censoring_status) ~ level, data = df)
    }else if(mode == 2){
      coxph(Surv(survival_days, censoring_status) ~ level.x * level.y, data = df)
    }
  }
}
# surv_km <- function(df, mode=1){
#   if(rv$depmap){
#     if(mode == 1){
#       survfit(Surv(dependency) ~ level, data = df)
#     }else if(mode == 2){
#       survfit(Surv(dependency) ~ level.x * level.y, data = df)
#     }
#   }else{
#     if(mode == 1){
#       survfit(Surv(survival_days, censoring_status) ~ level, data = df)
#     }else if(mode == 2){
#       survfit(Surv(survival_days, censoring_status) ~ level.x * level.y, data = df)
#     }
#   }
# }


#Function
get_info_most_significant_rna <- function(data, min, max, step, mode="g", df_o){
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)
  
  #initialize list of p value
  
  #p_df = c()
   p_df <- data.frame(integer(),
                        character(),
                    stringsAsFactors=FALSE)
   
  
  # initialize the most significant p value and df
  least_p_value <- 1; df_most_significant <- NULL
  
  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  if(mode == "g"){
    exp <-data[,2] %>% unlist(.) %>% unname(.)
  }else if(mode == "gs"){
    exp <- rowMeans(data[,-1],na.rm=T) %>% unlist(.) %>% unname(.)
  }
  
  # retrieve survival analysis df_o
  df_o <- df_o
  
  # the quantiles we will use to define the level of gene percentages
  quantiles <- quantile(exp, quantile_s, na.rm = T)
  
  for(i in seq_along(quantiles)){
    q <- quantiles[i]
    df <- generate_surv_df(df_o, patient_ids, exp, q)
    hr <- coef(summary(surv_cox(df)))[,2]
    
    # # test if there is significant difference between high and low level genes
    # if(rv$cox_km == "cox"){
    # surv_diff <- surv_cox(df)
    # p_diff <- coef(summary(surv_diff))[,5]
    # }else if(rv$cox_km == "km"){
    
    #if(rv$depmap){
    if(rvdepmap){
      surv_diff <- survdiff(Surv(dependency) ~ level, data = df)
    }else{
      surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
    }
    p_diff <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    
    #append current p value to the p value df
    new_row = c(p_diff,unlist(strsplit(names(quantiles[i]),split = '%',fixed=T)),quantiles[i],hr)
    p_df <- rbind(p_df,new_row)
    # }
    if(!is.na(p_diff)){
      if(p_diff <= least_p_value){
        least_p_value <- p_diff
        df_most_significant <- df
        least_hr <- hr
        cutoff_most_significant <- names(quantiles[i])
        
      }
    }
  }
  #Transform the p_df a little bit to make it work with the ggplot
  colnames(p_df) = c('p_value','quantile','expression','hr')
  p_df$p_value = as.numeric(p_df$p_value)
  p_df$quantile = as.numeric(p_df$quantile)
  p_df$expression = as.numeric(p_df$expression)
  p_df$hr = as.numeric(p_df$hr)
  # proceed only if enough data
  if(is.null(df_most_significant)){
    return(NULL)
  }else{
    results <- list(
      df = df_most_significant,
      cutoff = cutoff_most_significant,
      p_df = p_df
      ,hr = least_hr
    )
    return(results)
  }
}

res <- get_info_most_significant_rna(data,min = 0.2,max = 0.8, step = 0.01, df_o = df_o)

quantile_df$hr
# 
# graph <- ggplot(data=quantile_df, aes(x = quantile_df$quantile, y = quantile_df$p_value, group=1)) +
#   geom_line(linetype = "dashed")+
#   geom_point()+labs(x = "Quantile(%)", y = "P Value") +
#   scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75,1))+
#   scale_x_continuous(limits=c(20, 80), breaks=c(20,40,60,80))
# 
# graph <- ggplotly(graph) %>% 
#   layout(hovermode = "x unified")
# graph
#   
#custom_colorscale = brewer.pal(n = length(quantile_df$p_value), "YlOrRd")[1,length(quantile_df$p_value)]



col_scale_no <- c(0, 0.16666666, 0.33333333333, 0.5, 0.66666666, 1)
col_scale <- sequential_hcl(5, palette = "YlOrRd") %>% rev(.) %>% col2rgb()
col_scale <- sapply(1:ncol(col_scale), function(x) {x <- col_scale[,x]; paste0("rgb(",paste0(x, collapse = ", "),")")})
col_scale <- c("rgb(255, 255, 255)", col_scale)
col_scale <- lapply(1:length(col_scale), function(i) list(col_scale_no[i], col_scale[i]))


#Blue and Red Orange Yellow colorscale for Hazard Ratio graph
col_scale_hr_no <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
col_scale_temp <- sequential_hcl(5, palette = "blues3", rev = TRUE) %>% rev(.) %>%  col2rgb()
col_scale_temp<-col_scale_temp[,-ncol(col_scale_temp)] # delete last column

col_scale_hr <- sequential_hcl(1, palette = "YlOrRd") %>% rev(.) %>% col2rgb()
col_scale_hr <- cbind(col_scale_temp, col_scale_hr)
rm(col_scale_temp)
col_scale_hr <- sapply(1:ncol(col_scale_hr), function(x) {x <- col_scale_hr[,x]; paste0("rgb(",paste0(x, collapse = ", "),")")})
col_scale_hr <- c("rgb(255, 255, 255)", col_scale_hr)
col_scale_hr <- lapply(1:length(col_scale_hr), function(i) list(col_scale_hr_no[i], col_scale_hr[i]))
col_scale_hr



#FUNCTIONS ARE HERE----
single_plot <- function(quantile_df, index){
  #Due to bug in Plotly, I am using the solution mentioned in this website:
  #https://stackoverflow.com/questions/55251470/missing-data-when-supplying-a-dual-axis-multiple-traces-to-subplot
  #If index == 1, this is first graph and we need to pass in normal y axis
  if(index == 1){
    yaxes = c("y","y2")
    HR_Axis <- list(overlaying = yaxes[1],
                            side = "right", title = "Hazard Ratio")
    P_Axis <- list(side = "left", title = "P Value")
    }
  else{
    yaxes = c("y2","y3")
    HR_Axis <- list(overlaying = yaxes[2],
                    side = "right", title = "Hazard Ratio")
    P_Axis <- list(side = "left", title = "P Value")
  }
  
  fig <- plot_ly(quantile_df, x = quantile_df$quantile)
  fig <- fig %>% add_trace(y = ~quantile_df$p_value,type = 'scatter',#color =~p_value,
                           line = list(color = 'rgb(173,173,173)', width = 2),
                           yaxis = yaxes[1],
                           name = 'P Value',
                           marker=list(
                             color=~p_value,
                             colorscale=col_scale,
                             cmid = 0.5,
                             reversescale =TRUE
                           ),
                           text = quantile_df$expression,
                           name = '',mode = 'lines+markers', hovertemplate = paste(
                             "%{y:.3f}<br>",
                             "Quantile(in %) : %{x:.0f}<br>",
                             "ExpressionS : %{text:.3f}<br>"
                           )) %>%
    add_trace(y = quantile_df$hr, name = 'Hazard Ratio',mode = 'lines+markers',type = 'scatter',
              line = list(color = 'rgb(212, 235, 242)', width = 2),
              marker=list(
                symbol = 'diamond',
                color=quantile_df$hr,
                colorscale='RdBu',
                cmid = 1,
                reversescale =FALSE
              ),
              yaxis = yaxes[2]
    )
  #different layout setting for first graph and second graph
  if(index == 1){
    fig <- fig %>%
      layout(title = 'P Values of Different Quantiles',
             xaxis = list(title = 'Quantile(%)'),
             yaxis = P_Axis,
             yaxis2 = HR_Axis,
             hovermode = "x unified"
      )
  }
    else{
      fig <- fig %>%
        layout(title = 'P Values of Different Quantiles',
             xaxis = list(title = 'Quantile(%)'),
             yaxis2 = P_Axis,
             yaxis3 = HR_Axis_special,
             hovermode = "x unified"
      )
      }
  
  return(fig)
}

assemble_percentile_plot <- function(quantile_df_list){
  fig_list <- c()
  for(index in 1:length(quantile_df_list)){
    fig_list[[index]] <- single_plot(quantile_df_list[[index]], index = index)
  }
  
  return(fig_list)
}

#Subplot----

test <- c()
test[["quantile"]][[1]] <- res$p_df
test[["quantile"]][[2]] <- res$p_df



fig_list <- assemble_percentile_plot(test$quantile)
subplot(fig_list[[1]], fig_list[[2]])


# projectss <- unname(unlist(projects[categories]))
# projectss


economics$uempmed
p <- subplot(
  plot_ly(economics, x = ~date, y = ~uempmed) %>% 
    layout(annotations = list(x = 0.2 , y = 1.05, text = "AA", showarrow = F,
                              xref='paper', yref='paper'), 
           showlegend = FALSE),
  plot_ly(economics, x = date, y = ~unemploy) %>% 
    layout(annotations = list(x = 0.2 , y = 1.05, text = "AA", showarrow = F, 
                              xref='paper', yref='paper'), 
           showlegend = FALSE),showlegend = FALSE)



#START OF CIRCLIZE----

graph_df <- readRDS("~/Downloads/circle_graph.rds")

df_crispr <- read.csv("~/Downloads/DepMap/DepMap-CRISPR.csv")
df_drug <- read.csv("~/Downloads/DepMap/DepMap-Drug.csv")
df_rnai <- read.csv("~/Downloads/DepMap/DepMap-RNAi.csv")


graph_df$survival_days <- as.numeric(graph_df$survival_days)

max_days <- max(graph_df$survival_days,na.rm = TRUE)
graph_df$survival_days[is.na(graph_df$survival_days)]  = max_days
rm(max_days)


#There are some columns in df_drug that have typeof as "character"
test <- sapply(df_drug[,2:length(colnames(df_drug))], function(x) as.numeric(x))
depmap_to_graph_df <- function(df,row_mean,subproject,name = colnames(graph_df)){
  result <- data.frame(df[,1]) 
  result <- cbind(result,"DepMap")
  result <- cbind(result,row_mean)
  result <- cbind(result,subproject)
  result <- cbind(result,"DepMap")
  colnames(result) <- name
  return(result)
}

result <- depmap_to_graph_df(df_rnai, rowMeans(df_rnai[,2:length(colnames(df_rnai))],na.rm = TRUE), "DepMap-RNAi")
graph_df <- rbind(result, graph_df)
result <- depmap_to_graph_df(df_drug, rowMeans(test,na.rm = TRUE), "DepMap-Drug")
graph_df <- rbind(result, graph_df)
result <- depmap_to_graph_df(df_crispr, rowMeans(df_crispr[,2:length(colnames(df_crispr))],na.rm = TRUE), "DepMap-CRISPR")
graph_df <- rbind(result, graph_df)
rm(result)



test <- colMeans(test,na.rm = TRUE)

#This is needed by DepMap boxplot
#Follows the order of Crispr, Drug, RNAi
DepMap <- list()
DepMap[[1]] <- colMeans(df_crispr[,2:length(colnames(df_crispr))],na.rm = TRUE)
DepMap[[2]] <- test
DepMap[[3]] <- colMeans(df_rnai[,2:length(colnames(df_rnai))],na.rm = TRUE)
names(DepMap[[1]]) <- NULL
names(DepMap[[2]]) <- NULL
names(DepMap[[3]]) <- NULL


#df_inner----
df_inner <- data.frame()
# newrow <- c("DepMap","DepMap-CRISPR",nrow(df_crispr))
# df_inner <- rbind(df_inner,newrow)
# newrow <- c("DepMap","DepMap-Drug",nrow(df_drug))
# df_inner <- rbind(df_inner,newrow)
# newrow <- c("DepMap","DepMap-RNAi",nrow(df_rnai))
# df_inner <- rbind(df_inner,newrow)
# rm(newrow)

unique_subproject <- unique(graph_df$subproject)

for(ele in unique_subproject){
  temp <- graph_df[graph_df$subproject == ele,]
  count <- nrow(temp)
  bigproject <- unique(temp$bigproject)
  newrow <- c(bigproject,ele,count)
  df_inner <- rbind(df_inner,newrow)
}
colnames(df_inner) <- c("bigproject","subproject","count")
df_inner$count <- as.numeric(df_inner$count)


#nrows are recorded here
# newrow <- c("DepMap","DepMap",nrow(df_crispr),"DepMap-CRISPR","DepMap")
# graph_df <- rbind(newrow,graph_df)
# newrow <- c("DepMap","DepMap",nrow(df_drug),"DepMap-Drug","DepMap")
# graph_df <- rbind(newrow,graph_df)
# newrow <- c("DepMap","DepMap",nrow(df_rnai),"DepMap-RNAi","DepMap")
# graph_df <- rbind(newrow,graph_df)

#Delete for testing purpose
#graph_df <- graph_df[4:(nrow(graph_df)),]
# takes in an array and returns an array where the max is mapped to 100 and the samllest is mapped to 10




#inner circle with mapping----
categories <- c("DepMap","TARGET","TCGA")
categories_cols = brewer.pal(3, "Set2")

#run CSurvival one time to get this
projectss <- unname(unlist(projects[categories]))
projectss

#Assign mapping index for each category/sub-category
df_inner <- cbind(df_inner,c(1:3,1:9,1:33))
df_inner <- cbind(df_inner,rep(1:3, c(3,9,33)))
colnames(df_inner)[(length(colnames(df_inner)) - 1) : (length(colnames(df_inner)))] <- c("subproject_mapping","bigproject_mapping")

#add a mapping column to graph dataframe
graph_df <- cbind(graph_df,unlist(lapply(1:length(graph_df$subproject), function(i) return(df_inner[df_inner$subproject == graph_df$subproject[i],"subproject_mapping"]))))
colnames(graph_df)[length(colnames(graph_df))] <- "subproject_mapping"



make_box_data <- function(graph_df,df_inner, bigproject){
  result <- list()
  for(index in 1:max(df_inner[df_inner$bigproject == bigproject,"subproject_mapping"])){
    temp = graph_df$survival_days[graph_df$subproject_mapping == index & graph_df$bigproject == bigproject]
    result[[index]] = temp 
  }
  return(result)
}

TCGA <- make_box_data(graph_df,df_inner,"TCGA")
TARGET <- make_box_data(graph_df,df_inner,"TARGET")




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

calculate_ylim <- function(graph_df,DepMap,bigproject = c("DepMap","TARGET","TCGA"),survival_days = "survival_days"){
  box_ylim <- matrix(nrow = 3,ncol = 2)
  survival_days <- rep(survival_days,3)
  box_ylim[1,] = c(min(unlist(DepMap)),max(unlist(DepMap)))
  box_ylim[2,] = c(min(graph_df[graph_df$bigproject == bigproject[2],survival_days[2]]),max(graph_df[graph_df$bigproject == bigproject[2],survival_days[2]]))
  box_ylim[3,] = c(min(graph_df[graph_df$bigproject == bigproject[3],survival_days[3]]),max(graph_df[graph_df$bigproject == bigproject[3],survival_days[3]]))
  row.names(box_ylim) <- bigproject
  return(box_ylim)
}

max(unlist(TCGA))
box_ylim
box_ylim <- calculate_ylim(graph_df,DepMap = DepMap)


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


rm(result)
rm(test)
rm(df_crispr)
rm(df_drug)
rm(df_rnai)



