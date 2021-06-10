install.packages("plotly")
install.packages("tidyverse")
install.packages("htmlwidgets")
library(circlize)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(colorspace)

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

#TODO: ADD COLORSCALE FOR HR PLOT----
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







#Circular Graph----
library(circlize)
#TEST QUICK GUIDE
set.seed(999)
n = 1000
df = data.frame(sectors = sample(letters[1:8], n, replace = TRUE),x = rnorm(n), y = runif(n))


circos.par("track.height" = 0.1)
circos.initialize(df$sectors, x = df$x)

circos.track(df$sectors, y = df$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)

bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)
circos.clear()



sectors = letters[1:3]
set.seed(999)
n = 1000
df = data.frame(sectors = sample(letters[1:3], n, replace = TRUE),x = rnorm(n), y = runif(n))

circos.par(points.overflow.warning = FALSE)
circos.initialize(df$sectors, x = df$x)
circos.trackPlotRegion(df$sectors, y = df$y,
                       track.height = 0.25, panel.fun = function(x, y) {
                         circos.text(median(df$x),median(df$y),"ABC", facing = "bending.outside", niceFacing = TRUE,
                                     cex = 1.2)
                       })
circos.update(sector.index = "c", track.index = 1)
circos.text(CELL_META$xcenter, CELL_META$ycenter, "DepMap", col = "black")

circos.clear()
