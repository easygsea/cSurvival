install.packages("plotly")
install.packages("tidyverse")
install.packages("htmlwidgets")
library(plotly)
library(tidyverse)
library(htmlwidgets)

#SET UP----
#read data.rds
data <- readRDS("~/Downloads/rstudio-export/rna_data.rds")
df_o <- readRDS("~/Downloads/rstudio-export/df_o.rds")
rvdepmap <- readRDS("~/Downloads/rstudio-export/rv_depmap.rds")



#Helper function
# original survival df
original_surv_df <- function(patient_ids){
  df_o <- rv$df_survival
  df_o %>% dplyr::filter(patient_id %in% patient_ids)
  df_o[match(patient_ids,df_o$patient_id),] %>%
    dplyr::select(-gender)
}


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
    new_row = c(p_diff,unlist(strsplit(names(quantiles[i]),split = '%',fixed=T)),quantiles[i])
    p_df <- rbind(p_df,new_row)
    # }
    if(!is.na(p_diff)){
      if(p_diff <= least_p_value){
        least_p_value <- p_diff
        df_most_significant <- df
        cutoff_most_significant <- names(quantiles[i])
      }
    }
  }
  #Transform the p_df a little bit to make it work with the ggplot
  colnames(p_df) = c('p_value','quantile','expression')
  p_df$p_value = as.numeric(p_df$p_value)
  p_df$quantile = as.numeric(p_df$quantile)
  p_df$expression = as.numeric(p_df$expression)
  # proceed only if enough data
  if(is.null(df_most_significant)){
    return(NULL)
  }else{
    results <- list(
      df = df_most_significant,
      cutoff = cutoff_most_significant,
      p_df = p_df
    )
    return(results)
  }
}

res <- get_info_most_significant_rna(data,min = 0.2,max = 0.8, step = 0.01, df_o = df_o)
res$p_df


# 
# graph <- ggplot(data=res$p_df, aes(x = res$p_df$quantile, y = res$p_df$p_value, group=1)) +
#   geom_line(linetype = "dashed")+
#   geom_point()+labs(x = "Quantile(%)", y = "P Value") +
#   scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75,1))+
#   scale_x_continuous(limits=c(20, 80), breaks=c(20,40,60,80))
# 
# graph <- ggplotly(graph) %>% 
#   layout(hovermode = "x unified")
# graph
#   
#custom_colorscale = brewer.pal(n = length(res$p_df$p_value), "YlOrRd")[1,length(res$p_df$p_value)]

fig <- plot_ly(res$p_df, x = res$p_df$quantile, line = list(color = 'rgb(147,149,151)', width = 2))
fig <- fig %>% add_trace(y = ~res$p_df$p_value,type = 'scatter',#color =~p_value,
                         marker=list(
                           color=~p_value,
                           colorbar=list(
                             title='Colorbar'
                           ),
                           colorscale='YlOrRd',#custom_colorscale,
                           reversescale =FALSE
                         ),
                         text = res$p_df$expression,
                         
                         #colors = brewer.pal("YlOrRd "),
                         #rev('YlOrRd'),#brewer.pal(length(res$p_df$p_value),"YlOrRd "),
  name = '',mode = 'lines+markers', hovertemplate = paste(
  "P value is : %{y:.3f}<br>",
  "Quantile(in %) is: %{x:.0f}<br>",
  "Expression is : %{text:.3f}<br>"
)) %>%
layout(title = 'P Values of Different Quantiles',
       xaxis = list(title = 'Quantile(%)'),
       yaxis = list (title = 'P Value'))

fig


#hovertemplate = paste('P Value: $%{y:.2f}','<br>Quantile: %{x}<br>')



