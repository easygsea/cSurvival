library(ggplot2)
library(colorspace)
setwd("/Applications/Codes/cSurvival_dev/basic_scripts")
df <- readRDS("./highlight_plot/df.rds")
p <- readRDS("./highlight_plot/plot.rds")
annot_cells <- readRDS("./highlight_plot/annot_cell.rds")
pos <- readRDS("./highlight_plot/pos.rds")
cols <- readRDS("./highlight_plot/cols.rds")
lels_len <- length(levels(df$Level))


# change heximal colors to 90% transparency
addalpha <- function(colors, alpha=0.35) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}




    

      #CHANGE PART
  #Outside of if statement. Get color set up properly
  cols <- c(cols,addalpha(cols),darken(cols, 0.4))
  names(cols) <- c(levels(df$Level), paste0(levels(df$Level), ",Not highlighted"), paste0(levels(df$Level), ",Highlighted"))
  cols
  
  
  
  
  
  
#  names(cols) <- 
  #Inside of if statement. Get df$annotation set up properly

  
  df$Annotation <- paste0(df$Level,",",ifelse(df$Cell %in% annot_cells, "Highlighted","Not highlighted"))
  
  p <- p + 
          geom_jitter(position=pos, aes(colour = ifelse(df$Cell %in% annot_cells, "Highlighted", "Not highlighted"), size = ifelse(df$Cell %in% annot_cells, "Highlighted", "Not highlighted"))) +
          scale_fill_manual(values=c(cols,"#999999", "#E69F00")) + 
          scale_size_manual(values=c(3,1)) +
          #scale_shape_manual(values=c(8,16))
          labs(color = "", shape = "", size = "")
      
    ggplotly(p,tooltip = c("Line","x","y"))

    
    
    
    
    
    
    #CHECK VIOLIN
    
    if(length(input$annot_data_points)>0){
      if(input$annot_data_points != ""){
        df$Annotation <- paste0(df$mut_cat,",",ifelse(df$patient_id %in% input$annot_data_points, "Highlighted", "Not highlighted"))
        cat_name <- levels(df$mut_cat)
        cat_name <- cat_name[cat_name != "Other"]
        cols <- c("#00BFC4", "#F8766D", "#5a696a", addalpha("#00BFC4"), "#f41e0f", addalpha("#F8766D"))
        names(cols) <- c("Other",cat_name,"Other,Highlighted","Other,Not highlighted",paste0(cat_name,",Highlighted"),paste0(cat_name,",Not highlighted"))
      }else{
        df$Annotation <- df$mut_cat
        cols <- c("#00BFC4", "#F8766D")
      }
    }else{
      df$Annotation <- df$mut_cat
      cols <- c("#00BFC4", "#F8766D")
    }
    
    
    
    
    
    
    
    
    
#Experiment HEATMAP
min = 0.1
max = 0.9
step = 0.05

populate_quantile_df <- function(min, max, step, substitue_value = 1){
  temp <- lapply(quantile_s, function(i){
    rep(i, length(quantile_s))
  })
  temp <- unlist(temp)  
  temp_second_col <- rep(quantile_s, length(quantile_s))
  result <- data.frame(temp, temp_second_col, substitue_value)
  colnames(result) <- c("Q1","Q2", "p_value")
  result <- data.frame(sapply(result, function(x) as.numeric(as.character(x))))
  result$Q1 = result$Q1 * 100
  result$Q2 = result$Q2 * 100
  return(result)
}
quantile_s = seq(min, max, by = step)
typeof(result$Q1)

heatmap_df = populate_quantile_df(min, max, step, 1)



#draw heatmap function
heatmap_df <- readRDS("/Applications/Codes/cSurvival_dev/heatmap_df.rds")
col_scale_no <- c(0, 0.16666666, 0.33333333333, 0.5, 0.66666666, 1)
col_scale <- sequential_hcl(5, palette = "YlOrRd") %>% rev(.) %>% col2rgb()
col_scale <- sapply(1:ncol(col_scale), function(x) {x <- col_scale[,x]; paste0("rgb(",paste0(x, collapse = ", "),")")})
col_scale <- c("rgb(255, 255, 255)", col_scale)
col_scale <- lapply(1:length(col_scale), function(i) list(col_scale_no[i], col_scale[i]))
pvalue_heatmap <- function(heatmap_df){
  dat = heatmap_df
  dat$log_p_value <- -log10(dat$p_value)
  dat$log_p_value[is.na(dat$log_p_value)] <- 0
  
  req(length(dat$log_p_value)>0)
  
  fig <- plot_ly() %>%
    add_trace(data = dat, x = ~Q1, y = ~Q2, z = ~log_p_value, type = "heatmap",
              colorscale  = col_scale,zmax = max(dat$log_p_value),zmin=0,
              colorbar = list(
                title = list(text="-log10(P)", side = "right")
                ,len = 1.2
              ),
              text = ~p_value,
              hovertemplate = paste('<b>%{x}</b> vs <b>%{y}</b><br>',
                                    'Corresponding P-value: <b>%{text}</b>'
              )
    )%>%
    add_annotations(x = dat$Q1, y = dat$Q2,
                    text = as.character(round(dat$p_value,2)),
                    showarrow = FALSE, xref = 'x', yref = 'y', font=list(color='black', size = 10)
                    ,ax = 20, ay = -20)
    #TODO: Add y axis and x axis label
    # %>% 
    # layout(
    #   # title = "Pariwise comparisons",
    #   xaxis = list(title = paste0("Pariwise comparisons adjusted by ",mul_methods), showticklabels = T),
    #   yaxis = list(title = "", showticklabels = T)
    #   # ,margin = list(l=200)
    # )
  fig
}
pvalue_heatmap(heatmap_df)
# 
# #String manipulation:
# string = "ABCD|100 expression"
# 
# str_remove(string, fixed(" expression"))

    