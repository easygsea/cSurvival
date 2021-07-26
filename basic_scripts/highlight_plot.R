library(ggplot2)
df <- readRDS("./highlight_plot/df.rds")
p <- readRDS("./highlight_plot/plot.rds")
annot_cells <- readRDS("./highlight_plot/annot_cell.rds")
pos <- readRDS("./highlight_plot/pos.rds")
cols <- readRDS("./highlight_plot/cols.rds")

# change heximal colors to 90% transparency
addalpha <- function(colors, alpha=0.35) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

    

      #CHANGE PART
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