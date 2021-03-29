# extract gene expression/mutation data
extract_gene_data <- function(x, type){
  if(type == "rna"){
    g_ui_id <- paste0("g_",x)
    a_range <- 2:(length(rv[[paste0("genes",x)]])+1)
    col_to_drop <- a_range[input[[g_ui_id]] != rv[[paste0("genes",x)]]]
    
    data <- fread(paste0(rv$indir,"df_gene_scale.csv"),sep=",",header=T,drop = col_to_drop)
    return(data)
  }
}

# generate survival df
generate_surv_df <- function(patient_ids, exp, q){
  # generate the data from for model first
  gene_quantiles <- exp %>% 
    sapply(function(x) ifelse(x > q, "high", "low"))
  names(gene_quantiles) <- patient_ids
  
  # generate survival analysis df
  df <- rv$df_survival
  df$level <- gene_quantiles[match(df$patient_id,names(gene_quantiles))]
  df %>% dplyr::filter(level != "NULL")
  lels <- unique(df$level) %>% sort(.,decreasing = T)
  df$level <- factor(df$level, levels = lels)
  return(df)
}

# generate survival df for analysis
get_info_most_significant_rna <- function(data, min, max, step){
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)

  # initialize the most significant p value and model
  least_p_value <- 1
  quantile_most_significant <- NULL
  
  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  exp <-data[,2] %>% unlist(.) %>% unname(.)

  # the quantiles we will use to define the level of gene percentages
  quantiles <- quantile(exp, quantile_s)
  
  for(i in seq_along(quantiles)){
    q <- quantiles[i]
    df <- generate_surv_df(patient_ids, exp, q)

    # test if there is significant difference between high and low level genes
    surv_diff <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
    p_diff <- coef(summary(surv_diff))[,5]
    if(p_diff < least_p_value){
      # least_p_value = p_diff
      df_most_significant <- df
      cutoff_most_significant <- names(quantiles)[i]
    }
  }
  lels <- unique(df_most_significant$level) %>% sort(.,decreasing = T)
  df_most_significant$level <- factor(df_most_significant$level, levels = lels)
  
  results <- list(
    df = df_most_significant,
    cutoff = cutoff_most_significant
  )
  return(results)
}

# generate df if user-defined cutoffs
get_df_by_cutoff <- function(data, cutoff){
  cutoff <- cutoff/100
  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  exp <-data[,2] %>% unlist(.) %>% unname(.)
  
  # the quantile we will use to define the level of gene percentages
  q <- quantile(exp, cutoff)
  df <- generate_surv_df(patient_ids, exp, q)
}

# combine and generate interaction df

## Perform survival analysis
cal_surv_rna <- function(df, conf.int=T, surv.median.line="none", palette="jco"){
  # run KM
  km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df)
  km.surv <- ggsurvplot(km.fit, data=df, risk.table = TRUE, palette = palette)
  
  # create new df to seperate effects
  lels <- unique(df$level) %>% sort(.,decreasing = T)
  new_df <- with(df,data.frame(level = lels))
  
  # run Cox regression
  cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
  cox.fit <- survfit(cox_fit,newdata=new_df)
  cox.surv <- ggsurvplot(cox.fit,data=new_df,
                         title = "Survival Curves",
                         xlab = "Days",
                         ylab = "Survival probability",
                         conf.int=conf.int,
                         surv.median.line = surv.median.line,            # Add median survival lines
                         # legend.title = call_datatype(x),               # Change legend titles
                         legend.labs = lels,  # Change legend labels
                         palette = palette,                    # Use JCO journal color palette
                         risk.table.height = 0.3
                         # risk.table = T,                  # Add No at risk table
                         # cumevents = TRUE,                   # Add cumulative No of events table
                         # tables.height = 0.15,               # Specify tables height
                         # tables.theme = theme_cleantable(),  # Clean theme for tables
                         # tables.y.text = FALSE               # Hide tables y axis text
  )
  
  # add KM table to Cox table
  cox.surv$table <- km.surv$table
  fig <- cox.surv
  
  results <- list(
    stats = summary(cox_fit),
    fig = fig
  )
  return(results)
}
