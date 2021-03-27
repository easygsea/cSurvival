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
    df <- generate_surv_df(patient_ids, exp, quantiles[i])

    # test if there is significant difference between high and low level genes
    surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
    p_diff <- 1- pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    if(p_diff < least_p_value){
      least_p_value = p_diff
      df_most_significant <- df
      quantile_most_significant <- names(quantiles)[i]
      cutoff_most_significant <- quantiles[i]
    }
  }
  results <- list(
    df = df_most_significant,
    pval = least_p_value,
    quantile = quantile_most_significant,
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

## Perform survival analysis
cal_surv_rna <- function(df){
  km_fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df)
  return(km_fit)
}
