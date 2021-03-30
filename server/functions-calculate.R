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
    sapply(function(x) ifelse(x > q, "High", "Low"))
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
cal_surv_rna <- 
  function(
    df
  ){
    # # 1. KM # #
    # summary statistics
    km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = data)
    p.km <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
    
    # run KM
    km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df)
    
    # # 2. Cox # #
    # create new df to seperate effects in Cox regression
    lels <- levels(df$level)
    new_df <- with(df,data.frame(level = lels))
    
    # run Cox regression
    cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
    
    # summary statistics
    cox.stats <- summary(cox_fit)
    hr.cox <- cox.stats$coefficients[,2]
    p.cox <- cox.stats$coefficients[,5]

    # run Cox survival analysis
    cox.fit <- survfit(cox_fit,newdata=new_df)

    # save df, fit, and statistics
    results <- list(
      km = list(
        df = df,
        fit = km.fit,
        stats = km.stats
        ,lels = lels
        ,hr = "NA"
        ,p = format(as.numeric(p.km), scientific = T, digits = 3)
      )
      ,cox = list(
        df = new_df,
        fit = cox.fit,
        stats = cox.stats
        ,lels = lels
        ,hr = round(as.numeric(hr.cox), 2)
        ,p = format(as.numeric(p.cox), scientific = T, digits = 3)
      )
    )
    return(results)
  }

## Plot survival results
plot_surv <- 
  function(
    res, mode=rv$cox_km
    , title=NULL
    , risk.table = rv$risk_table, cumevents = rv$cum_table, ncensor.plot = FALSE # parameters for KM mode
    , conf.int=rv$confi, conf.int.style = rv$confi_opt# "ribbon" "step"
    # , surv.median.line="none" # "hv", "h", "v"
    , palette="jco"
    , base_size=20
  ){
    df <- res[[mode]][["df"]]
    fit <- res[[mode]][["fit"]]
    lels <- res[[mode]][["lels"]]
    
    # median survival lines
    if(is.null(rv$median)){
      surv.median.line="none"
    }else if(length(rv$median)==2){
      surv.median.line="hv"
    }else{
      surv.median.line=rv$median
    }
    
    if(mode == "km"){
      fig <- ggsurvplot(fit, data=df, 
                        title = title,
                        xlab = "Days",
                        ylab = "Survival probability",
                        conf.int=conf.int, conf.int.style=conf.int.style,
                        risk.table = risk.table, 
                        cumevents = cumevents,                   # Add cumulative No of events table
                        ncensor.plot = ncensor.plot,
                        palette = palette,
                        legend.title = "",               # Change legend titles
                        legend.labs = lels,  # Change legend labels
                        ggtheme = theme_survminer(
                          base_size = base_size,
                          font.main = c((base_size + 2), "plain", "black"),
                          font.submain = c(base_size, "plain", "black"),
                          font.x = c(base_size, "plain", "black"),
                          font.y = c(base_size, "plain", "black"),
                          font.caption = c(base_size, "plain", "black"),
                          font.tickslab = c((base_size - 2), "plain", "black"),
                          # legend = c("top", "bottom", "left", "right", "none"),
                          font.legend = c(base_size, "plain", "black")
                        ),
                        tables.height = 0.15,               # Specify tables height
                        tables.theme = theme_cleantable(),  # Clean theme for tables
                        tables.y.text = FALSE               # Hide tables y axis text
      )
      
      # adjust Cox table size
      base_size2 <- base_size
      fig$table <- fig$table + theme_cleantable(
        base_size = base_size2,
        font.main = c(base_size2, "plain", "black"),
        font.submain = c(base_size2, "plain", "black"),
        font.caption = c(base_size2, "plain", "black"),
        font.tickslab = c((base_size2 - 2), "plain", "black"),
        # legend = c("top", "bottom", "left", "right", "none"),
        font.legend = c(base_size2, "plain", "black")
      )
    }else if(mode == "cox"){
      fig <- ggsurvplot(fit,data=df,
                        title = title,
                        xlab = "Days",
                        ylab = "Survival probability",
                        conf.int=conf.int, conf.int.style=conf.int.style,
                        surv.median.line = surv.median.line,            # Add median survival lines
                        legend.title = "",               # Change legend titles
                        legend.labs = lels,  # Change legend labels
                        ggtheme = theme_survminer(
                          base_size = base_size,
                          font.main = c((base_size + 2), "plain", "black"),
                          font.submain = c(base_size, "plain", "black"),
                          font.x = c(base_size, "plain", "black"),
                          font.y = c(base_size, "plain", "black"),
                          font.caption = c(base_size, "plain", "black"),
                          font.tickslab = c((base_size - 2), "plain", "black"),
                          # legend = c("top", "bottom", "left", "right", "none"),
                          font.legend = c(base_size, "plain", "black")
                        ),
                        palette = palette                    # Use JCO journal color palette
      )
    }
    
    return(fig)
  }
