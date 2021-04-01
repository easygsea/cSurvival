# extract gene expression/mutation data
extract_gene_data <- function(x, type){
  df_file <- list(
    "rna" = "df_gene_scale.csv"
    ,"lib" = "df_gene_scale.csv"
    ,"manual" = "df_gene_scale.csv"
  )
  # all genes in selected project
  a_range <- 2:(length(rv[[paste0("genes",x)]])+1)
  
  if(type == "rna"){
    all_genes <- rv[[paste0("genes",x)]]
    # selected gene
    g_ui_id <- paste0("g_",x)
    genes <- input[[g_ui_id]]
  }else if(type == "snv"){
    all_genes <- rv[[paste0("genes",x)]]
    # selected gene
    g_ui_id <- paste0("g_",x)
    genes <- input[[g_ui_id]]
    
    # detect mutation method
    if(rv$tcga){
      snv_id <- paste0("snv_method_",x)
      if(is.null(input[[snv_id]])){
        method <- rv[[snv_id]]
      }else{
        method <- input[[snv_id]]
      }
      df_file <- c(
        df_file
        ,"snv" = paste0("df_snv_class_",method,".csv")
      )
    }else{
      df_file <- c(
        df_file
        ,"snv" = paste0("df_snv_class.csv")
      )
    }
    
  }else if(type == "lib"){
    all_genes <- sapply(rv[[paste0("genes",x)]], function(x) toupper(strsplit(x,"\\|")[[1]][1])) %>% unname(.)
    genes <- toupper(rv[[paste0("gs_genes_",x)]])
  }else if(type == "manual"){
    all_genes <- sapply(rv[[paste0("genes",x)]], function(x) toupper(strsplit(x,"\\|")[[1]][1])) %>% unname(.)
    genes <- toupper(rv[[paste0("gs_m_",x)]])
  }
  
  # infile
  infile <- paste0(rv$indir,df_file[[type]])
  
  # # method 1 fread drop columns
  col_to_drop <- a_range[!all_genes %in% genes]
  data <- fread(infile,sep=",",header=T,drop = col_to_drop)
  
  # # # method 2 fread essential columns
  # ofile <- paste0(rv$indir,"tmp.csv")
  # unlink(ofile)
  # col_to_keep <- a_range[input[[g_ui_id]] == rv[[paste0("genes",x)]]]
  # system(paste0("cut -d',' -f1,",col_to_keep," ",infile," > ",ofile))
  # data <- fread(ofile,sep=",",header=T)

  return(data)
}

# original survival df
original_surv_df <- function(patient_ids){
  df_o <- rv$df_survival
  df_o %>% dplyr::filter(patient_id %in% patient_ids)
  df_o[match(patient_ids,df_o$patient_id),]
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
  return(df)
}

# generate survival df for analysis
get_info_most_significant_rna <- function(data, min, max, step, mode="g"){
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)
  
  # initialize the most significant p value and model
  least_p_value <- 1
  quantile_most_significant <- NULL
  
  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  if(mode == "g"){
    exp <-data[,2] %>% unlist(.) %>% unname(.)
  }else if(mode == "gs"){
    exp <- rowMeans(data[,-1]) %>% unlist(.) %>% unname(.)
  }
  
  # retrieve survival analysis df_o
  df_o <- original_surv_df(patient_ids)

  # the quantiles we will use to define the level of gene percentages
  quantiles <- quantile(exp, quantile_s)

  for(i in seq_along(quantiles)){
    q <- quantiles[i]
    df <- generate_surv_df(df_o, patient_ids, exp, q)

    # # test if there is significant difference between high and low level genes
    # if(rv$cox_km == "cox"){
      surv_diff <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
      p_diff <- coef(summary(surv_diff))[,5]
    # }else if(rv$cox_km == "km"){
    #   surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
    #   p_diff <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    # }
      
    if(p_diff <= least_p_value){
      least_p_value <- p_diff
      df_most_significant <- df
      cutoff_most_significant <- names(quantiles[i])
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
  exp <- data[,2] %>% unlist(.) %>% unname(.)
  
  # retrieve survival analysis df_o
  df_o <- original_surv_df(patient_ids)
  
  # the quantile we will use to define the level of gene percentages
  q <- quantile(exp, cutoff)
  df <- generate_surv_df(df_o, patient_ids, exp, q)
}

# generate df if SNV mutation data
get_df_snv <- function(data, nons){
  nons <- tolower(nons)
  # extract patients' IDs and mutation data
  patient_ids <- data$patient_id
  mutations <- data[,2] %>% unlist(.) %>% unname(.)
  # check if nonsyn
  mm <- mutations[!is.na(mutations)]
  mm <- sapply(mm, function(x){
    x <- strsplit(x, "\\|")[[1]]
    if(any(tolower(x) %in% nons)){
      "Nonsynonymous"
    }else{
      "Synonymous"
    }
  })
  mutations[!is.na(mutations)] <- mm
  mutations[is.na(mutations)] <- "Synonymous"
  data[,2] <- mutations
  
  # rename columns
  col_names <- colnames(data)
  col_names[length(col_names)] <- "level"
  colnames(data) <- col_names
  
  # reset levels
  lels <- unique(data$level) %>% sort(.,decreasing = T)
  data$level <- factor(data$level, levels = lels)

  # retrieve survival analysis df_o
  df_o <- original_surv_df(patient_ids)
  
  df <- df_o %>% inner_join(data, by = "patient_id")
}

# combine and generate interaction df

## Perform survival analysis
cal_surv_rna <- 
  function(
    df,n
    ,p.adjust.method = "hommel"
  ){
    # # 1. KM # #
    # summary statistics
    km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
    p.km <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
    
    # run KM
    km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df)
    
    # # 2. Cox # #
    # create new df to seperate effects in Cox regression
    lels <- levels(df$level)
    if(n == 1){
      new_df <- with(df,data.frame(level = lels))
      
      # run Cox regression
      cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
      # summary statistics
      cox.stats <- summary(cox_fit)
    }else if(n == 2){
      km2 <- pairwise_survdiff(Surv(survival_days, censoring_status) ~ level, data = df, p.adjust.method = p.adjust.method)
      km.stats <- list(km.stats,km2)
      
      lels_x <- levels(df$`level.x`)
      lels_y <- levels(df$`level.y`)
      new_df <- with(df,data.frame(level = lels, level.x = lels_x, level.y = lels_y))
      
      # run Cox regression
      cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level.x * level.y, data = df)
      # summary statistics
      cox.stats <- summary(cox_fit)
      # run Cox regression for visualization purpose
      cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
    }

    hr.cox <- sapply(cox.stats$coefficients[,2], function(x){
      round(as.numeric(x), 2)
    }) %>% paste0(.,collapse = ", ")
    p.cox <- sapply(cox.stats$coefficients[,5], function(x){
      format(as.numeric(x), scientific = T, digits = 3)
    }) %>% paste0(.,collapse = ", ")

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
        ,hr = hr.cox
        ,p = p.cox
      )
    )
    return(results)
  }

## Plot survival results
plot_surv <- 
  function(
    res, mode=rv$cox_km, two_rows="one"
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
                        surv.median.line = surv.median.line,            # Add median survival lines
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
      
      # adjust KM table size
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
    
    if(two_rows=="all"){
      fig <- fig + guides(col = guide_legend(nrow=2,byrow=TRUE))
    }
    return(fig)
  }
