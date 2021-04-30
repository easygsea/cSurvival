# determine if plot_type is for a survival curve
if_surv <- function(plot_type=rv$plot_type){
  plot_type == "all" | suppressWarnings(!is.na(as.numeric(plot_type))) | plot_type == "gender"
}

#==============================================#
####        Survival analysis functions     ####
#==============================================#

# extract gene expression/mutation data
extract_gene_data <- function(x, type){
  g_ui_id <- paste0("g_",x)
  
  df_file <- list(
    "rna" = "df_gene.csv"
    ,"lib" = "df_gene_scale.csv"
    ,"manual" = "df_gene_scale.csv"
    ,"cnv" = "df_cnv.csv"
    ,"mir" = "df_mir.csv"
  )
  # # all genes in selected project
  # a_range <- 2:(length(rv[[paste0("genes",x)]])+1)
  
  if(type == "rna"){
    # original file that stores raw FPKM values
    all_file <- paste0(rv$indir,"df_gene.csv")[[1]]
    # # read in all genes
    # all_genes <- fread(all_file,sep=",",header=T,nrows=0) %>% names(.)
    # # all_genes <- sapply(all_genes, function(x){
    # #   x <- strsplit(x,"\\|")[[1]]
    # #   if(length(x) == 1){x}else{x[-1]}
    # # }) %>% unname(.)
    # all_genes <- all_genes[-1]
    # 
    # # change a_range
    # a_range <- 2:(length(all_genes)+1)
    
    # selected gene
    genes <- input[[g_ui_id]]
    # extract ENSG info
    genes <- strsplit(genes,"\\|")[[1]]
    if(length(genes) == 1){genes <- genes}else{genes <- tail(genes,n=1)}
  }else if(type == "snv"){
    # selected gene
    genes <- input[[g_ui_id]]
    
    # detect mutation method
    if(rv$tcga){
      snv_id <- paste0("snv_method_",x)
      method <- ifelse_rv(snv_id)
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
    genes <- rv[[paste0("genes",x)]][all_genes %in% genes]
  }else if(type == "manual"){
    all_genes <- sapply(rv[[paste0("genes",x)]], function(x) toupper(strsplit(x,"\\|")[[1]][1])) %>% unname(.)
    genes <- toupper(rv[[paste0("gs_m_",x)]])
    genes <- rv[[paste0("genes",x)]][all_genes %in% genes]
  }else if(type == "cnv" | type == "mir"){
    # all_genes <- rv[[paste0("genes",x)]]
    genes <- input[[g_ui_id]]
  }
  
  # infile
  infiles <- paste0(rv$indir,df_file[[type]])
  
  # # method 1 fread drop columns
  # col_to_drop <- a_range[!all_genes %in% genes]
  l <- lapply(infiles,function(y){
    fread(y,sep=",",header=T,select = c("patient_id", genes))
  })
  data <- rbindlist(l, use.names = T)

  # # # method 2 fread essential columns
  # ofile <- paste0(rv$indir,"tmp.csv")
  # unlink(ofile)
  # col_to_keep <- a_range[input[[g_ui_id]] == rv[[paste0("genes",x)]]]
  # system(paste0("cut -d',' -f1,",col_to_keep," ",infile," > ",ofile))
  # data <- fread(ofile,sep=",",header=T)
  
  # save original expression or mutation data, if applicable
  if(type == "rna" | type == "mir"){
    # save original expression data
    rv[[paste0("exprs_",x)]] <- data
  }else if(type == "snv"){
    muts <- data[,2] %>% unlist(.)
    names(muts) <- data$patient_id
    muts <- muts[!is.na(muts)]
    rv[[paste0("mutations_",x)]] <- muts
  }
  if(type == "lib" | type == "manual"){
    # save original FPKM data
    if(x == 1){rv[[paste0("exprs_",x)]] <- data}
    
    # z score transform expression values
    n_col <- ncol(data)
    exp_scale <- apply(data[,2:n_col], 2, scale)
    data <- cbind(data[,1,drop=F],exp_scale)
  }
  
  return(data)
}

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

# generate survival df for analysis
get_info_most_significant_rna <- function(data, min, max, step, mode="g"){
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)
  
  # initialize the most significant p value
  least_p_value <- 1

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
    if(!is.na(p_diff)){
      if(p_diff <= least_p_value){
        least_p_value <- p_diff
        df_most_significant <- df
        cutoff_most_significant <- names(quantiles[i])
      }
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

# generate df if CNV copy number data
get_info_most_significant_cnv <- function(data, mode){
  # initialize the most significant p value
  least_p_value <- 1

  # extract patients' IDs and expression values
  patient_ids <- data$patient_id

  # retrieve survival analysis df_o
  df_o <- original_surv_df(patient_ids)

  # copy numer data
  exp <-data[,2] %>% unlist(.) %>% unname(.)
  
  # threshold for different programs
  if(rv$tcga){
    threshold <- 0
  }else{
    threshold <- 2
  }

  # if(mode != "both"){
    # loop between loss and gain, if not both
    if(mode == "auto"){
      if(rv$tcga){cats <- c(-1,1)}else{cats <- c(1,3)}
      names(cats) <- c("Loss","Gain")
    }else if(mode == "gain"){
      if(rv$tcga){cats <- 1}else{cats <- 3}
      names(cats) <- "Gain"
    }else if(mode == "loss"){
      if(rv$tcga){cats <- -1}else{cats <- 1}
      names(cats) <- "Loss"
    }
    
    
    for(cat in seq_along(cats)){
      i <- as.numeric(cats[[cat]])
      if(i > threshold){
        lells <- ifelse(exp >= i, "Gain", "Other")
      }else{
        lells <- ifelse(exp <= i, "Loss", "Other")
      }
      lels <- unique(lells) %>% sort(.,decreasing = T)
      df <- df_o
      df$level <- factor(lells, levels = lels)
      
      # # test if there is significant difference between high and low level genes
      # if(rv$cox_km == "cox"){
      surv_diff <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
      p_diff <- summary(surv_diff)$logtest[3] #coef(summary(surv_diff))[,5]
      # }else if(rv$cox_km == "km"){
      #   surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
      #   p_diff <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      # }
      if(p_diff <= least_p_value){
        least_p_value <- p_diff
        df_most_significant <- df
        cutoff_most_significant <- names(cats[cat])
      }
    }
  # }
  
  # # additional analysis if to look at both gain and loss
  # if(mode == "auto" | mode == "both"){
  #   lells <- ifelse(exp > threshold, "Gain", ifelse(exp < threshold, "Loss", "Other"))
  #   lels <- unique(lells) %>% sort(.,decreasing = T)
  #   df <- df_o
  #   df$level <- factor(lells, levels = lels)
  #   surv_diff <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
  #   p_diff <- summary(surv_diff)$logtest[3] #min(coef(summary(surv_diff))[,5]) 
  #   if(p_diff <= least_p_value){
  #     least_p_value <- p_diff
  #     df_most_significant <- df
  #     cutoff_most_significant <- "Gain & Loss"
  #   }
  # }
  
  results <- list(
    df = df_most_significant,
    cutoff = cutoff_most_significant
  )
  return(results)
}

## Perform survival analysis
cal_surv_rna <- 
  function(
    df,n
    ,p.adjust.method = rv[["km_mul"]]
  ){
    # # 1. KM # #
    # summary statistics
    km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
    p.km <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
    
    # run KM
    km.fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df)
    
    # # 2. Cox # #
    if(n == 1){
      # run Cox regression
      cox_fit <- coxph(Surv(survival_days, censoring_status) ~ level, data = df)
      # summary statistics
      cox.stats <- summary(cox_fit)
    }else if(n == 2){
      km2 <- try(pairwise_survdiff(Surv(survival_days, censoring_status) ~ level, data = df, p.adjust.method = p.adjust.method))
      if(inherits(km2, "try-error")) {
        rv$try_error <- rv$try_error + 1
        shinyalert(paste0("The selected two genes/loci have exactly the same data."
                          ," If CNV, you might have selected genes from the same cytoband."
                          ," Please contrast genes from different cytobands."))
      }
      req(!inherits(km2, "try-error")) #require to be no error to proceed the following codes
      
      km.stats <- list(km.stats,km2)

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
    lels_x <- levels(df$`level.x`)
    lels_y <- levels(df$`level.y`)
    # lels <- apply(expand.grid(lels_x,lels_y),1,paste0,collapse="_")
    # create new df to seperate effects in Cox regression
    lels <- levels(df$level)
    new_df <- with(df,data.frame(level = lels))
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
    req(!is.null(fit))
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
      fig <- fig + guides(col = guide_legend(nrow=rv$variable_nr,byrow=TRUE))
    }else if(two_rows=="gender"){
      fig <- fig + guides(col = guide_legend(nrow=2,byrow=TRUE))
    }
    return(fig)
  }

#==============================================#
####            eVITTA functions            ####
#==============================================#
# the gene expression table
de_dfgene <- function(){
  # patients
  patients <- rv[["df_gender"]][["patient_id"]]
  
  # read in data
  infiles <- paste0(rv$indir,"df_gene.csv")
  
  l <- lapply(infiles, function(x){
    info <- fread(x,sep=",",header=T)
    return(info)
  })
  
  # whole gene expression table
  df_gene <- rbind_common(l) %>% dplyr::filter(patient_id %in% patients)
  
  # remaining patient ids
  patients <- df_gene$patient_id
  # the genes
  genes <- colnames(df_gene)[-1]
  
  incProgress(amount = 0.1, message = wait_msg("Formatting gene expression matrix for DE analysis..."))
  
  # transpose
  df_gene <- transpose(df_gene) %>% .[-1,] %>% as.data.frame()
  # reassign patient ids
  colnames(df_gene) <- patients
  rownames(df_gene) <- genes
  
  incProgress(amount = 0.2, message = wait_msg("Converting gene IDs..."))
  
  # id conversion
  # create individual tables using org.Hs
  egENS <- toTable(org.Hs.egENSEMBL)
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  
  # bind the tables
  id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
  
  # extract the gene ids from df_gene and find their names from id conversion table
  gene_ids <- as_tibble_col(rownames(df_gene), column_name = "ensembl_id")
  gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")
  
  # convert gene ids
  df_gene <- df_gene %>% dplyr::mutate(ensembl_id=rownames(df_gene)) %>%
    left_join(gene_ids_table, by = "ensembl_id")
  
  # remove unnecessary columns
  df_gene <- df_gene %>%
    dplyr::distinct(symbol,.keep_all = TRUE) %>%
    dplyr::filter(!is.na(symbol))
  
  genes <- df_gene$symbol
  
  # convert to numeric matrix
  df_gene <- df_gene %>%
    dplyr::select(-c(ensembl_id, gene_id, symbol)) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix(.)

  rownames(df_gene) <- genes
  
  return(df_gene)
}
