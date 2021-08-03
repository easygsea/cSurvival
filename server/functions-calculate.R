# determine if plot_type is for a survival curve
if_surv <- function(plot_type=rv$plot_type){
  plot_type == "all" | suppressWarnings(!is.na(as.numeric(plot_type))) | plot_type == "gender"
}
# determine if to plot forest plot
if_forest <- function(){
  gp <- rv$risk_gpr
  rv$cox_km == "cox" & (!(gp == "All" & (rv$plot_type == "all" | rv$plot_type == "gender")))
}
# determine if risk group need to be defined
if_any_exp <- function(){
  cats <- sapply(1:rv$variable_n, function(x){input[[paste0("cat_",x)]]}) %>% unique()
  dbs <- sapply(1:rv$variable_n, function(x){input[[paste0("db_",x)]]}) %>% unique()
  if_g <- length(cats)==1 & ("g" %in% cats)
  if_exp <- !((length(dbs)==1 & (("snv" %in% dbs) | (!rv$depmap & any(c("snv","cnv") %in% dbs))))|(length(dbs)>1 & !rv$depmap & all(c("snv","cnv") %in% dbs)))
  ("gs" %in% cats) | (if_g & if_exp)
}
#==============================================#
####        Survival analysis functions     ####
#==============================================#

# extract gene expression/mutation data
extract_gene_data <- function(x, type){
  g_ui_id <- paste0("g_",x)

  df_file <- list(
    "rna" = "df_gene.csv"
    ,"lib" = "df_gene.csv"
    ,"manual" = "df_gene_scale.csv"
    ,"cnv" = "df_cnv.csv"
    ,"mir" = "df_mir.csv"
    ,"met" = "df_met.csv"
    ,"pro" = "df_proteomic.csv"
    ,"rrpa" = "df_rrpa.csv"
    ,"crispr" = "DepMap-CRISPR.csv"
    ,"rni" = "DepMap-RNAi.csv"
    ,"drug" = "DepMap-Drug.csv"
  )
  # # all genes in selected project
  # a_range <- 2:(length(rv[[paste0("genes",x)]])+1)
  all_genes <- sapply(rv[[paste0("genes",x)]], function(x) toupper(strsplit(x,"\\|")[[1]][1])) %>% unname(.)

  if(type == "rna"){
    # # original file that stores raw FPKM values
    # all_file <- paste0(rv$indir,"df_gene.csv")[[1]]
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

    # # selected gene
    genes <- input[[g_ui_id]]
    # if(rv$target){
    #   # extract ENSG info
    #   genes <- strsplit(genes,"\\|")[[1]]
    #   if(length(genes) == 1){genes <- genes}else{genes <- tail(genes,n=1)}
    # }
    # # selected gene for normalization purpose
    g_ui_norm_id <- paste0("gnorm_",x); g_ui_norm_g_id <- paste0("gnorm_g_",x); g_ui_norm_gs_lib_id <- paste0("gnorm_gs_lib_",x)
    if(input[[g_ui_norm_id]] == "g"){
      genes_n <- input[[g_ui_norm_g_id]]
    }else if(input[[g_ui_norm_id]] == "gs"){
      genes_n <- toupper(rv[[paste0("gnorm_gs_genes_",x)]])
      genes_n <- rv[[paste0("genes",x)]][all_genes %in% genes_n]
    }
  }else if(type == "snv"){
    # selected gene
    genes <- input[[g_ui_id]]

    # detect mutation method
    if(rv$tcga){
      # snv_id <- paste0("snv_method_",x)
      # method <- ifelse_rv(snv_id)

      # df_file[["snv"]] = paste0("df_snv_class_",method,".csv")
      df_file[["snv"]] = paste0("df_snv_class_977.csv")

    }else{
      df_file[["snv"]] = paste0("df_snv_class.csv")
    }

  }else if(type == "lib"){
    genes <- toupper(rv[[paste0("gs_genes_",x)]])
    genes <- rv[[paste0("genes",x)]][all_genes %in% genes]

    g_ui_norm_id <- paste0("gnorm_",x); g_ui_norm_g_id <- paste0("gnorm_g_",x); g_ui_norm_gs_lib_id <- paste0("gnorm_gs_lib_",x)
    if(input[[g_ui_norm_id]] == "g"){
      genes_n <- input[[g_ui_norm_g_id]]
    }else if(input[[g_ui_norm_id]] == "gs"){
      genes_n <- toupper(rv[[paste0("gnorm_gs_genes_",x)]])
      genes_n <- rv[[paste0("genes",x)]][all_genes %in% genes_n]
    }
  }else if(type == "manual"){
    genes <- toupper(rv[[paste0("gs_m_",x)]])
    genes <- rv[[paste0("genes",x)]][all_genes %in% genes]
    rv[[paste0("gs_m_len",x)]] <- length(genes)
  }else if(type == "cnv" | type == "mir" | type == "met" | type == "pro" | type == "rrpa" | type == "crispr" | type == "rnai" | type == "drug"){
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
  # # read in normalization genes
  if(exists("genes_n")){
    l_n <- lapply(infiles,function(y){
      fread(y,sep=",",header=T,select = c("patient_id", genes_n))
    })
  }else{
    l_n <- NULL
  }
  # if(type == "snv" & rv$tcga){
  #   data <- Reduce(
  #     function(x, y) inner_join(x, dplyr::select(y, patient_id, genes), by = "patient_id"),
  #     l
  #   )
  #   uni_mode <- ifelse_rv(paste0("snv_uni_",x))
  #   data[[genes]] <- apply(data[ ,2:ncol(data)] , 1 , function(x) common_mut(x,mode=uni_mode))
  #   data <- data %>% dplyr::select(patient_id, genes)
  # }else{
    data <- rbindlist(l, use.names = T)
    if(!is.null(l_n)){
      data_n <- rbindlist(l_n, use.names = T)
    }else{
      data_n <- NULL
    }
    # ## remove NA data
    # if(!(type == "lib" | type == "manual")){
    #   data_s <- unlist(data[,2])
    #   data <- data[!is.na(data_s) & data_s != "",]
    # }
  # }

  # if depmap, filter patient ID
  if(rv$depmap){
    data <- data %>% dplyr::filter(patient_id %in% input$ccle_cells)
    if(!is.null(data_n)){
      data_n <- data_n %>% dplyr::filter(patient_id %in% input$ccle_cells)
    }
  }
  # # # method 2 fread essential columns
  # ofile <- paste0(rv$indir,"tmp.csv")
  # unlink(ofile)
  # col_to_keep <- a_range[input[[g_ui_id]] == rv[[paste0("genes",x)]]]
  # system(paste0("cut -d',' -f1,",col_to_keep," ",infile," > ",ofile))
  # data <- fread(ofile,sep=",",header=T)

  # save original expression or mutation data, if applicable
  if(type == "rna" | type == "mir" | type == "met" | type == "rrpa" | type == "crispr" | type == "rnai" | type == "drug"){
    # calculate normalized expression, if applicable
    if(!is.null(data_n)){
      if(input[[g_ui_norm_id]]=="gs"){
        # z score transform expression values
        n_col_n <- ncol(data_n)
        exp_scale_n <- apply(data_n[,2:n_col_n], 2, scale)
        data_n <- cbind(data_n[,1,drop=F],exp_scale_n) %>%
          mutate(Mean=rowMeans(exp_scale_n,na.rm=T)) %>%
          dplyr::select(patient_id, Mean)
        # data_n[["Mean"]] <- pnorm(data_n[["Mean"]])
      }
      data <- dplyr::left_join(data, data_n, by = "patient_id")
      data[,2] <- data[,2] / data[,3]
      data <- data[,1:2]
    }
    # save original expression data
    rv[[paste0("exprs_",x)]] <- data
  }else if(type == "snv"){
    muts <- data[,2] %>% unlist(.)
    names(muts) <- data$patient_id
    # muts <- muts[!is.na(muts)]
    rv[[paste0("mutations_",x)]] <- muts
  }else if(type == "lib" | type == "manual"){
    # z score transform expression values
    n_col <- ncol(data)
    exp_scale <- apply(data[,2:n_col], 2, scale)
    data <- cbind(data[,1,drop=F],exp_scale) %>%
      mutate(Mean=rowMeans(exp_scale,na.rm=T)) %>%
      dplyr::select(patient_id, Mean)

    # calculate normalized expression, if applicable
    if(!is.null(data_n)){
      if(input[[g_ui_norm_id]]=="gs"){
        # z score transform expression values
        n_col_n <- ncol(data_n)
        exp_scale_n <- apply(data_n[,2:n_col_n], 2, scale)
        data_n <- cbind(data_n[,1,drop=F],exp_scale_n) %>%
          mutate(Mean=rowMeans(exp_scale_n,na.rm=T)) %>%
          dplyr::select(patient_id, Mean)
        # data_n[["Mean"]] <- pnorm(data_n[["Mean"]])
      }
      data <- dplyr::left_join(data, data_n, by = "patient_id")
      data[,2] <- data[,2] / data[,3]
      data <- data[,1:2]
    }

    # save mean scaled FPKM data
    rv[[paste0("exprs_",x)]] <- data
  }

  return(data)
}

# original survival df
original_surv_df <- function(patient_ids){
  df_o <- rv$df_survival
  df_o %>% dplyr::filter(patient_id %in% patient_ids)
  df_o[match(patient_ids,df_o$patient_id),] %>%
    dplyr::select(-gender) %>%
    dplyr::filter(!is.na(patient_id))
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
surv_cox <- function(df, mode=1){
  # if(rv$depmap){
  #   if(mode == 1){
  #     coxph(Surv(dependency) ~ level, data = df)
  #   }else if(mode == 2){
  #     coxph(Surv(dependency) ~ level.x + level.y, data = df)
  #   }
  # }else{
    if(mode == 1){
      coxph(Surv(survival_days, censoring_status) ~ level, data = df)
    }else if(mode == 2){
      coxph(Surv(survival_days, censoring_status) ~ level.x * level.y, data = df)
    }
  # }
}
surv_km <- function(df){
  survfit(Surv(survival_days, censoring_status) ~ level, data = df)
}

transform_p_df <- function(p_df){
  p_df <- rbindlist(p_df)
  colnames(p_df) = c('p_value','quantile','expression','hr')
  p_df$p_value = as.numeric(p_df$p_value)
  p_df$quantile = as.numeric(p_df$quantile)
  p_df$expression = as.numeric(p_df$expression)
  p_df$hr = as.numeric(p_df$hr)
  return(p_df)
}

# function to check if expression-like
exp_yyy <- function(extract_mode){
  extract_mode != "snv" & ((!rv$depmap & extract_mode != "cnv") | rv$depmap)
}

# function to re-assign groups based on user's selection of a risk subgroup
assign_gp <- function(df,gp){
  lels_tmp <- unique(df$level)
  if(grepl("Gain",gp)){
    gp_hl <- strsplit(gp,"_")[[1]]; gp_hl_i <- which(gp_hl %in% c("High","Low"))
    gp_hl <- gp_hl[gp_hl_i]
    gain_or_loss <- strsplit(lels_tmp,"_") %>% unlist() %>% unique
    gain_or_loss <- gain_or_loss[!gain_or_loss %in% c("High","Low","Other")]
    if(gp_hl_i == 1){
      gain_or_loss_gp <- paste0(gp_hl,"_",gain_or_loss)
    }else if(gp_hl_i == 2){
      gain_or_loss_gp <- paste0(gain_or_loss,"_",gp_hl)
    }
  }else if(gp != "All"){
    gain_or_loss_gp <- gp
  }
  if(gain_or_loss_gp %in% lels_tmp){
    gain_or_loss_gp_other <- lels_tmp[lels_tmp != gain_or_loss_gp] %>% paste0(.,collapse = ", ")
    df$level <- ifelse(df$level == gain_or_loss_gp,gain_or_loss_gp,gain_or_loss_gp_other)
    df$level <- factor(df$level)
    df$level <- relevel(df$level, ref = gain_or_loss_gp_other)
  }else{
    # return NULL is user-selected group not found
    df <- NULL
  }
  return(df)
}

get_info_most_significant_rna <- function(data, min, max, step, num=1, data2=NULL, min2=NULL, max2=NULL, step2=NULL, gp=rv$risk_gp,cat=""){
  nCores <- detectCores() - 1
  # convert RVs into static variables
  depmap_T <- rv$depmap; p_kc <- rv$min_p_kc; gps <- rv$risk_gps
  print(gp)
  perc_min <- rv$min_gp_size / 100
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)

  # retrieve survival analysis df_o
  df_o <- original_surv_df(data$patient_id)
  if(num > 1){
    df_o2 <- original_surv_df(data2$patient_id)
    df_o <- dplyr::filter(df_o, patient_id %in% df_o2$patient_id)
    df_o2 <- dplyr::filter(df_o2, patient_id %in% df_o$patient_id)
    data <- data %>% dplyr::filter(patient_id %in% df_o$patient_id) %>% dplyr::filter(patient_id %in% df_o2$patient_id)
    data2 <- data2 %>% dplyr::filter(patient_id %in% df_o$patient_id) %>% dplyr::filter(patient_id %in% df_o2$patient_id)
  }else if(num == 1){
    # remove NA patient ids
    data <- data %>% dplyr::filter(patient_id %in% df_o$patient_id)
    colnames(data) <- c("patient_id","exp")
    if(cat != ""){
      # data: exp-like continuous data
      # data2: mut-like categorical data
      # combined into data: columns = patient_id, exp, mut
      colnames(data2) <- c("patient_id","mut")
      df_o2 <- original_surv_df(data2$patient_id)
      data2 <- data2 %>% dplyr::filter(patient_id %in% df_o$patient_id) %>% dplyr::filter(patient_id %in% df_o2$patient_id)
      data <- data %>% dplyr::inner_join(data2,by="patient_id")
      df_o <- dplyr::filter(df_o, patient_id %in% data$patient_id)
    }
  }

  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  exp <-data[,2] %>% unlist(.) %>% unname(.)
  # the quantiles we will use to define the level of gene percentages
  quantiles <- quantile(exp, quantile_s, na.rm = T)

  if(num > 1){
    quantile_s2 = seq(min2, max2, by = step2)
    patient_ids2 <- data2$patient_id
    exp2 <-data2[,2] %>% unlist(.) %>% unname(.)
    quantiles2 <- quantile(exp2, quantile_s2, na.rm = T)
    
    n_min_r <- perc_min * nrow(data2)

    rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
      q <- quantiles[i]
      df <- generate_surv_df(df_o, patient_ids, exp, q)

      rrr2 <- mclapply(seq_along(quantiles2),mc.cores = nCores,function(j){
        q2 <- quantiles2[j]
        # system(sprintf('echo "\n%s"', q2))
        df2 <- generate_surv_df(df_o2, patient_ids2, exp2, q2)
        df_list <- list(df,df2)
        # generate interaction df
        df_combined <- Reduce(
          function(x, y) inner_join(x, dplyr::select(y, patient_id, level), by = "patient_id"),
          df_list
        )
        x_y <- c("x","y")[1:length(df_list)]
        df_combined[["level"]] <- apply(df_combined %>% dplyr::select(paste0("level.",x_y)),1,paste0,collapse="_")
        # combine subgroups if indicated
        if(gp != "All"){
          other_gp <- paste0(gps[!gps %in% c(gp,"All")],collapse = ", ")
          df_combined[["level"]] <- ifelse(df_combined[["level"]] == gp,gp,other_gp)
          df_combined[["level"]] <- factor(df_combined[["level"]],levels = c(other_gp,gp))
        }

        # determine if meet min % samples requirement
        n_min <- min(table(df_combined[["level"]]))
        if(n_min < n_min_r){
          results <- NULL
        }else{
          # # test if there is significant difference between high and low level genes
          if(depmap_T){
            surv_diff <- kruskal.test(dependency ~ level, data = df)
            hr <- NA
            p_diff <- surv_diff$p.value #summary(surv_diff)[[1]][[5]][1]
          }else{
            if(p_kc == "km"){
              # surv_diff <- surv_km(df_combined)
              km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df_combined)
              p_diff <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
            }else if(p_kc == "cox"){
              if(gp == "All"){
                surv_diff <- surv_cox(df_combined,mode = 2)
              }else{
                surv_diff <- surv_cox(df_combined,mode = 1)
              }
              p_diff <- summary(surv_diff)$logtest[3] #coefficients[,5]
            }
          }

          if(!is.na(p_diff)){
            # #append current p value to the p value df
            new_row = c(p_diff,unlist(strsplit(names(quantiles2[j]),split = '%',fixed=T)),quantiles2[j],NA)
            results <- list(new_row,p_diff,df_combined,hr,names(quantiles2[j]))
            names(results) <- c("new_row","least_p_value","df_most_significant","least_hr","cutoff_most_significant")
          }else{
            results <- NULL
          }
        }
        return(results)
      })

      rrr2 <- Filter(Negate(is.null), rrr2)
      if(is.null(rrr2)){
        return(NULL)
      }else{
        # P tracking record on variable 2 j
        p_df2 <- lapply(rrr2, function(x) {data.frame(t(data.frame(x[["new_row"]])))})
        p_df2 <- try(transform_p_df(p_df2))
        # system(sprintf('echo "p_df2: %s hihi\n"', head(p_df2)))
        if(inherits(p_df2, "try-error")){
          return(NULL)
        }else{
          ## Determine the min-P point at percentile i in variable 1
          # Find the minimum P-value
          pvals <- sapply(rrr2, function(x) x[["least_p_value"]])
          pvals_i <- which.min(pvals)[[1]]
          res <- rrr2[[pvals_i]]
          least_p_value0 <- res[["least_p_value"]]
          df_most_significant0 <- res[["df_most_significant"]]
          cutoff_most_significant0 <- res[["cutoff_most_significant"]]
          
          new_row1 = c(least_p_value0,unlist(strsplit(names(quantiles[i]),split = '%',fixed=T)),quantiles[i],NA)
          
          # Proceed only if enough data
          if(is.null(df_most_significant0)){
            return(NULL)
          }else{
            results <- list(
              new_row = new_row1,
              least_p_value = least_p_value0,
              df_most_significant = df_most_significant0,
              least_hr = NA,
              cutoff_most_significant = c(names(quantiles[i]),cutoff_most_significant0)
              ,p_df = p_df2
            )
            return(results)
          }
        }
      }
    })
  }else{
    if(cat != ""){
      cat_s <- strsplit(cat,"_")[[1]]
      cat_si <- which(tolower(cat_s) %in% c("high","low"))
    }

    rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
      q <- quantiles[i]
      df <- generate_surv_df(df_o, patient_ids, exp, q)
      
      if(cat != ""){
        if(cat == "All"){
          df$level <- paste0(data[["mut"]],"_",df$level)
        }else if(cat_si == 1){
          df[["level.x"]] <- factor(df[["level"]]); df[["level.x"]] <- relevel(df[["level.x"]], ref = "Low")
          df[["level.y"]] <- factor(data[["mut"]]); df[["level.y"]] <- relevel(df[["level.y"]], ref = "Other")
          df$level <- paste0(df$level,"_",data[["mut"]])
          df <- assign_gp(df,cat)
        }else if(cat_si == 2){
          df[["level.x"]] <- factor(data[["mut"]]); df[["level.x"]] <- relevel(df[["level.x"]], ref = "Other")
          df[["level.y"]] <- factor(df[["level"]]); df[["level.y"]] <- relevel(df[["level.y"]], ref = "Low")
          df$level <- paste0(data[["mut"]],"_",df$level)
          df <- assign_gp(df,cat)
        }
      }
      
      if(is.null(df)){
        return(NULL)
      }else{
        # # test if there is significant difference between high and low level genes
        if(depmap_T){
          # surv_diff <- survdiff(Surv(dependency) ~ level, data = df)
          if(cat == "All"){
            surv_diff <- kruskal.test(dependency ~ level, data = df)
          }else{
            surv_diff <- wilcox.test(dependency ~ level, data = df)
          }
          hr <- NA
          p_diff <- surv_diff$p.value
        }else{
          if(cat==""){
            surv_diff <- surv_cox(df)
            hr <- coef(summary(surv_diff))[,2]
          }else if(cat=="All"){
            surv_diff <- surv_cox(df)
            hr <- NA
          }else{
            surv_diff <- surv_cox(df,mode=2)
            hr <- NA
          }
          if(p_kc == "km"){
            p_diff <- summary(surv_diff)$sctest[3]
          }else if(p_kc == "cox"){
            p_diff <- summary(surv_diff)$logtest[3] #coefficients[,5]
          }
        }
        
        if(!is.na(p_diff)){
          # #append current p value to the p value df
          new_row = c(p_diff,unlist(strsplit(names(quantiles[i]),split = '%',fixed=T)),quantiles[i],hr)
          # p_df <- rbind(p_df,new_row)
          # q_value <- as.numeric(sub("%", "", names(q)))
          # if((q_value <= 50 & p_diff <= least_p_value)|(q_value > 50 & p_diff < least_p_value)){
          #   least_p_value <- p_diff
          #   df_most_significant <- df
          #   least_hr <- hr
          #   cutoff_most_significant <- names(quantiles[i])
          # }
          results <- list(new_row,p_diff,df,hr,names(quantiles[i]))
          names(results) <- c("new_row","least_p_value","df_most_significant","least_hr","cutoff_most_significant")
          return(results)
        }
      }
    })
  }

  rrr <- Filter(Negate(is.null), rrr)
  if(length(rrr)>0){
    #Find the minimum P-value
    pvals <- sapply(rrr, function(x) x[["least_p_value"]])
    pvals_i <- which.min(pvals)[[1]]
    res <- rrr[[pvals_i]]
    df_most_significant <- res[["df_most_significant"]]
    # assign levels in order
    df_most_significant$level <- factor(df_most_significant$level)
    lels <- levels(df_most_significant$level)
    if(cat == "All"){
      lel_name <- lels[grepl("Low",lels)]
      lel_name <- lel_name[grepl("Other",lel_name)]
      df_most_significant$level <- relevel(df_most_significant$level, ref = lel_name)
    }
    cutoff_most_significant <- res[["cutoff_most_significant"]]
    least_hr <- res[["least_hr"]]
    
    #Transform the p_df a little bit to make it work with the ggplot
    if(num > 1){
      p_df1 <- lapply(rrr, function(x) {data.frame(t(data.frame(x[[1]])))})
      p_df1 <- transform_p_df(p_df1)
      p_df2 <- res[["p_df"]]
      p_df <- list(p_df1,p_df2)
    }else{
      p_df <- lapply(rrr, function(x) {data.frame(t(data.frame(x[[1]])))})
      p_df <- transform_p_df(p_df)
    }
    
    # proceed only if enough data
    if(is.null(df_most_significant)){
      return(NULL)
    }else{
      results <- list(
        df = df_most_significant,
        cutoff = cutoff_most_significant
        ,hr = least_hr
        ,p_df = p_df
        
      )
      return(results)
    }
  }else{
    return(NULL)
  }
}

# generate df if user-defined cutoffs
get_df_by_cutoff <- function(data, cutoff){
  cutoff <- cutoff/100

  # retrieve survival analysis df_o
  df_o <- original_surv_df(data$patient_id)
  # remove NA patient ids
  data <- data %>% dplyr::filter(patient_id %in% df_o$patient_id)
  # extract patients' IDs and expression values
  patient_ids <- data$patient_id
  exp <- data[,2] %>% unlist(.) %>% unname(.)

  # the quantile we will use to define the level of gene percentages
  q <- quantile(exp, cutoff)
  df <- generate_surv_df(df_o, patient_ids, exp, q)
}

# generate df if SNV mutation data
get_df_snv <- function(data, nons, syns){
  # retrieve survival analysis df_o
  df_o <- original_surv_df(data$patient_id)
  # remove NA patient ids
  data <- data %>% dplyr::filter(patient_id %in% df_o$patient_id)

  nons <- tolower(nons)
  syns <- tolower(syns)
  # extract patients' IDs and mutation data
  patient_ids <- data$patient_id
  mutations <- data[,2] %>% unlist(.) %>% unname(.)
  # check if nonsyn
  mutations[is.na(mutations)] <- "WT" #"Synonymous"
  mutations[mutations == ""] <- "WT" #"Synonymous"
  mm <- mutations #[!is.na(mutations)]
  mm <- sapply(mm, function(x){
    x <- strsplit(x, "\\|")[[1]]
    if(any(tolower(x) %in% nons)){
      "Mutated" #"Nonsynonymous"
    }else if(any(tolower(x) %in% syns)){
      "Other" #"Synonymous"
    }else{
      NA
    }
  })
  mutations <- mm #[!is.na(mutations)]
  data[,2] <- mutations

  # rename columns
  col_names <- colnames(data)
  col_names[length(col_names)] <- "level"
  colnames(data) <- col_names

  # filter data
  data <- data %>% dplyr::filter(!is.na(level))
  # reset levels
  lels <- unique(data$level) %>% sort(.,decreasing = T)
  data$level <- factor(data$level, levels = lels)

  df <- df_o %>% inner_join(data, by = "patient_id")
}

# generate df if CNV copy number data
get_info_most_significant_cnv <- function(data, mode){
  # initialize the most significant p value
  least_p_value <- 1

  # retrieve survival analysis df_o
  df_o <- original_surv_df(data$patient_id)
  # remove NA patient ids
  data <- data %>% dplyr::filter(patient_id %in% df_o$patient_id)

  # extract patients' IDs and expression values
  patient_ids <- data$patient_id

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
      surv_diff <- surv_cox(df)
      p_diff <- summary(surv_diff)$logtest[3] #coef(summary(surv_diff))[,5]
      # }else if(rv$cox_km == "km"){
      #   surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
      #   p_diff <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
      # }
      if(!is.na(p_diff)){
        if(p_diff <= least_p_value){
          least_p_value <- p_diff
          df_most_significant <- df
          cutoff_most_significant <- names(cats[cat])
        }
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

## multiple test p value correction
correct_p <- function(p_diff,min,max,step){
  p_z <- qnorm(1 - p_diff/2) # (1-Pmin/2)-quantile of the standard normal distribution
  # p_dens <- dnorm(p_z) # probability density funciton
  if(is.infinite(p_z)){
    p_z <- 10
  }
  # if(is.infinite(p_z)){
  #   p_diff_adj <- p_diff
  # }else{
  #   p_diff_adj <- p_dens * (p_z - (1/p_z)) * log((max*(1-min))/((1-max)*min)) + 4 * p_dens/p_z
  # }
  p_zz <- p_z ^ 2
  p_acc_2 <- 0
  p_acc_3 <- 0
  quantiles <- seq(min,max,step)
  for(i in 1:(length(quantiles)-1)){
    qqs <- 1 - (quantiles[i]*(1-quantiles[i+1]))/((1-quantiles[i])*quantiles[i+1])
    qqs <- sqrt(qqs)
    p_acc_2 <- p_acc_2 + qqs
    p_acc_3 <- p_acc_3 + qqs^3
    # p_acc <- p_acc + sqrt(qqs) - ((p_zz / 4 - 1) * (sqrt(qqs))^3 / 6)
  }
  p_diff_adj <- exp(-(p_zz)/2) / pi * (p_acc_2 - (p_zz / 4 - 1) * p_acc_3 / 6)
  return(p_diff_adj)
}

## Perform survival analysis
cal_surv_rna <-
  function(
    df,n,min,max,step
    ,p.adjust.method = rv[["km_mul"]]
    ,iter_mode=T
    ,gp=rv$risk_gp
  ){
    lels <- levels(df$level)
    # # 1. KM # #
    # summary statistics
    if(rv$depmap){
      # km.stats <- survdiff(Surv(dependency) ~ level, data = df)
      if(n == 1){
        surv_diff <- wilcox.test(dependency ~ level, data = df)
        p_diff <- surv_diff$p.value
      }else if(n == 2){
        if(gp=="All"){
          surv_diff <- kruskal.test(dependency ~ level, data = df)
        }else{
          surv_diff <- wilcox.test(dependency ~ level, data = df)
        }
        p_diff <- surv_diff$p.value#summary(surv_diff)[[1]][[5]][1]
        # surv_diff <- pairwise.t.test(df$dependency, df$level, p.adjust.method = "hommel")
      }

      if(iter_mode){
        p_diff_adj <- correct_p(p_diff,min,max,step)
      }else{
        p_diff_adj <- NULL
      }

      results = list(
        df = df
        ,fit = surv_diff
        # ,stats = km.stats
        ,lels = lels
        # ,hr = "NA"
        ,p = p_diff
        ,p.adj = p_diff_adj
      )
    }else{
      km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df)
      p.km <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)

      # run KM
      km.fit <- surv_km(df)

      # # 2. Cox # #
      if(n == 1){
        # run Cox regression
        cox_fit <- surv_cox(df)
        # summary statistics
        cox.stats <- summary(cox_fit)
      }else if(n == 2){
        # # if(rv$depmap){
        # #   km2 <- try(pairwise_survdiff(Surv(dependency) ~ level, data = df, p.adjust.method = p.adjust.method))
        # # }else{
        #   km2 <- try(pairwise_survdiff(Surv(survival_days, censoring_status) ~ level, data = df, p.adjust.method = p.adjust.method))
        # # }
        # if(inherits(km2, "try-error")) {
        #   rv$try_error <- rv$try_error + 1
        #   shinyalert(paste0("The selected two genes/loci have exactly the same data."
        #                     ," If CNV, you might have selected genes from the same cytoband."
        #                     ," Please contrast genes from different cytobands."))
        # }
        # req(!inherits(km2, "try-error")) #require to be no error to proceed the following codes
        #
        # # if(iter_mode){
        # #   km2$p.value <- paste0(res.km$p.value," (",correct_p(km2$p.value,min,max,step),")")
        # # }
        # km.stats <- list(km2,km.stats)

        # run Cox regression
        if(gp=="All"){
          cox_fit1 <- surv_cox(df, mode=2)
        }else{
          cox_fit1 <- surv_cox(df, mode=1)
        }
        # summary statistics
        cox.stats <- summary(cox_fit1)
        # run Cox regression for visualization purpose
        cox_fit <- surv_cox(df)
      }

      hr.cox <- sapply(cox.stats$coefficients[,2], function(x){
        round(as.numeric(x), 2)
      }) %>% paste0(.,collapse = ", ")
      p.cox <- as.numeric(cox.stats$logtest[3])
      #   sapply(cox.stats$coefficients[,5], function(x){
      #   as.numeric(x)
      # })

      # multiple p correction
      if(iter_mode){
        p.km.adj <- correct_p(p.km,min,max,step)
        p.cox.adj <- correct_p(cox.stats$logtest[3],min,max,step)
        #   sapply(cox.stats$coefficients[,5], function(x){
        #   correct_p(as.numeric(x),min,max,step)
        # })
      }else{
        p.km.adj <- NULL
        p.cox.adj <- NULL
      }


      # # run Cox survival analysis
      # lels_x <- levels(df$`level.x`)
      # lels_y <- levels(df$`level.y`)
      # # lels <- apply(expand.grid(lels_x,lels_y),1,paste0,collapse="_")
      # create new df to seperate effects in Cox regression
      new_df <- with(df,data.frame(level = lels))
      cox.fit <- survfit(cox_fit,newdata=new_df)

      # save df, fit, and statistics
      if(!exists("cox_fit1")){cox_fit1 <- cox_fit}
      results <- list(
        km = list(
          df = df,
          fit = km.fit,
          stats = km.stats
          ,lels = lels
          ,hr = "NA"
          ,p = p.km
          ,p.adj = p.km.adj
        )
        ,cox = list(
          df = new_df,
          fit = cox.fit,
          stats = cox.stats
          ,lels = lels
          ,hr = hr.cox
          ,p = p.cox
          ,p.adj = p.cox.adj
          ,cox_fit = cox_fit1
          ,cox_df = df
        )
      )
    }

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
    , palette=rv$palette
    , base_size=20
  ){
    df <- res[[mode]][["df"]]
    # extract data
    fit <- res[[mode]][["fit"]]
    lels <- res[[mode]][["lels"]]
    req(!is.null(fit))
    # adjust colors, if applicable
    if(palette == "br"){
      n_lel <- length(lels)
      col_alt <- c("#939597","#F0A1BF") # grey pink
      if(n_lel > 2){
        palette <- c("black",col_alt[1:(n_lel-2)],"red")
      }else{
        palette <- c("black","red")
      }
    }
    # median survival lines
    if(is.null(rv$median)){
      surv.median.line="none"
    }else if(length(rv$median)==2){
      surv.median.line="hv"
    }else{
      surv.median.line=rv$median
    }

    if(!rv$depmap){
      # x-axis time intervals
      if(rv$ymd == "d"){
        xscale <- 1; xlab <- "Days"; breaktime <- rv$ymd_int_d
      }else if(rv$ymd == "m"){
        xscale <- "d_m"; xlab <- "Months"; breaktime <- 30.4375 * rv$ymd_int_m # 608.75 # 30.4375 * 20 = 365.25
      }else if(rv$ymd == "y"){
        xscale <- "d_y"; xlab <- "Years"; breaktime <- 365.25 * rv$ymd_int_y
      }
    }else{
      dep_name <- dependency_names() %>% firstlower()
      xscale <- 1; xlab <- paste0("Exponential of ",dep_name); breaktime <- NULL
    }

    if(mode == "km"){
      fig <- ggsurvplot(fit, data=df,
                        title = title,
                        xlab = xlab,
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
                        ,xscale = xscale
                        ,break.time.by = breaktime
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
                        xlab = xlab,
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
                        ,xscale = xscale
                        ,break.time.by = breaktime
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
####            DepMap functions            ####
#==============================================#
translate_cells <- function(patient_ids){
  df <- rv$depmap_ccle
  cells <- tibble(patient_id = patient_ids) %>% dplyr::inner_join(df, by = "patient_id") %>% .[["CCLE_Name"]]
}

retrieve_dens_df <- function(){
  df <- rv[["res"]][["df"]]
  req(!is.null(df[["dependency"]]))
  # if(rv$project == "DepMap-Drug"){
  #   df[["dependency"]] <- log2(df[["dependency"]])
  # }else{
  #   df[["dependency"]] <- log10(df[["dependency"]])
  # }
  colnames(df)[2] <- dependency_names()
  colnames(df)[ncol(df)] <- "Level"
  df[["Cell"]] <- paste0(translate_cells(df$patient_id),"|",df$patient_id)
  # df[["Cell"]] <- paste0(gsub("(^.*?)_(.*)","\\1",translate_cells(df$patient_id)),"|",df$patient_id)
  df <- df %>% dplyr::select(-patient_id)
  rv[["dens_df"]] <- df
  rv$annot_cells_y <- "yes"
  return(df)
}
#==============================================#
####            eVITTA functions            ####
#==============================================#
# the gene expression table
de_dfgene <- function(){
  # patients
  if(rv$variable_nr == 1){
    patients <- rv[["df_gender"]][["patient_id"]]
  }else{
    patients <- rv[["df_all"]][["patient_id"]]
  }

  # read in data
  infiles <- paste0(rv$indir,"df_gene.csv")

  l <- lapply(infiles, function(x){
    info <- fread(x,sep=",",header=T)
    return(info)
  })

  # whole gene expression table
  df_gene <- rbind_common(l) %>% dplyr::filter(patient_id %in% patients)

  # filter ? genes
  if(rv$depmap | rv$tcga){
    df_gene <- df_gene %>% dplyr::select(-starts_with("\\?"))
  }

  # remaining patient ids
  patients <- df_gene$patient_id
  # the genes
  genes <- colnames(df_gene)[-1]

  incProgress(amount = 0.1, message = wait_msg("Formatting gene expression matrix for DE analysis..."))

  # transpose
  df_gene <- data.table::transpose(df_gene) %>% .[-1,] %>% as.data.frame()
  # reassign patient ids
  colnames(df_gene) <- patients
  rownames(df_gene) <- genes

  incProgress(amount = 0.2, message = wait_msg("Converting gene IDs..."))

  if(rv$target){
    # id conversion
    egSYMBOL <- toTable(org.Hs.egSYMBOL)

    # create individual tables using org.Hs
    egENS <- toTable(org.Hs.egENSEMBL)

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

    df_gene <- df_gene %>%
      dplyr::select(-c(ensembl_id, gene_id, symbol))
  }else{
    genes <- sapply(genes, function(x) str_split(x, "\\|")[[1]][1])
  }

  # convert to numeric matrix
  df_gene <- df_gene %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix(.)

  rownames(df_gene) <- genes

  return(df_gene)
}

#==============================================#
####            General functions            ####
#==============================================#
data_summary <- function(x,k=rv$violin_k) {
  m <- mean(x)
  ymin <- m - k * sd(x)
  ymax <- m + k * sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
