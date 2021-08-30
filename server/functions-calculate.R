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
    ,"manual" = "df_gene.csv"
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
    data <- dplyr::arrange(data, patient_id)
    if(!is.null(l_n)){
      data_n <- rbindlist(l_n, use.names = T)
      data_n <- dplyr::arrange(data_n, patient_id)
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
  if(type == "rna" | type == "mir" | type == "met" | type == "rrpa" | type == "crispr" | type == "rnai" | type == "drug" | (type == "cnv" & rv$depmap)){
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
  df_o <- df_o[match(patient_ids,df_o$patient_id),] %>%
    dplyr::select(-gender) %>%
    dplyr::filter(!is.na(patient_id))
  arrange(df_o, patient_id)
}

# generate survival df
generate_surv_df <- function(df, patient_ids, exp, q){
  # generate the data from for model first
  gene_quantiles <- ifelse(exp > q, "High", "Low")
  # names(gene_quantiles) <- patient_ids

  # # generate survival analysis df
  # df$level <- gene_quantiles[match(df$patient_id,names(gene_quantiles))]
  df$level <- gene_quantiles
  lels <- unique(df$level) %>% sort(.,decreasing = T)
  df$level <- factor(df$level, levels = lels)

  # df <- df %>% dplyr::filter(!is.na(patient_id))
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
# function to check if minimum P-value searching
exp_iter_yyy <- function(x){
  iter_id <- paste0("iter_",x)
  ifelse(is.null(input[[iter_id]]), T, input[[iter_id]] == "iter")
}
# function to determine if any iteration to perform
if_any_iter <- function(n=rv$variable_n){
  T
  # all(sapply(1:n,function(x) exp_yyy(input_mode(x)) & exp_iter_yyy(x)))
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
#Helper function to populate a dataframe made from quantiles with p value column default set to 1
populate_quantile_df <- function(quantile_s, min, max, step, substitue_value = 1){
  temp <- lapply(quantile_s, function(i){
    rep(i, length(quantile_s))
  })
  temp <- unlist(temp)
  temp_second_col <- rep(quantile_s, length(quantile_s))
  result <- data.frame(temp, temp_second_col, substitue_value)
  colnames(result) <- c("Q1","Q2", "p_value")
  #convert as numeric dataframe
  result <- data.frame(sapply(result, function(x) as.numeric(as.character(x))))
  result$Q1 = result$Q1 * 100
  result$Q2 = result$Q2 * 100
  return(result)
}
#Helper function to update a heatmap dataframe's columns based on a index mapping df
update_heatmap_column <- function(target_df,index_df,column_name = "p_value"){
  result <- target_df
  result[[column_name]][index_df$index] <- as.numeric(index_df[[column_name]])
  result <- result[[column_name]]
  return(result)
}
#takes in rrr(a list of multiple attributes including many dfs) and produce a dataframe with assembled df for column "heatmap_new_row"
assemble_new_rows <- function(rrr){
  heatmap_new_rows <- lapply(rrr, function(x) {data.frame(t(data.frame(x[["heatmap_new_row"]])))})
  heatmap_new_rows <- do.call("rbind", heatmap_new_rows)
  rownames(heatmap_new_rows) <- NULL
  return(heatmap_new_rows)
}
#create a tibble based on heatmap df
create_tracker <- function(heatmap_df, val = NA){
  res <- heatmap_df[,1:2]
  res$p_value <- val
  res$hr <- val
  return(res)
}



# the function to assign groups in ONE-GENE-like when a risk group is selected in 2-gene analysis
assign_df_levels <- function(df, data, cat, cat_si){
  if(cat != ""){
    if(cat == "All"){
      df[["level.x"]] <- factor(df[["level"]]); df[["level.x"]] <- relevel(df[["level.x"]], ref = "Low")
      df[["level.y"]] <- factor(data[["mut"]]); df[["level.y"]] <- relevel(df[["level.y"]], ref = sort(levels(df[["level.y"]]),decreasing = T)[1])
      df$level <- paste0(data[["mut"]],"_",df$level)
    }else if(cat_si == 1){
      df[["level.x"]] <- factor(df[["level"]]); df[["level.x"]] <- relevel(df[["level.x"]], ref = "Low")
      df[["level.y"]] <- factor(data[["mut"]]); df[["level.y"]] <- relevel(df[["level.y"]], ref = sort(levels(df[["level.y"]]),decreasing = T)[1])
      df$level <- paste0(df$level,"_",data[["mut"]])
      df <- assign_gp(df,cat)
    }else if(cat_si == 2){
      df[["level.x"]] <- factor(data[["mut"]]); df[["level.x"]] <- relevel(df[["level.x"]], ref = sort(levels(df[["level.x"]]),decreasing = T)[1])
      df[["level.y"]] <- factor(df[["level"]]); df[["level.y"]] <- relevel(df[["level.y"]], ref = "Low")
      df$level <- paste0(data[["mut"]],"_",df$level)
      df <- assign_gp(df,cat)
    }
  }
  return(df)
}

# function to find the minimum P value in a list
find_minP_res <- function(rrr2){
  pvals <- sapply(rrr2, function(x) x[["least_p_value"]])
  pvals_i <- which.min(pvals)[[1]]
  rrr2[[pvals_i]]
}

# function to fit survival curves onto ONE GENE
one_gene_cox <- function(df,cat,q,depmap_T,p_kc,new_row_T=T){
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
      if(cat=="All"){
        surv_diff <- surv_cox(df,mode=2)
        hr <- NA
      }else{
        surv_diff <- surv_cox(df)
        hr <- coef(summary(surv_diff))[,2]
      }
      if(p_kc == "km"){
        p_diff <- summary(surv_diff)$sctest[3]
      }else if(p_kc == "cox"){
        p_diff <- summary(surv_diff)$logtest[3] #coefficients[,5]
      }
    }

    if(!is.na(p_diff)){
      # p_df <- rbind(p_df,new_row)
      # q_value <- as.numeric(sub("%", "", names(q)))
      # if((q_value <= 50 & p_diff <= least_p_value)|(q_value > 50 & p_diff < least_p_value)){
      #   least_p_value <- p_diff
      #   df_most_significant <- df
      #   least_hr <- hr
      #   cutoff_most_significant <- names(quantiles[i])
      # }
      if(new_row_T){
        # #append current p value to the p value df
        new_row = c(p_diff,gsub("%$","",names(q)),q,hr)
        results <- list(new_row,p_diff,df,hr,names(q))
        names(results) <- c("new_row","least_p_value","df_most_significant","least_hr","cutoff_most_significant")
      }else{
        results <- list(p_diff,df,hr,names(q))
        names(results) <- c("least_p_value","df_most_significant","least_hr","cutoff_most_significant")
      }
      return(results)
    }
  }
}

# function to loop quantiles2 based on q in quantiles
two_gene_cox_inner <- function(
  q, q2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, heatmap_df
){
  # system(sprintf('echo "\n%s"', q2)
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
    df_combined[["level"]] <- ifelse(df_combined[["level"]] == gp,gp,other_gp)
    df_combined[["level"]] <- factor(df_combined[["level"]],levels = c(other_gp,gp))
  }

  # determine if meet min % samples requirement
  n_mins <- table(df[["level"]]); n_min <- min(n_mins)
  return_null <- F
  # if(gp == "All"){if(length(n_mins) < 4){return_null <- T}}
  if(n_min < n_min_r){return_null <- T}
  if(return_null){
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
        hr <- NA
      }else if(p_kc == "cox"){
        if(gp == "All"){
          surv_diff <- surv_cox(df_combined,mode = 2)
          hr <- NA
        }else{
          surv_diff <- surv_cox(df_combined,mode = 1)
          hr <- coef(summary(surv_diff))[,2]
        }
        p_diff <- summary(surv_diff)$logtest[3] #coefficients[,5]
      }
    }

    if(!is.na(p_diff)){
      # # #append current p value to the p value df
      new_row = c(p_diff,gsub("%$","",names(q2)),q2,NA)
      #TODO: optimize!
      quantile_1  = unlist(strsplit(names(q),split = '%',fixed=T))
      quantile_2  = unlist(strsplit(names(q2),split = '%',fixed=T))
      #start of heatmp section:
      heatmap_df_index = which((heatmap_df$Q1 == as.numeric(quantile_1))&(heatmap_df$Q2 == as.numeric(quantile_2)))
      #print(hr)
      #heatmap_new_row = c(round(heatmap_df_index,0),p_diff)#,hr)
      #names(heatmap_new_row) <- c("index","p_value")#,"hr")
      heatmap_new_row = c(round(heatmap_df_index,0),p_diff,hr)
      names(heatmap_new_row) <- c("index","p_value","hr")
      results <- list(new_row,p_diff,df_combined,hr,names(q2),heatmap_new_row)
      names(results) <- c("new_row","least_p_value","df_most_significant","least_hr","cutoff_most_significant","heatmap_new_row")
    }else{
      results <- NULL
    }
  }
  return(results)
}

two_gene_cox <- function(
  q, quantiles2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores, new_row_T=T,heatmap_df
){
  #TODO CHANGE BACK TO mclapply
  #rrr2 <- mclapply(seq_along(quantiles2),mc.cores = nCores,function(j){
  rrr2 <- lapply(seq_along(quantiles2),function(j){
    q2 <- quantiles2[j]
    two_gene_cox_inner(q, q2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T,heatmap_df = heatmap_df)
  })

  rrr2 <- Filter(Negate(is.null), rrr2)
  if(length(rrr2)<1){
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
      res <- find_minP_res(rrr2)
      least_p_value0 <- res[["least_p_value"]]
      df_most_significant0 <- res[["df_most_significant"]]
      cutoff_most_significant0 <- res[["cutoff_most_significant"]]
      heatmap_new_rows <- assemble_new_rows(rrr2)
      #   lapply(rrr2, function(x) {
      #   data.frame(t(data.frame(x[["heatmap_new_row"]])))})
      # heatmap_new_rows <- do.call("rbind", heatmap_new_rows)

      # Proceed only if enough data
      if(is.null(df_most_significant0)){
        return(NULL)
      }else{
        if(new_row_T){
          new_row1 = c(least_p_value0,gsub("%$","",names(q)),q,NA)
          results <- list(
            new_row = new_row1,
            least_p_value = least_p_value0,
            df_most_significant = df_most_significant0,
            least_hr = NA,
            cutoff_most_significant = c(names(q),cutoff_most_significant0)
            ,p_df = p_df2
            ,heatmap_new_rows = heatmap_new_rows
          )
        }else{
          results <- list(
            least_p_value = least_p_value0,
            df_most_significant = df_most_significant0,
            least_hr = NA,
            cutoff_most_significant = c(names(q),cutoff_most_significant0)
            ,p_df = p_df2
            #not add new lines
            #,heatmap_new_rows = heatmap_new_rows
          )
        }
        return(results)
      }
    }
  }
}

# heuristic search with two genes
two_gene_heuristic <- function(
  quantiles, quantiles2,
  df_o, patient_ids, exp,
  df_o2, patient_ids2, exp2
  ,gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores, heatmap_df, include_tracking = FALSE
){

  # create a tracking dataframe
  i_len <- length(quantiles); j_len <- length(quantiles2)
  df_tracking <- matrix(NA, i_len, j_len)
  colnames(df_tracking) <- 1:i_len#names(quantiles2)
  rownames(df_tracking) <- 1:j_len#names(quantiles)
  update_tracker <- create_tracker(heatmap_df)

  # start from the median quantile
  i <- floor(i_len/2); q <- quantiles[i]
  # if not enough data, skip and render users an error msg
  if(q == 0 | length(unique(quantiles2))==1){return(NULL);next;}
  # proceed only if enough data
  df <- generate_surv_df(df_o, patient_ids, exp, q)

  # mark initial iterated points
  df_tracking[,i] <- 1

  # loop quantiles2 using q
  rrr2 <- two_gene_cox(q, quantiles2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores,heatmap_df = heatmap_df)
  # anchor the optimized start point of heuristic searching
  q2 <- rrr2[["cutoff_most_significant"]][2]
  # if not enough data, skip and render users an error msg
  if(is.null(q2)){return(NULL);next;}
  j <- which(names(quantiles2) == q2)
  names(j) <- q2; q2 <- quantiles2[j]
  rrr2[["cutoff_most_significant"]][1] <- names(q)
  init_min_p <- rrr2[["least_p_value"]]
  View(rrr2)
  # start surrounding searching
  final_min_p <- init_min_p
  while(i > 1 & i < i_len & j > 1 & j < j_len){
    names(i) <- names(quantiles[i])
    names(j) <- names(quantiles2[j])
    a <- c(i-1,j); b <- c(i,j-1); c <- c(i+1,j); d <- c(i,j+1)
    comb <- list(a,b,c,d)
    if(include_tracking){
      print(comb)
    }
    # fix regression models via parallel processing
    #rrr_sr <- mclapply(1:4,mc.cores = nCores,function(k){
    rrr_sr <- lapply(1:4,function(k){
      ij_k <- comb[[k]]; i_k <- ij_k[1]; j_k <- ij_k[2]
      # skip if tracked
      if(!is.na(df_tracking[j_k,i_k])){
        return(NULL)
      }else{
        # fit cox regression
        q <- quantiles[i_k]; q2 <- quantiles2[j_k]
        df <- generate_surv_df(df_o, patient_ids, exp, q)
        results <- two_gene_cox_inner(q, q2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, heatmap_df = heatmap_df)
        if(!is.null(results)){
          results[["cutoff_most_significant"]] <- names(c(q,q2))
          names(ij_k) <- names(c(q,q2))
          ij_k <- list(ij_k); names(ij_k) <- "ij"
          results <- append(results,ij_k)#,heatmap_info)
        }
        return(results)
      }
    })

    if(include_tracking){
      View(df_tracking)
    }
    # the new P, if any
    rrr_sr <- Filter(Negate(is.null), rrr_sr)

    # size of rrr_sr may vary so it is best to update values here
    # assign 1 to tracked quantile combination
    View(rrr_sr)
    #create a cutoffs df in order to replace cell values in update_tracker

    heatmap_new_rows <- lapply(rrr_sr, function(x) {data.frame(t(data.frame(x[["heatmap_new_row"]])))})
    heatmap_new_rows <- do.call("rbind", heatmap_new_rows)

    #heatmap_new_rows <- do.call("rbind",heatmap_new_rows)
    print(heatmap_new_rows)

    ijs <- unlist(comb)
    df_tracking[ijs[2],ijs[1]] <- 1; df_tracking[ijs[4],ijs[3]] <- 1; df_tracking[ijs[6],ijs[5]] <- 1; df_tracking[ijs[8],ijs[7]] <- 1
    # the new P, if any
    rrr_sr <- Filter(Negate(is.null), rrr_sr)
    if(length(rrr_sr)>0){

      #update heatmap info
      res <- find_minP_res(rrr_sr)
      p_min <- res[["least_p_value"]]
      ij <- res[["ij"]]; i <- ij[1]; j <- ij[2]
      # mark new min P
      if(p_min < final_min_p){
        final_min_p <- p_min
        rrr <- rrr_sr
      }else{
        break
      }
    }else{
      break
    }
  }
  #End of while loop
  # assign to rrr
  if(final_min_p >= init_min_p){
    rrr <- list()
    rrr[[1]] <- rrr2
  }
  #if need to include tracking df
  if(include_tracking){
    return(list(rrr,update_tracker))
  }else{
    return(rrr)
  }

}




get_info_most_significant_rna <-
  function(
    data, min, max, step,
    num=1, data2=NULL, min2=NULL, max2=NULL, step2=NULL,
    gp=rv$risk_gp,cat="",
    search_mode=rv$search_mode, n_perm=rv$n_perm
  ){
  nCores <- detectCores() - 2; if(nCores > 12) nCores <- 12
  # convert RVs into static variables
  depmap_T <- rv$depmap; p_kc <- rv$min_p_kc; gps <- rv$risk_gps
  if(gp != "All"){
    other_gp <- paste0(gps[!gps %in% c(gp,"All")],collapse = ", ")
  }else{
    other_gp <- NULL
  }
  print(gp)
  perc_min <- rv$min_gp_size / 100
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)
  #initiate heatmap dataframe according to quantile_s
  heatmap_df <- populate_quantile_df(quantile_s,min, max, step, 1)
  heatmap_df <- data.frame(sapply(heatmap_df, function(x) as.numeric(as.character(x))))
  heatmap_df$annotation <- NA
  heatmap_df$hr <- NA

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
  # permutation matrix
  set.seed(0)
  n_samples <- nrow(df_o)
  idx.mat <- matrix(NA, n_samples, n_perm)
  for(ii in 1:n_perm) idx.mat[,ii] <- sample(1:n_samples)

  # ---- TWO GENES ----
  if(num > 1){
    quantile_s2 = seq(min2, max2, by = step2)
    patient_ids2 <- data2$patient_id
    exp2 <-data2[,2] %>% unlist(.) %>% unname(.)
    quantiles2 <- quantile(exp2, quantile_s2, na.rm = T)

    n_min_r <- perc_min * nrow(data2)
    # rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
    #   q <- quantiles[i]
    #   df <- generate_surv_df(df_o, patient_ids, exp, q)
    #
    #   rrr2 <- mclapply(seq_along(quantiles2),mc.cores = nCores,function(j){
    #     q2 <- quantiles2[j]
    #     # system(sprintf('echo "\n%s"', q2))
    #     df2 <- generate_surv_df(df_o2, patient_ids2, exp2, q2)
    #     df_list <- list(df,df2)
    #     # generate interaction df
    #     df_combined <- Reduce(
    #       function(x, y) inner_join(x, dplyr::select(y, patient_id, level), by = "patient_id"),
    #       df_list
    #     )
    #     x_y <- c("x","y")[1:length(df_list)]
    #     df_combined[["level"]] <- apply(df_combined %>% dplyr::select(paste0("level.",x_y)),1,paste0,collapse="_")
    #     # combine subgroups if indicated
    #     if(gp != "All"){
    #       other_gp <- paste0(gps[!gps %in% c(gp,"All")],collapse = ", ")
    #       df_combined[["level"]] <- ifelse(df_combined[["level"]] == gp,gp,other_gp)
    #       df_combined[["level"]] <- factor(df_combined[["level"]],levels = c(other_gp,gp))
    #     }
    #
    #     # determine if meet min % samples requirement
    #     n_min <- min(table(df_combined[["level"]]))
    #     if(n_min < n_min_r){
    #       results <- NULL
    #     }else{
    #       # # test if there is significant difference between high and low level genes
    #       if(depmap_T){
    #         surv_diff <- kruskal.test(dependency ~ level, data = df)
    #         hr <- NA
    #         p_diff <- surv_diff$p.value #summary(surv_diff)[[1]][[5]][1]
    #       }else{
    #         if(p_kc == "km"){
    #           # surv_diff <- surv_km(df_combined)
    #           km.stats <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df_combined)
    #           p_diff <- 1 - pchisq(km.stats$chisq, length(km.stats$n) - 1)
    #         }else if(p_kc == "cox"){
    #           if(gp == "All"){
    #             surv_diff <- surv_cox(df_combined,mode = 2)
    #           }else{
    #             surv_diff <- surv_cox(df_combined,mode = 1)
    #           }
    #           p_diff <- summary(surv_diff)$logtest[3] #coefficients[,5]
    #         }
    #       }
    #
    #       if(!is.na(p_diff)){
    #         # #append current p value to the p value df
    #         quantile_1  = unlist(strsplit(names(quantiles[i]),split = '%',fixed=T))
    #         quantile_2  = unlist(strsplit(names(quantiles2[j]),split = '%',fixed=T))
    #
    #         new_row = c(p_diff,quantile_2,quantiles2[j],NA)
    #         #start of heatmp section:
    #         heatmap_df_index = which((heatmap_df$Q1 == as.numeric(quantile_1))&(heatmap_df$Q2 == as.numeric(quantile_2)))
    #         heatmap_new_row = c(round(heatmap_df_index,0),p_diff)
    #         names(heatmap_new_row) <- c("index","p_value")
    #         #end of heatmap
    #
    #         results <- list(new_row,p_diff,df_combined,hr,names(quantiles2[j]),heatmap_new_row)
    #         names(results) <- c("new_row","least_p_value","df_most_significant","least_hr","cutoff_most_significant","heatmap_new_row")
    #       }else{
    #         results <- NULL
    #       }
    #     }
    #
    #     return(results)
    #   })
    #   #View(rrr2)
    #   #print(rrr2)
    #   #print("CALCULATED RRR2")
    #   rrr2 <- Filter(Negate(is.null), rrr2)
    #   if(is.null(rrr2)){
    #     return(NULL)
    #   }else{
    #
    #     # P tracking record on variable 2 j
    #     p_df2 <- lapply(rrr2, function(x) {data.frame(t(data.frame(x[["new_row"]])))})
    #     p_df2 <- try(transform_p_df(p_df2))
    #     # system(sprintf('echo "p_df2: %s hihi\n"', head(p_df2)))
    #     if(inherits(p_df2, "try-error")){
    #       return(NULL)
    #     }else{
    #       ## Determine the min-P point at percentile i in variable 1
    #       # Find the minimum P-value
    #       pvals <- sapply(rrr2, function(x) x[["least_p_value"]])
    #       pvals_i <- which.min(pvals)[[1]]
    #       res <- rrr2[[pvals_i]]
    #       least_p_value0 <- res[["least_p_value"]]
    #       df_most_significant0 <- res[["df_most_significant"]]
    #       cutoff_most_significant0 <- res[["cutoff_most_significant"]]
    #       heatmap_new_rows <- lapply(rrr2, function(x) {data.frame(t(data.frame(x[["heatmap_new_row"]])))})
    #       heatmap_new_rows <- do.call("rbind", heatmap_new_rows)
    #       new_row1 = c(least_p_value0,unlist(strsplit(names(quantiles[i]),split = '%',fixed=T)),quantiles[i],NA)
    #
    #       # Feel free to run this command to see what is inside of heatmap_new_rows variable
    #       # View(heatmap_new_rows)
    #
    #       # Proceed only if enough data
    #       if(is.null(df_most_significant0)){
    #         return(NULL)
    #       }else{
    #         results <- list(
    #           new_row = new_row1,
    #           least_p_value = least_p_value0,
    #           df_most_significant = df_most_significant0,
    #           least_hr = NA,
    #           cutoff_most_significant = c(names(quantiles[i]),cutoff_most_significant0)
    #           ,p_df = p_df2
    #           ,heatmap_new_rows = heatmap_new_rows
    #         )
    #         return(results)
    #       }
    #     }
    #   }
    # })

    #THIS IS THE NEW METHOD:
    # ---- TWO GENES exhaustive ----

    if(search_mode == "exhaustive"){
      #TODO: CHANGE BACK
      rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
      #rrr <- lapply(seq_along(quantiles),function(i){
        q <- quantiles[i]
        df <- generate_surv_df(df_o, patient_ids, exp, q)

        two_gene_cox(q, quantiles2, df_o2, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores, heatmap_df = heatmap_df)
      })

    # ---- TWO GENES heuristic ----
    }else if(search_mode == "heuristic"){
      rrr <- two_gene_heuristic(quantiles, quantiles2,
                         df_o, patient_ids, exp,
                         df_o2, patient_ids2, exp2
                         ,gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores, heatmap_df = heatmap_df, include_tracking = TRUE)
      df_tracking <- rrr[[2]]
      rrr <- rrr[[1]]
    }
  # ---- ONE GENE ----
  }else{
    if(cat != ""){
      cat_s <- strsplit(cat,"_")[[1]]
      cat_si <- which(tolower(cat_s) %in% c("high","low"))
      n_min_r <- perc_min * nrow(data)
    }else{
      n_min_r <- NULL
    }
    rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
      q <- quantiles[i]
      df <- generate_surv_df(df_o, patient_ids, exp, q)
      df <- assign_df_levels(df, data, cat, cat_si)
      if(!is.null(n_min_r)){
        # determine if meet min % samples requirement
        n_mins <- table(df[["level"]]); n_min <- min(n_mins)
        return_null <- F
        if(cat == "All"){if(length(n_mins) < 4){return_null <- T}}
        if(n_min < n_min_r){return_null <- T}
        if(return_null){
          return(NULL)
        }else{
          one_gene_cox(df,cat,q,depmap_T,p_kc)
        }
      }else{
        one_gene_cox(df,cat,q,depmap_T,p_kc)
      }
    })
  }

  rrr <- Filter(Negate(is.null), rrr)
  if(length(rrr)>0){
    #Find the minimum P-value
    pvals <- sapply(rrr, function(x) x[["least_p_value"]])
    pvals_i <- which.min(pvals)[[1]]
    #prepared index mapping df for heatmapdf
    heatmap_mapping_df <- lapply(rrr, function(x) {data.frame(x[["heatmap_new_rows"]])})
    heatmap_mapping_df <- rbindlist(heatmap_mapping_df)

    res <- rrr[[pvals_i]]
    df_most_significant <- res[["df_most_significant"]]
    # assign levels in order
    df_most_significant$level <- factor(df_most_significant$level)
    lels <- levels(df_most_significant$level)
    if(cat == "All"){
      lel_name <- lels[grepl("Low",lels)]
      lel_name <- lel_name[grepl("Other",lel_name)]
      if (identical(lel_name,character(0))) lel_name <- "Low_Low"
      df_most_significant$level <- relevel(df_most_significant$level, ref = lel_name)
    }
    cutoff_most_significant <- res[["cutoff_most_significant"]]
    least_hr <- res[["least_hr"]]
    least_p <- res[["least_p_value"]]

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
    # ---- PERMUTATION ----
    if(least_p < 0.1){
      if(num == 1){
        rrr_perm <- mclapply(1:n_perm,mc.cores = nCores,function(ii){
          df_o_new <- df_o[idx.mat[,ii],]
          df_o_new$patient_id <- patient_ids

          ppp <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
            q <- quantiles[i]
            df <- generate_surv_df(df_o_new, patient_ids, exp, q)
            df <- assign_df_levels(df, data, cat, cat_si)
            results <- one_gene_cox(df,cat,q,depmap_T,p_kc,new_row_T=F)
            return(results[["least_p_value"]])
          })

          ppp <- Filter(Negate(is.null), ppp)
          return(min(unlist(ppp)))
        })
      }else if(num > 1){
        if(search_mode == "heuristic"){

          rrr_perm <- mclapply(1:n_perm,mc.cores = nCores,function(ii){
          #rrr_perm <- lapply(1:n_perm,function(ii){
            df_o_new <- df_o[idx.mat[,ii],]
            df_o_new$patient_id <- patient_ids
            rrr <- two_gene_heuristic(quantiles, quantiles2,
                                      df_o_new, patient_ids, exp,
                                      df_o_new, patient_ids2, exp2
                                      ,gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores,heatmap_df = heatmap_df)
            rrr <- Filter(Negate(is.null), rrr)
            res <- find_minP_res(rrr)
            heatmap_new_rows <- assemble_new_rows(rrr)
            list(res[["least_p_value"]],heatmap_new_rows)
          })
        }else if(search_mode == "exhaustive"){
          rrr_perm <- mclapply(1:n_perm,mc.cores = nCores,function(ii){
            df_o_new <- df_o[idx.mat[,ii],]
            df_o_new$patient_id <- patient_ids
            #TODO:CHANGE BACK TO MCLAPPLY
            rrr <- mclapply(seq_along(quantiles),mc.cores = nCores,function(i){
            #rrr <- lapply(seq_along(quantiles),function(i){
              q <- quantiles[i]
              df <- generate_surv_df(df_o_new, patient_ids, exp, q)
              two_gene_cox(q, quantiles2, df_o_new, patient_ids2, exp2, df, gp, gps, other_gp, n_min_r, p_kc, depmap_T, nCores,heatmap_df = heatmap_df)
            })
            rrr <- Filter(Negate(is.null), rrr)
            res <- find_minP_res(rrr)
            res[["least_p_value"]]
          })
        }
      }

      # permutation adjusted P value
      rrr_perm <- Filter(Negate(is.null), rrr_perm)
      if(length(rrr_perm)>0){
        #assemble here
        # index_df <- lapply(rrr_perm, `[[`, 2)
        # index_df <- do.call("rbind", index_df)

        #Find all the P-value
        pvals_perm <- unlist(lapply(rrr_perm, `[[`, 1))# unlist(rrr_perm)
        p_adj <- sum(pvals_perm <= least_p) / n_perm
        least_error <- 1/n_perm
        if(p_adj < least_error) p_adj <- paste0("< ",least_error)
      }else{
        p_adj <- NULL
      }
    }else{
      p_adj <- NULL
    }


    # proceed only if enough data
    if(is.null(df_most_significant)){
      return(NULL)
    }else{
      #update the original heatmap_df here
      #annotation column is for showing the texts on top of heatmap, so it needs a separate column
      heatmap_mapping_df$annotation <- heatmap_mapping_df$p_value
      heatmap_df$p_value <- update_heatmap_column(target_df = heatmap_df,index_df = heatmap_mapping_df)
      heatmap_df$annotation <- update_heatmap_column(target_df = heatmap_df,index_df = heatmap_mapping_df, column_name = "annotation")
      heatmap_df$hr <- update_heatmap_column(target_df = heatmap_df,index_df = heatmap_mapping_df, column_name = "hr")


      results <- list(
        df = df_most_significant,
        cutoff = cutoff_most_significant
        ,hr = least_hr
        ,p_df = p_df
        ,heatmap_df = heatmap_df
        ,p.adj = p_adj
      )
      return(results)
    }
  }else{
    return(NULL)
  }
}

# perform interaction survival analysis on continuous-categorical combinations
cal_conti_cat_interaction <- function(x,gp_r,df_list){
  results <- get_info_most_significant_rna(rv[["dataF1"]], rv[["minF1"]], rv[["maxF1"]], rv[["stepF1"]], data2=rv[["dataF2"]], cat=gp_r)
  rv[["padj_perm"]] <- results[["p.adj"]]
  if("p_df" %in% names(results)){
    rv[["quantile_graph"]][[1]] <- results[["p_df"]]
  }
  # extract most significant df
  df <- results[["df"]]
  df_1 <- df
  df_lels <- as.character(df[["level"]])
  if(gp_r == "All"){
    df_1[["level"]] <- sapply(df_lels, function(x) strsplit(x,"_")[[1]][2])
    df[["level"]] <- sapply(df_lels, function(x) strsplit(x,"_")[[1]][1])
    # assign levels
    df_1[["level"]] <- factor(df_1[["level"]]); df[["level"]] <- factor(df[["level"]])
    df_1[["level"]] <- relevel(df_1[["level"]], ref = "Low")
    df[["level"]] <- relevel(df[["level"]], ref = sort(levels(df[["level"]]),decreasing = T)[1])
  }else{
    df_1[["level"]] <- df[["level.x"]]
    df[["level"]] <- df[["level.y"]]
  }

  df_list[[1]] <- df_1; rv[["df_1"]] <- df_1
  rv[["cutoff_1"]] <- paste0("<b>",results[["cutoff"]],"</b>")
  rv[["cutoff_all"]] <- paste0("#",x,": ",rv[["cutoff_1"]])
  # no of cases in each group
  lels_1 <- levels(df_1$level)
  rv[["lels_1"]] <- lapply(lels_1, function(x) table(df_1$level == x)["TRUE"] %>% unname(.))
  names(rv[["lels_1"]]) <- lels_1
  # perform survival analysis
  rv[["cox_1"]] <- cal_surv_rna(df_1,1,rv[["minF1"]], rv[["maxF1"]], rv[["stepF1"]],iter_mode=rv[["iter_mode_value1"]])
  return(list(df,df_list,results))
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

      # if(iter_mode){
      #   p_diff_adj <- correct_p(p_diff,min,max,step)
      # }else{
      #   p_diff_adj <- NULL
      # }

      results = list(
        df = df
        ,fit = surv_diff
        # ,stats = km.stats
        ,lels = lels
        # ,hr = "NA"
        ,p = p_diff
        ,p.adj = rv[["padj_perm"]]
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

      # # multiple p correction
      # if(iter_mode){
      #   p.km.adj <- correct_p(p.km,min,max,step)
      #   p.cox.adj <- correct_p(cox.stats$logtest[3],min,max,step)
      #   #   sapply(cox.stats$coefficients[,5], function(x){
      #   #   correct_p(as.numeric(x),min,max,step)
      #   # })
      # }else{
      #   p.km.adj <- NULL
      #   p.cox.adj <- NULL
      # }


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
          ,p.adj = rv[["padj_perm"]]
        )
        ,cox = list(
          df = new_df,
          fit = cox.fit,
          stats = cox.stats
          ,lels = lels
          ,hr = hr.cox
          ,p = p.cox
          ,p.adj = rv[["padj_perm"]]
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
  cells <- tibble(patient_id = patient_ids) %>% dplyr::inner_join(df, by = "patient_id")
  cells <- paste0(cells[["CCLE_Name"]],"|",cells[["gender"]])
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
