observeEvent(input$btn_jump_to_geo,{
  rv$gsea_done <- ""
  # --------------- perform differential expression analysis ---------
  print("hihihi")
  withProgress(value = 0.1, message = wait_msg("Compiling data for DE analysis..."),{
    # a) gene expression matrix
    df_gene <- de_dfgene()
    
    # # run DE analysis
    incProgress(amount = 0.1, message = "Filtering samples according to study design...")
    # # design matrix
    if(rv$variable_nr == 1){
      df_design <- rv[["df_gender"]]
    }else{
      df_design <- rv[["df_all"]]
    }
    
    df_gene <- df_gene[ , (colnames(df_gene) %in% df_design$patient_id)]

    # # run DE analysis
    incProgress(amount = 0.1, message = "Filtering lowly expressed genes...")
    
    # 1) create dgelist
    y <- DGEList(counts=df_gene)
    
    # 2) filter low expressing genes
    # filter according to design matrix
    min_n <- min(table(df_design$level))
    keep <- rowSums(y$counts>1) >= min_n
    y <- y[keep,,keep.lib.sizes=TRUE]
    
    incProgress(amount = 0.1, message = "Creating design matrix...")

    # filter df_design according to count table
    df_design <- df_design %>% dplyr::filter(patient_id %in% colnames(y$counts))

    # create the design model
    design <- model.matrix(~level.x*level.y, df_design)
    
    incProgress(amount = 0.2, message = wait_msg("Performing DE analysis..."))
    
    # 4) voom directly on counts, if data are very noisy, as would be used for microarray
    v <- voom(y, design, plot=F, normalize.method="quantile")
    
    # 5) DEG analysis
    fit <- lmFit(v, design)
    fit <- eBayes(fit,trend=TRUE, robust=TRUE)
    
    # results
    results <- decideTests(fit)

    incProgress(amount = 0.1, message = "Exporting results...")
    
    # available coefficients
    coefs <- colnames(design)[-1]
    
    # name the coefficients
    lels1 <- levels(df_design$level.x) %>% rev()
    lels2 <- levels(df_design$level.y) %>% rev()
    names(coefs) <- c(
      paste0(lels1, collapse = " vs. "),
      paste0(lels2, collapse = " vs. "),
      paste0("(",paste0(lels2, collapse = " vs. "),") vs. (",paste0(lels1, collapse = " vs. "),")")
    )
    
    # export DEG table
    degss = lapply(coefs, function(x){
      topTable(fit, coef=x,sort.by="P",number=Inf)
    })
    names(degss) <- coefs
    
    # saveRDS(degss,"degss.rds")
    # saveRDS(coefs,"coefs.rds")
    # save the variables for easyGEO
    rv$variables_for_geo[['degss']] <- save_csurvival_variable(degss)
    rv$variables_for_geo[['coefs']] <- save_csurvival_variable(coefs)
    # print(rv$variables_for_geo)
  })
  
  rv$gsea_done <- "yes"
},ignoreInit = T)