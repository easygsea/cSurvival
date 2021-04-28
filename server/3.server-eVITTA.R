observeEvent(list(rv$plot_type,input$evitta),{
  req(rv$plot_type == "gsea")
  x <- rv$evitta <- input$evitta
  # --------------- perform differential expression analysis ---------
  if(x == "geo" & rv$gsea_done==""){
    print("hihihi")
    print(getwd())
    withProgress(value = 0.1, message = wait_msg("Compiling data for DE analysis..."),{
      # a) gene expression matrix
      df_gene <- de_dfgene()

      #    saveRDS(df_gene,"df_gene_all.rds")
      incProgress(amount = 0.1, message = "Creating design matrix...")
      # b) design matrix
      if(rv$variable_nr == 1){
        df_design <- rv[["df_gender"]]
      }else{
        df_design <- rv[["df_all"]]
      }
      # saveRDS(df_design,"df_design_all.rds")
      # create the design model
      design <- model.matrix(~level.x*level.y, df_design)

      # # run DE analysis
      incProgress(amount = 0.1, message = "Filtering lowly expressed genes...")

      # 1) create dgelist
      y <- DGEList(counts=df_gene)

      # 2) filter low expressing genes
      min_n <- min(table(df_design$level))
      keep <- rowSums(y$counts>1) >= min_n
      y <- y[keep,,keep.lib.sizes=TRUE]

      incProgress(amount = 0.2, message = wait_msg("Performing DE analysis..."))

      # 3) voom directly on counts, if data are very noisy, as would be used for microarray
      v <- voom(y, design, plot=F, normalize.method="quantile")

      # 4) DEG analysis
      fit <- lmFit(v, design)
      fit <- eBayes(fit,trend=TRUE, robust=TRUE)

      # results
      results <- decideTests(fit)
      print(summary(results))

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
      print(head(degss[[1]]))
      # saveRDS(degss,"degss.rds")
      # saveRDS(coefs,"coefs.rds")
      # save the variables for easyGEO
      rv$variables_for_geo[['degss']] <- save_csurvival_variable(degss)
      rv$variables_for_geo[['coefs']] <- save_csurvival_variable(coefs)
      print(rv$variables_for_geo)
    })

    rv$gsea_done <- "yes"
  }
},ignoreInit = T)