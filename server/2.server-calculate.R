# start analysis when users click the btn
observeEvent(input$confirm,{
  if(input$project == ""){
    shinyalert("Please select a project to begin your analysis")
  }else{
    withProgress(value = 1, message = "Performing analysis. Please wait a minute ...",{
      lapply(1:rv$variable_n, function(x){
        g_ui_id <- paste0("g_",x)
        a_range <- 2:(length(rv[[paste0("genes",x)]])+1)
        col_to_drop <- a_range[input[[g_ui_id]] != rv[[paste0("genes",x)]]]

        data <- fread(paste0(rv$indir,"df_gene_scale.csv"),sep=",",header=T,drop = col_to_drop)
        print(x)
        print(head(data))
      })
    })
  }
})