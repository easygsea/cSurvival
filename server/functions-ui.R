#======================================================================#
####                     general UI functions                       ####
#======================================================================#
# radiotooltip
radioTooltip <- function(id, choice, title, placement = "top", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

#======================================================================#
####          function to generate dynamic 1 2, or more G/GS UIs     ####
#======================================================================#
extract_color <- function(x, cols=bcols){
  bcol_n <- x %% 2; if(bcol_n == 0){bcol_n <- 2}
  return(cols[bcol_n])
}
plot_ui <- function(n){
  # automatically adjust column width according to # of analysis selected
  if(n == 1){col_w <- 12}else{col_w <- 6}
  
  # create the UI list
  ui <- lapply(1:n, function(x){
    # color number for the wellpanel
    col <- extract_color(x)
    # category to analyze
    cat_id <- paste0("cat_",x); cat_id_q <- paste0(cat_id,"_q")
    # type of db to analyze
    db_id <- paste0("db_",x); db_id_q <- paste0(db_id,"_q")
    # gene to analyze
    g_ui_id <- paste0("g_",x); g_ui_id_q <- paste0(g_ui_id,"_q")
    # gene set to analyze
    gs_mode_id <- paste0("gs_mode_",x); gs_mode_id_q <- paste0(gs_mode_id,"_q")
    gs_db_id <- paste0("gs_db_",x); gs_db_id_q <- paste0(gs_db_id,"_q")
    gs_lib_id <- paste0("gs_l_",x); gs_lib_id_q <- paste0(gs_lib_id,"_q")
    gs_lib_genes_id <- paste0("gs_lgs_",x); gs_lib_genes_id_q <- paste0(gs_lib_genes_id,"_q") # verbatimTextOutput
    gs_gene_id <- paste0("gs_lg_",x); gs_gene_id_q <- paste0(gs_gene_id,"_q") # gene to search
    # manual gene input
    gs_manual_id <- paste0("gs_m_",x); gs_manual_id_q <- paste0(gs_manual_id,"_q")
    gs_genes_id <- paste0("gs_mg_",x)

    # the UI
    column(
      col_w,
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", col, "; border: .5px solid #fff;"),
        h4(paste0("Analysis #",x), align = "center"),
        prettyRadioButtons(
          cat_id,
          label = HTML(paste0(x,".1. Select data category:",add_help(cat_id_q))),
          choices = c("Gene" = "g", "Gene set" = "gs")
          ,selected = rv[[cat_id]]
          ,status = "danger"
          ,icon = icon("check")
          ,shape = "curve"
          # ,plain = T
          # ,outline = T
          ,animation = "jelly"
          ,inline = T
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='g'", cat_id),
          radioGroupButtons(
            inputId = db_id,
            label = HTML(paste0(x,".2. Select type of molecular data:",add_help(db_id_q))),
            choices = c("Expression"="rna", 
                        "SNV"="snv",
                        "CNV"="cnv",
                        "miRNA"="mir",
                        "Methylation"="met"),
            status = "danger",
            selected = rv[[db_id]],
            checkIcon = list(
              yes = tags$i(class = "fa fa-check-square", 
                           style = "color: white"),
              no = tags$i(class = "fa fa-square-o", 
                          style = "color: white"))
          )
          ,selectizeInput(
            g_ui_id,
            HTML(paste0(x,".3. Select your gene of interest:",add_help(g_ui_id_q)))
            ,choices=c()
            ,selected=rv[[g_ui_id]]
            ,width = "100%"
          )
        )
        ,conditionalPanel(
          condition = sprintf("input.%s=='gs'", cat_id),
          radioButtons(
            gs_mode_id,
            HTML(paste0(x,".2. Mode of analysis:",add_help(gs_mode_id_q)))
            ,choices = c("Library"="lib","Manual"="manual")
            ,selected = rv[[gs_mode_id]]
            ,inline = T
          )
          ,conditionalPanel(
            condition = sprintf("input.%s=='lib'", gs_mode_id),
            fluidRow(
              column(
                4,
                pickerInput(
                  gs_db_id,
                  HTML(paste0(x,".3a. Select database:"),add_help(gs_db_id_q))
                  ,choices=gmt_dbs
                  ,selected=rv[[gs_db_id]]
                  ,options = list(
                    `live-search` = TRUE
                    , title = 'Nothing selected'
                    # ,onInitialize = I('function() { this.setValue(""); }')
                  )
                )
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  pickerInput(
                    gs_lib_id,
                    HTML(paste0(x,".3b. Select your gene set of interest:"),add_help(gs_lib_id_q))
                    ,choices=names(rv[[paste0("gmts",x)]])
                    ,selected=rv[[gs_lib_id]]
                    ,options = list(
                      `live-search` = TRUE
                      , title = 'Nothing selected'
                      # ,onInitialize = I('function() { this.setValue(""); }')
                    )
                  )
                )
                
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  # selectizeInput(
                  #   gs_gene_id,
                  #   HTML(paste0("Filter gene sets that comprise a gene:",add_help(gs_gene_id_q)))
                  #   ,choices=""
                  #   ,options = list(
                  #     placeholder = 'Type to search ...'
                  #     ,onInitialize = I('function() { this.setValue(""); }')
                  #   )
                  # )
                  searchInput(
                    gs_gene_id,
                    HTML(paste0("Filter gene sets that comprise a gene:", add_help(gs_gene_id_q))),
                    placeholder = "Enter a gene here in HUGO symbol format",
                    btnSearch = icon("search"),
                    btnReset = icon("remove"),
                    width = "100%"
                  )
                  
                )
              )
            )
            ,conditionalPanel(
              condition = sprintf("input.%s != ''", gs_lib_id),
              fluidRow(
                column(
                  12,
                  # div(id=paste0("vtxt_anchor",x))
                  verbatimTextOutput(gs_lib_genes_id)
                )
              )
            )
          )
          ,conditionalPanel(
            condition = sprintf("input.%s=='manual'", gs_mode_id),
            textAreaInput(
              gs_manual_id,
              HTML(paste0(x,".3. Enter your genes:"),add_help(gs_manual_id_q))
              ,value = ""
              ,placeholder = "Type to enter..."
            )
          )
        )
        
        # tooltip for data category
        ,bsTooltip(cat_id_q, HTML("<b>Gene</b>: To study if the expression level, mutational status, copy number, or methylation level of a gene correlates with poorer/better prognosis.<br><b>Gene set</b>: to study if the average expression level of a gene set, e.g. genes in the same pathway, TF targets, drug targets, miRNA targets, interacting proteins, or user-defined list of genes.")
                   ,placement = "right")
        ,bsTooltip(db_id_q, HTML("To study if cancer prognosis is associated with a gene\\'s expression level, mutational status, copy number variation; a microRNA\\'s expression; or the methylation level of a DNA segment.")
                   ,placement = "right")
        ,bsTooltip(g_ui_id_q, HTML("Search and select. We currently support analysis with Entrez ID, HUGO symbol, or Ensembl gene ID")
                   ,placement = "right")
        ,bsTooltip(gs_mode_id_q, HTML("Select <b>Library</b> to analyze a pathway, a biological process, a cellular location, a transcriptional factor, a drug, or a gene\\'s interacting partners.<br>Alternatively, select <b>Manual</b> to enter your own list of genes.")
                   ,placement = "right")
        ,bsTooltip(gs_db_id_q, HTML("For full list of available gene set databases, visit <b>easyGSEA User Guide</b>.")
                   ,placement = "right")
        ,bsTooltip(gs_lib_id_q, HTML("Search for keywords (e.g. glycolysis, chemokine, tor signaling) and select the one of interest.")
                   ,placement = "right")
        ,bsTooltip(gs_gene_id_q, HTML("To filter out gene sets that contains your gene of interest, in HUGO symbol format. Click search icon to search or hit \\'Enter\\'.")
                   ,placement = "top")
        ,bsTooltip(gs_manual_id_q,HTML("Newline-, space- or comma-delimited.")
                   ,placement = "right")
        ,radioTooltip(id = db_id, choice = "rna", title = HTML("Gene expression level quantified by RNA-seq"))
        ,radioTooltip(id = db_id, choice = "snv", title = HTML("Simple Nucleotide Variation (SNV)"))
        ,radioTooltip(id = db_id, choice = "cnv", title = HTML("Copy Number Variation (CNV)"))
        ,radioTooltip(id = db_id, choice = "mir", title = HTML("microRNA expression level"))
        ,radioTooltip(id = db_id, choice = "met", title = HTML("DNA methylation level"))
      )
      
    )
    
  })
  
  # create the UI
  do.call(tagList, ui)
}

#======================================================================#
####                    Update rv according to Input               ####
#======================================================================#
# # update these into rv when selections change
update_all <- function(){
  for(x in 1:rv$variable_n){
    lst <- list(
      # category to analyze
      cat_id <- paste0("cat_",x)
      # type of db to analyze
      ,db_id <- paste0("db_",x)
      # gene to analyze
      ,g_ui_id <- paste0("g_",x)
      # gene set to analyze
      ,gs_mode_id <- paste0("gs_mode_",x)
      ,gs_db_id <- paste0("gs_db_",x)
      ,gs_lib_id <- paste0("gs_l_",x)
      ,gs_lib_genes_id <- paste0("gs_lgs_",x)
      ,gs_gene_id <- paste0("gs_lg_",x)
      # manual gene input
      ,gs_manual_id <- paste0("gs_m_",x)
      ,gs_genes_id <- paste0("gs_mg_",x)
      ,lower_id <- paste0("lower_",x)
      ,higher_id <- paste0("higher_",x)
      ,step_id <- paste0("step_",x)
    )
    
    updateRV(lst)
    
    # req(input[[gs_db_id]] != "")
    # req(rv[[gs_lib_id]] != "")
    # req(rv[[paste0("gmts",x)]])
    # req(length(rv[[paste0("gmts",x)]])>0)
    # 
    # output[[gs_lib_genes_id]] <- renderText({
    #   genes <- rv[[paste0("gmts",x)]][[input[[gs_lib_id]]]]
    #   paste0("(n=",length(genes),") ", paste0(genes, collapse = " "))
    # })
  }
}

init_rvs <- function(){
  lapply(1:rv$variable_n, function(x){
    lower_id <- paste0("lower_",x)
    higher_id <- paste0("higher_",x)
    step_id <- paste0("step_",x)
    cat_id <- paste0("cat_",x)
    db_id <- paste0("db_",x)
    
    rv[[lower_id]] <- .15
    rv[[higher_id]] <- .15
    rv[[step_id]] <- .01
    rv[[cat_id]] <- "g"
    rv[[db_id]] <- "rna"
  })
}

#======================================================================#
####           function to generate dynamic run parameter UI     ####
#======================================================================#
plot_run_ui <- function(n){
  if(n == 1){col_w <- 12}else{col_w <- 6}

  lapply(1:n, function(x){

    lower_id <- paste0("lower_",x); lower_id_q <- paste0(lower_id,"_q",x)
    higher_id <- paste0("higher_",x); higher_id_q <- paste0(higher_id,"_q",x)
    step_id <- paste0("step_",x); step_id_q <- paste0(step_id,"_q",x)
    col <- extract_color(x)
    
    column(
      col_w,align="center",
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", col, "; border: .5px solid #fff;"),
        if(rv[[paste0("cat_",x)]] == "g" & rv[[paste0("db_",x)]] != "snv"){
          div(
            h4(paste0("Run parameters for Analysis #",x), align = "center"),
            sliderTextInput(
              lower_id,
              label = HTML(paste0("Lower threshold:",add_help(lower_id_q)))
              ,selected = rv[[lower_id]]
              ,choices = c(.05, .1, .15, .2, .25, .3, .35, .4, .45)
              ,grid = TRUE
            )
            ,sliderTextInput(
              higher_id,
              label = HTML(paste0("Higher threshold:",add_help(higher_id_q)))
              ,selected = rv[[higher_id]]
              ,choices = c(.05, .1, .15, .2, .25, .3, .35, .4, .45)
              ,grid = TRUE
            )
            ,sliderTextInput(
              step_id,
              label = HTML(paste0("Step size:",add_help(step_id_q)))
              ,selected = rv[[step_id]]
              ,choices = c(.01, .02, .03, .05, .1, .15, .2, .25)
              ,grid = TRUE
            )
            ,if(x == 1){
              bsButton("toall", "Apply to all", style = "warning")
            }
            ,bsTooltip(lower_id_q,HTML("The percentile to start analysis.")
                       ,placement = "right")
            ,bsTooltip(higher_id_q,HTML("The percentile to end analysis.")
                       ,placement = "right")
            ,bsTooltip(step_id_q,HTML("Step size to iterate to find the optimum threshold to separate high from low expressions.")
                       ,placement = "right")
          )
        }else{
          h4(paste0("Mutation analysis does not need parameter adjustment (Analysis #)",x), align = "center")
        }
        
      )
    )
  })
}