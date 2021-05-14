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

# a link icon, need to wrap in HTML
# example: HTML(paste0("Visit eVITTA at ",link_icon("evitta_link","https://tau.cmmt.ubc.ca/eVITTA/")))
link_icon <- function(id, link, title="Click to visit", icon="fas fa-external-link-alt", color="#00c0ef", style=""){
  sprintf(
    '<a href="%s" target="_blank"><i class="%s" id="%s" style = "color:%s"></i></a>'
  ,link,icon,id,color)
}

# # add a gear button
add_gear <- function(
  id, left="6.9em", top="1em", title="Click for advanced run parameters", up = F, width = "80%"
  , placement="top"
){
  div(
    style=sprintf("position: relative; align: center; left: %s; top: %s;",left, top),
    dropdownButton(
      circle = TRUE, status = "info",
      size = "xs",
      icon = icon("gear"),# class = "opt"),
      up = up, width = width,
      tooltip = tooltipOptions(title = title, placement = placement),
      
      fluidRow(
        uiOutput(id)
      )
    )
  )
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
    gs_lib_genes_id <- paste0("gs_lgs_",x) # verbatimTextOutput on gs genes
    gs_gene_id <- paste0("gs_lg_",x); gs_gene_id_q <- paste0(gs_gene_id,"_q") # gene to search
    gs_gene_genes_id <- paste0("gs_lgg_",x) # verbatimTextOutput on input genes to filter GS
    # manual gene input
    gs_manual_id <- paste0("gs_m_",x); gs_manual_id_q <- paste0(gs_manual_id,"_q")
    gs_manual_btn_id <- paste0("add_btn_",x)
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
          choices = c("Gene or locus" = "g", "Gene set (GS)" = "gs")
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
            choices = data_types(),
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
            ,options = list(
              placeholder = g_placeholder
              ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[g_ui_id]]))
            )
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
                selectizeInput(
                  gs_db_id,
                  HTML(paste0(x,".3a. Select database:"),add_help(gs_db_id_q))
                  ,choices=gmt_dbs
                  ,selected=rv[[gs_db_id]]
                  ,options = list(
                    # `live-search` = TRUE,
                    placeholder = 'Type to search ...'
                    ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_db_id]]))
                  )
                )
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  selectizeInput(
                    gs_lib_id,
                    HTML(paste0(x,".3b. Select gene set (GS):"),add_help(gs_lib_id_q))
                    ,choices=names(rv[[paste0("gmts_tmp",x)]])
                    ,selected=rv[[gs_lib_id]]
                    ,options = list(
                      # `live-search` = TRUE,
                      placeholder = rv[[paste0("gs_placeholder",x)]]
                      ,onInitialize = I(sprintf('function() { this.setValue("%s"); }',rv[[gs_lib_id]]))
                    )
                  )
                )
                
              )
              ,conditionalPanel(
                condition = sprintf("input.%s != ''", gs_db_id),
                column(
                  4,
                  searchInput(
                    gs_gene_id,
                    HTML(paste0("(Optional) gene set filtering by gene combination:", add_help(gs_gene_id_q))),
                    value = rv[[gs_gene_id]],
                    placeholder = "Enter genes in HUGO symbol format",
                    btnSearch = icon("search"),
                    btnReset = icon("remove"),
                    width = "100%"
                  )
                  
                )
              )
            )
            ,fluidRow(
              column(
                    8,
                    conditionalPanel(
                      condition = sprintf("input.%s != ''", gs_lib_id),
                      span(verbatimTextOutput(gs_lib_genes_id), style = rv$verbTxtStyle1)
                    )
              )
              ,column(
                  4,
                  conditionalPanel(
                    condition = sprintf("input.%s != ''", gs_gene_id),
                    span(verbatimTextOutput(gs_gene_genes_id), style = rv$verbTxtStyle2)
                )
                
              )
            )
          )
          ,conditionalPanel(
            condition = sprintf("input.%s=='manual'", gs_mode_id),
            textAreaInput(
              gs_manual_id,
              HTML(paste0(x,".3. Enter your genes:"),add_help(gs_manual_id_q))
              ,value = rv[[gs_manual_id]]
              ,placeholder = "Type to enter..."
            )
            ,span(verbatimTextOutput(gs_genes_id), style = paste0(rv$verbTxtStyle1))
            ,bsButton(gs_manual_btn_id,tags$strong("Submit"),style = "warning")
          )
        )
        
        # tooltip for data category
        ,bsTooltip(cat_id_q, HTML("<b>Gene or locus</b>: To study if the expression level, mutational status, copy number, or methylation level of a gene or locus correlates with poorer/better survival.<br><b>Gene set</b>: To study if the average expression level of a gene set correlates with cancer survival, e.g. genes in the same pathway, TF targets, drug targets, miRNA targets, interacting proteins, or user-defined list of genes.")
                   ,placement = "right")
        ,bsTooltip(db_id_q, HTML("To study if cancer survival is associated with a gene\\'s expression level, mutational status, copy number variation; a microRNA\\'s expression; or the methylation level of a DNA segment.")
                   ,placement = "right")
        ,bsTooltip(g_ui_id_q, HTML(paste0("Search and select. If a gene or locus is not found, try its alias names."
                                          ," If still not found, it means its expression/alteration is barely detected in the selected cancer project."
                                          ,"Or, if pan-cancer analysis, its expression/alteration is not detected in all selected projects."))
                   ,placement = "right")
        ,bsTooltip(gs_mode_id_q, HTML("Select <b>Library</b> to analyze a pathway, a biological process, a cellular location, a transcriptional factor, a drug, or a gene\\'s interacting partners.<br>Alternatively, select <b>Manual</b> to enter your own list of genes.")
                   ,placement = "right")
        ,bsTooltip(gs_db_id_q, HTML("For full list of available gene set databases, visit <b>easyGSEA User Guide</b>.")
                   ,placement = "right")
        ,bsTooltip(gs_lib_id_q, HTML("Search for keywords (e.g. glycolysis, chemokine, tor signaling) and select the one of interest.")
                   ,placement = "right")
        ,bsTooltip(gs_gene_id_q, HTML(paste0("To filter out gene sets that contains your gene(s) of interest, in HUGO symbol format, delimited by \"&\" (and) or \"|\" (or). | is evaluated before &. Example: MYC&TP53|BCL2&BRCA1|BRCA2 is evaluated as MYC&(TP53|BCL2)&(BRCA1|BRCA2), which means MYC and (TP53 or BCL2) and (BRCA1 or BRCA2)."
                                             ," A maximum of 10 genes are supported."))
                   ,placement = "top")
        ,bsTooltip(gs_manual_id_q,HTML("Newline-, space- or comma-delimited")
                   ,placement = "right")
        ,radioTooltip(id = db_id, choice = "rna", title = HTML("Gene expression level quantified by RNA-seq"))
        ,radioTooltip(id = db_id, choice = "snv", title = HTML("Simple Nucleotide Variation (SNV)"))
        ,radioTooltip(id = db_id, choice = "cnv", title = HTML("Copy Number Variation (CNV)"))
        ,radioTooltip(id = db_id, choice = "mir", title = HTML("microRNA expression level"))
        ,radioTooltip(id = db_id, choice = "met", title = HTML("DNA methylation level"))
        ,radioTooltip(id = db_id, choice = "rrpa", title = HTML("Reverse-phase protein array (RPPA)"))
      )
      
    )
    
  })
  
  # create the UI
  do.call(tagList, ui)
}

#======================================================================#
####           function to generate dynamic run parameter UI     ####
#======================================================================#
plot_run_ui <- function(n){
  if(n == 1){col_w <- 12}else{col_w <- 6}
  
  ui <- lapply(1:n, function(x){
    
    iter_id <- paste0("iter_",x); iter_id_q <- paste0(iter_id,"_q")
    
    lower_id <- paste0("lower_",x); lower_id_q <- paste0(lower_id,"_q")
    higher_id <- paste0("upper_",x); higher_id_q <- paste0(higher_id,"_q")
    step_id <- paste0("step_",x); step_id_q <- paste0(step_id,"_q")
    
    clow_id <- paste0("clow_",x); clow_id_q <- paste0(clow_id,"_q")

    col <- extract_color(x)
    
    snv_id <- paste0("snv_method_",x); snv_id_q <- paste0(snv_id,"_q")
    snv_uni_id <- paste0("snv_uni_",x); snv_uni_id_q <- paste0(snv_uni_id,"_q")
    non_id <- paste0("nonsynonymous_",x); non_id_q <- paste0(non_id,"_q")
    syn_id <- paste0("synonymous_",x); syn_id_q <- paste0(syn_id,"_q")
    
    cnv_id <- paste0("cnv_par_",x); cnv_id_q <- paste0(cnv_id,"_q")
    
    cat_id <- paste0("cat_",x); db_id <- paste0("db_",x)
    check_inputs <- function(){
      if(is.null(input[[cat_id]]) & is.null(input[[db_id]])){
        return(rv[[cat_id]] == "g" & rv[[db_id]] != "snv")
      }else if(is.null(input[[cat_id]])){
        return(rv[[cat_id]] == "g" & input[[db_id]] != "snv")
      }else if(is.null(input[[db_id]])){
        return(input[[cat_id]] == "g" & rv[[db_id]] != "snv")
      }else{
        return((input[[cat_id]] == "g" & input[[db_id]] != "snv") | input[[cat_id]] == "gs")
      }
    }
    datatype <- call_datatype(x)
    
    # preselected SNV callers
    snv_pre_a <- ifelse(
      length(rv[[snv_id]]) == 1,
      paste0('"',rv[[snv_id]],'"'),
      paste0("['",paste0(rv[[snv_id]],collapse = "','"),"']")
    )
    
    # render the UI
    column(
      col_w,align="center",
      # tags$hr(style="border: .5px solid lightgrey; margin-top: 0.5em; margin-bottom: 0.5em;"),
      wellPanel(
        style = paste0("background-color: ", col, "; border: .5px solid #fff;"),
        h4(paste0("Advanced run parameters for Analysis #",x), align = "center"),
        h4(paste0("(",datatype,")")),
        tags$hr(style="border-color: #c2bfb5;"),
        if(check_inputs() & ifelse_rv(db_id) != "cnv"){
          div(
            radioGroupButtons(
              inputId = iter_id,
              label = HTML(paste0("Select method to stratify high- and low- ",datatype," groups"),add_help(iter_id_q)),
              choiceNames = c("Dynamic iteration", "Manual cutoffs"),
              choiceValues = c("iter","manual"),
              selected = rv[[iter_id]],
              size = "sm",
              checkIcon = list(
                yes = icon("check-square"),
                no = icon("square-o")
              ),
              # status = "primary",
              direction = "horizontal"
            ),
            conditionalPanel(
              sprintf('input.%s == "iter"',iter_id),
              sliderTextInput(
                lower_id,
                label = HTML(paste0("Start percentile:",add_help(lower_id_q)))
                ,selected = rv[[lower_id]]
                ,choices = c(.05, .1, .15, .2, .25, .3, .35, .4, .45, .5)
                ,grid = TRUE
              )
              ,sliderTextInput(
                higher_id,
                label = HTML(paste0("End percentile:",add_help(higher_id_q)))
                ,selected = rv[[higher_id]]
                ,choices = c(.5, .55, .6, .65, .7, .75, .8, .85, .9, .95)
                ,grid = TRUE
              )
              ,sliderTextInput(
                step_id,
                label = HTML(paste0("Step size:",add_help(step_id_q)))
                ,selected = rv[[step_id]]
                ,choices = c(.01, .02, .03, .05, .1, .15, .2, .25)
                ,grid = TRUE
              )
            )
            ,conditionalPanel(
              sprintf('input.%s == "manual"',iter_id),
              sliderInput(
                clow_id,
                HTML(paste0("Cutoff percentile:",add_help(clow_id_q))),
                value = rv[[clow_id]],
                min = 10,max = 90,step=.5
              )
            )
            
            ,if(x == 1 & rv$variable_n > 1){
              if(req_filter_on(paste0("db_",2:rv$variable_n),filter="snv",target="input")){
                bsButton("toall", strong("Apply to all"), style = "warning")
              }
            }
            ,bsTooltip(iter_id_q,HTML(paste0("<b>Dynamic iteration</b>: Determine the optimal cutoff by searching for the percentile yielding the lowest P-value"
                                             ,"<br><b>Manual cutoffs</b>: Manually enter the cutoffs for high- and low-",datatype," groups"))
                       ,placement = "top")
            ,bsTooltip(lower_id_q,HTML("The percentile to start iteration")
                       ,placement = "right")
            ,bsTooltip(higher_id_q,HTML("The percentile to end iteration")
                       ,placement = "right")
            ,bsTooltip(step_id_q,HTML(paste0("Step size to iterate to find the optimum cutoff"))
                       ,placement = "right")
            ,bsTooltip(clow_id_q,HTML("Cases &le; the cutoff will be classified as low, while those &gt; the cutoff will be classified as high")
                       ,placement = "right")
          )
        }else if(input[[cat_id]] == "g" & input[[db_id]] == "cnv"){
          div(
            radioGroupButtons(
              inputId = cnv_id,
              label = HTML(paste0("Select group to analyze:"),add_help(cnv_id_q)),
              choiceNames = c("Automatic", "Copy number gain", "Copy number loss"), #,"Copy number gain and loss"
              choiceValues = c("auto","gain","loss"), #,"both"
              selected = rv[[cnv_id]],
              size = "sm",
              checkIcon = list(
                yes = icon("check-square"),
                no = icon("square-o")
              ),
              # status = "primary",
              direction = "horizontal"
            )
            ,bsTooltip(cnv_id_q,HTML(paste0("<b>Automatic</b>: Automatically determines whether copy number gain and/or loss results in more significant survival difference"
                                            ,"<br><b>Copy number gain</b>: To compare cases with copy number gain with the rest of the population"
                                            ,"<br><b>Copy number loss</b>: To compare cases with copy number loss with the rest of the population"
                                            ,"<br><b>Copy number gain and loss</b>: To compare cases with copy number gain and loss, respectively, with the rest of the population"
            ))
            ,placement = "top")
          )
        }else{
          
          div(
            if(rv$tcga){
              div(
                # mutation caller options
                selectizeInput(
                  snv_id
                  ,label = HTML(paste0("Select somatic mutation caller(s):",add_help(snv_id_q)))
                  ,choices = snv_algorithms
                  ,selected = rv[[snv_id]]
                  ,multiple = T
                  ,options = list(
                    # `live-search` = TRUE,
                    placeholder = 'Type to search ...'
                    ,onInitialize = I(sprintf('function() { this.setValue(%s); }',snv_pre_a)))
                )
                ,bsTooltip(snv_id_q
                           ,HTML("Algorithm(s) used for calling somatic nucleotide variations. Multiple selections are allowed")
                                 ,placement = "top")
                # intersect or union
                ,conditionalPanel(
                  sprintf('input.%s.length > 1',snv_id),
                  radioGroupButtons(
                    inputId = snv_uni_id,
                    label = HTML(paste0("Select method to handle results by different callers: ",add_help(snv_uni_id_q))),
                    choices = c(
                      "Intersect" = "int"
                      ,"Union" = "uni"
                    ),
                    selected = rv[[snv_uni_id]],
                    size = "sm",
                    checkIcon = list(
                      yes = icon("check-square"),
                      no = icon("square-o")
                    ),
                    # status = "primary",
                    direction = "horizontal"
                  )
                )
                ,bsTooltip(snv_uni_id_q,HTML(paste0(
                  "When multiple callers are selected, select <b>Intersect</b> to extract consensus results for analysis, or <b>Union</b> to find a positive hit in any algorithm."
                )),placement = "top")
              )
            }
            
            # non-silent variants classifications
            ,selectizeInput(
              non_id
              ,label = HTML(paste0("Select variants of interest (Mutated): ",add_help(non_id_q)))
              ,choices = variant_types
              ,selected = rv[[non_id]]
              ,multiple = T
            )
            ,bsTooltip(non_id_q,HTML("By default, variants to be classified as High/Moderate variant consequences are selected. Adjust to suit the purpose of your study. For more information about variant classification: http://uswest.ensembl.org/Help/Glossary?id=535")
                       ,placement = "top")
            
            # silent variants classifications
            ,selectizeInput(
              syn_id
              ,label = HTML(paste0("Select control group (Other): ",add_help(syn_id_q)))
              ,choices = variant_types
              ,selected = rv[[syn_id]]
              ,multiple = T
            )
            ,bsTooltip(syn_id_q,HTML("By default, wild-type (WT) and variants to be classified as Low/No variant consequences are selected. Adjust to suit the purpose of your study. For more information, visit http://uswest.ensembl.org/Help/Glossary?id=535")
                       ,placement = "top")
            
            ,if(x == 1 & rv$variable_n > 1){
              if(req_filter_on(paste0("db_",2:rv$variable_n),filter="snv",target="input",mode="unequal")){
                bsButton("toall_m", strong("Apply to all"), style = "warning")
              }
            }
          )
        }
      )
    )
  })

  do.call(tagList, ui)
  
}