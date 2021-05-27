# available databases
db_dir <- "https://tau.cmmt.ubc.ca/cSurvival/project_data/"
observeEvent(input$db_download,{
  showModal(
    modalDialog(
      title = h3(HTML("cSurvival source data")),
      fluidRow(
        column(
          12, style = "font-size:120%;",
          p(paste0("cSurvival aims to provide reliable and integrated resources for both clinicians and experimental biologists to evaluate prognostic values of biomarkers and investigate interactions between genes, loci, pathways/gene sets, and drugs in cancers."))
          ,p(paste0("For each dataset, we controlled for quality by flagging problematic samples (e.g. tumor tissue origin incorrect, unacceptable prior treatment, prior malignancy, does not meet study protocol, subject withdrew consent, failed QC)",
                    ", extracted primary tumor data (Sample Type Code 03 for TCGA-LAML and TARGET-AML; 01 or 06 for TCGA-SKCM; 01 for others) to suit the purpose of survival analysis, and transformed the data into a standardized format for customizable and reproducible studies."))
          ,br()
        ),
        tabBox(id="db_panel",
          width = 12,
          tabPanel(
            strong_h4("TCGA PanCanAtlas"),
            red_title("Source data:")
            ,tags$li(HTML("Liu, J., Lichtenberg, T., Hoadley, K.A., Poisson, L.M., Lazar, A.J., Cherniack, A.D., Kovatich, A.J., Benz, C.C., Levine, D.A., Lee, A.V. and Omberg, L., 2018. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. Cell, 173(2), pp.400-416. <a href='https://gdc.cancer.gov/node/977' target='_blank'>https://gdc.cancer.gov/node/977</a>"))
            ,tags$li(HTML("Hoadley, K.A., Yau, C., Hinoue, T., Wolf, D.M., Lazar, A.J., Drill, E., Shen, R., Taylor, A.M., Cherniack, A.D., Thorsson, V. and Akbani, R., 2018. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. Cell, 173(2), pp.291-304. <a href='https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018' target='_blank'>https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018</a>"))
            ,br(),default_hr(),br()
            ,red_title("Flagged data:")
            ,tags$li(HTML(paste0("<b>Flagged cases</b>",add_help("flagged_case_q")," with annotations: ",dlink(paste0(db_dir,"977/flagged_cases.tsv"),"TCGA_flagged_cases.tsv"))))
            ,bsTooltip("flagged_case_q",HTML(flagged_exp),placement = "right",options = list(container = "body"))
            ,tags$li(HTML(paste0("<b>Flagged samples</b>",add_help("flagged_sample_q")," with annotations: ",dlink(paste0(db_dir,"977/all_bad.tsv"),"TCGA_flagged_samples.tsv"))))
            ,bsTooltip("flagged_sample_q",HTML(paste0(
              "Samples annotated as <b>Do_not_use</b> by Hoadley et al, Cell, 2018."
            )),placement = "right",options = list(container = "body"))
            ,br()
            # ------------- 1. curated TCGA data --------------
            ,default_hr()
            ,br()
            ,red_title("cSurvival curated clinical and omics data by project:")
            ,tagList(lapply(projects[["TCGA"]], function(x){
              x_name <- projects[["TCGA"]]
              x_name <- names(x_name)[x_name == x]
              div(
                grey_h4(x_name)
                ,tags$li(HTML(paste0("<b>Clinical outcome</b>: ",dlink(paste0(db_dir,x,"/df_survival_o.csv"),paste0(x,"_clinical.csv")))))
                ,tags$li(HTML(paste0("<b>Gene expression (RNA-seq, upper quantile normalized RESM)</b>: ",dlink(paste0(db_dir,x,"/df_gene.csv"),paste0(x,"_expression.csv")))))
                # ,tags$li(HTML(paste0("<b>Z-score-transformed gene expression (RNA-seq)</b>: ",dlink(paste0(db_dir,x,"/df_gene_scale.csv"),paste0(x,"_expression_scaled.csv")))))
                ,tags$li(HTML(paste0("<b>Simple nucleotide variation (SNV)</b>: ",dlink(paste0(db_dir,x,"/df_snv_class_977.csv"),paste0(x,"_mutation.csv")))))
                ,tags$li(HTML(paste0("<b>Copy number variation (CNV)</b>: ",dlink(paste0(db_dir,x,"/df_cnv.csv"),paste0(x,"_copy_number.csv")))))
                ,tags$li(HTML(paste0("<b>MicroRNA expression</b>: ",dlink(paste0(db_dir,x,"/df_mir.csv"),paste0(x,"_mir.csv")))))
                ,tags$li(HTML(paste0("<b>Methylation levels</b>: ",dlink(paste0(db_dir,x,"/df_met.csv"),paste0(x,"_methylation.csv")))))
                ,br()
                # ,bsTooltip("clinical_o_q",HTML(paste0(
                #   "<b>OS</b>: overall survival (OS) status, 0=alive, 1=dead;"
                #   ," <b>OS.time</b>: OS time."
                #   ,"<br><b>DSS</b>: disease-specific survival (DSS) status, 0=alive, 1=dead;"
                #   ," <b>DSS.time</b>: DSS time."
                #   ,"<br><b>DFI</b>: disease-free interval (DFI) status, 0=no recurrence, 1=recurrence;"
                #   ," <b>DFI.time</b>: DFI time"
                #   ,"<br><b>PFI</b>: progression-free interval (PFI) status, 0=no progression, 1=progression;"
                #   ," <b>PFI.time</b>: PFI time"
                # )),placement = "right")
              )
            }))
          )
          ,tabPanel(
            strong_h4("TARGET NCI GDC"),
            red_title("Source data:")
            ,tags$li(HTML("Heath, A.P., Ferretti, V., Agrawal, S., An, M., Angelakos, J.C., Arya, R., Bajari, R., Baqar, B., Barnowski, J.H., Burt, J. and Catton, A., 2021. The NCI Genomic Data Commons. Nature Genetics, 53(3), pp.257-262. <a href='https://portal.gdc.cancer.gov' target='_blank'>https://portal.gdc.cancer.gov</a>"))
            ,br(),default_hr(),br()
            # ------------- 3. curated TARGET data --------------
            ,red_title("cSurvival curated clinical and omics data by project:")
            ,tagList(lapply(projects[["TARGET"]], function(x){
              x_name <- projects[["TARGET"]]
              x_name <- names(x_name)[x_name == x]
              div(
                grey_h4(x_name)
                ,tags$li(HTML(paste0("<b>Clinical outcome</b>: ",dlink(paste0(db_dir,x,"/df_survival.csv"),paste0(x,"_clinical.csv")))))
                ,tags$li(HTML(paste0("<b>Gene expression (RNA-seq)</b>: ",dlink(paste0(db_dir,x,"/df_gene.csv"),paste0(x,"_expression.csv")))))
                # ,tags$li(HTML(paste0("<b>Z-score-transformed gene expression (RNA-seq)</b>: ",dlink(paste0(db_dir,x,"/df_gene_scale.csv"),paste0(x,"_expression_scaled.csv")))))
                ,if(x == "TARGET-ALL-P3" | x == "TARGET-NBL" | x == "TARGET-AML" | x == "TARGET-WT"){
                  tags$li(HTML(paste0("<b>Simple nucleotide variation (SNV)</b>: ",dlink(paste0(db_dir,x,"/df_snv_class.csv"),paste0(x,"_mutation.csv")))))
                }
                ,if(x == "TARGET-CCSK" | x == "TARGET-AML" | x == "TARGET-ALL-P2" | x == "TARGET-OS"){
                  tags$li(HTML(paste0("<b>Copy number variation (CNV)</b>: ",dlink(paste0(db_dir,x,"/df_cnv.csv"),paste0(x,"_copy_number.csv")))))
                }
                ,if(x == "TARGET-ALL-P3" | x == "TARGET-AML" | x == "TARGET-ALL-P2" | x == "TARGET-RT" | x == "TARGET-WT"){
                  tags$li(HTML(paste0("<b>MicroRNA expression</b>: ",dlink(paste0(db_dir,x,"/df_mir.csv"),paste0(x,"_mir.csv")))))
                }
                ,br()
              )
            }))
          )
          ,tabPanel(
            strong_h4("DepMap cell line"),
            red_title("Source data:")
            ,tags$li(HTML("DepMap, Broad (2021): DepMap 21Q2 Public. figshare. Dataset. <a href='https://doi.org/10.6084/m9.figshare.14541774.v2' target='_blank'>https://doi.org/10.6084/m9.figshare.14541774.v2</a>"))
            ,tags$li(HTML("Ghandi, M., Huang, F.W., Jané-Valbuena, J., Kryukov, G.V., Lo, C.C., McDonald, E.R., Barretina, J., Gelfand, E.T., Bielski, C.M., Li, H. and Hu, K., 2019. Next-generation characterization of the cancer cell line encyclopedia. Nature, 569(7757), pp.503-508."))
            ,tags$li(HTML("Meyers, R.M., Bryan, J.G., McFarland, J.M., Weir, B.A., Sizemore, A.E., Xu, H., Dharia, N.V., Montgomery, P.G., Cowley, G.S., Pantel, S. and Goodale, A., 2017. Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens in cancer cells. Nature genetics, 49(12), pp.1779-1784."))
            ,tags$li(HTML("Dempster, J.M., Rossen, J., Kazachkova, M., Pan, J., Kugener, G., Root, D.E. and Tsherniak, A., 2019. Extracting biological insights from the project achilles genome-scale CRISPR screens in cancer cell lines. BioRxiv, p.720243."))
            ,tags$li(HTML("Tsherniak, A., Vazquez, F., Montgomery, P.G., Weir, B.A., Kryukov, G., Cowley, G.S., Gill, S., Harrington, W.F., Pantel, S., Krill-Burger, J.M. and Meyers, R.M., 2017. Defining a cancer dependency map. Cell, 170(3), pp.564-576."))
            ,br(),default_hr(),br()
            # ------------- 3. curated DepMap data --------------
            ,red_title("cSurvival reformated data:")
            ,div(
              grey_h4("Gene perturbation effect")
              ,tags$li(HTML(paste0("<b>CRISPR-Cas9 gene knockout effect</b>: ",dlink(paste0(db_dir,"DepMap/DepMap-CRISPR.csv"),paste0("DepMap-CRISPR.csv")))))
              ,tags$li(HTML(paste0("<b>RNAi gene knockdown effect</b>: ",dlink(paste0(db_dir,"DepMap/DepMap-RNAi.csv"),paste0("DepMap-RNAi.csv")))))
              ,br(),grey_h4("Drug sensitivity")
              ,tags$li(HTML(paste0("<b>Drug sensitivity</b>: ",dlink(paste0(db_dir,"DepMap/DepMap-Drug.csv"),paste0("DepMap-RNAi.csv")))))
              ,br(),grey_h4("Cell line omics data")
              ,tags$li(HTML(paste0("<b>Gene expression (RNA-seq, TPM)</b>: ",dlink(paste0(db_dir,"DepMap/df_gene.csv"),paste0("DepMap_gene_expression.csv")))))
              ,tags$li(HTML(paste0("<b>Simple nucleotide variation (SNV)</b>: ",dlink(paste0(db_dir,"DepMap/df_snv_class.csv"),paste0("DepMap_mutation.csv")))))
              ,tags$li(HTML(paste0("<b>Copy number variation (CNV)</b>: ",dlink(paste0(db_dir,"DepMap/df_cnv.csv"),paste0("DepMap_copy_number.csv")))))
              ,tags$li(HTML(paste0("<b>Protemoics</b>: ",dlink(paste0(db_dir,"DepMap/df_proteomic.csv"),paste0("DepMap_proteomics.csv")))))
            )
          )
          ,tabPanel(
            strong_h4("eVITTA gene set libraries"),
            red_title("Source data:")
            ,tags$li(HTML("Cheng, X., Yan, J., Liu, Y., Wang, J. and Taubert, S., 2021. eVITTA: a web-based visualization and inference toolbox for transcriptome analysis. Nucleic Acids Research. <a href='https://tau.cmmt.ubc.ca/eVITTA' target='_blank'>https://tau.cmmt.ubc.ca/eVITTA</a>"))
            ,br(),default_hr(),br()
            # ------------- 4. eVITTA GSs --------------
            ,red_title("eVITTA gene set (GS) libraries by category:")
            ,tagList(lapply(names(gmt_dbs), function(db_cat){
              
              div(
                grey_h4(db_cat),
                tagList(
                  lapply(gmt_dbs[[db_cat]], function(x){
                    tags$li(HTML(paste0("<b>",x,"</b>:",dlink(retrieve_gmt_path(x),basename(retrieve_gmt_path(x))))))
                  })
                )
              )
            }))
          )
        )
      )
      ,size = "l",
      easyClose = TRUE
      ,footer = modalButton("Dismiss")
    )
  )
})

observeEvent(input$project_tcga,{req(input$project_tcga != "");rv$project_tcga <- input$project_tcga})
