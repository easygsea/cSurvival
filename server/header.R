# available databases
observeEvent(input$db_download,{
  showModal(
    modalDialog(
      title = h3(HTML("cSurvival source data")),
      fluidRow(
        style = "font-size:120%;",
        column(
          12,
          p(paste0("cSurvival aims to provide reliable and integrated resources for both clinicians and experimental biologists to evaluate prognostic values of biomarkers and inspect interactions between genes, loci, pathways/gene sets, and drugs in cancers."))
          ,p(paste0("For each dataset, we controlled for quality by removing flagged (e.g. tumor tissue origin incorrect, unacceptable prior treatment, prior malignancy, does not meet study protocol, subject withdrew consent, failed QC) samples",
                    ", extracted primary tumor data (Sample Type Code 01; for TCGA-SKCM, 01 or 06) to suit the purpose of survival analysis, and transformed the data into a standardized format for customizable and reproducible studies."))
          ,br()
        ),
        tabBox(
          width = 12,
          tabPanel(
            strong("TCGA PanCanAtlas"),
            strong(style = "color:#CF5C78;","Source data:")
            ,tags$li(HTML("Liu, J., Lichtenberg, T., Hoadley, K.A., Poisson, L.M., Lazar, A.J., Cherniack, A.D., Kovatich, A.J., Benz, C.C., Levine, D.A., Lee, A.V. and Omberg, L., 2018. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. Cell, 173(2), pp.400-416. <a href='https://gdc.cancer.gov/node/977' target='_blank'>https://gdc.cancer.gov/node/977</a>"))
            ,tags$li(HTML("Hoadley, K.A., Yau, C., Hinoue, T., Wolf, D.M., Lazar, A.J., Drill, E., Shen, R., Taylor, A.M., Cherniack, A.D., Thorsson, V. and Akbani, R., 2018. Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. Cell, 173(2), pp.291-304. <a href='https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018' target='_blank'>https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018</a>"))
            ,br()
            ,strong(style = "color:#CF5C78;","cSurvival curated data:")
            ,tags$li(HTML(paste0("<b>Flagged</b> (e.g. redaction) <b>cases with annotations</b>: ",dlink(paste0(pro_dir,"977/flagged_cases.tsv"),"flagged_cases.tsv"))))
            ,br()
            ,selectizeInput(
              "project_tcga",
              div(style = "color:#CF5C78;", "Select a project to download data:")
              ,choices = projects[grepl("TCGA",names(projects))][[1]]
              ,width = "100%"
              ,options = list(
                `live-search` = TRUE,
                placeholder = "Type to search ..."
                ,onInitialize = I(sprintf('function() { this.setValue(%s); }',"TCGA-ACC")) #['TCGA-LUAD','TCGA-LUSC']
              )
            )
          )
          ,tabPanel(
            strong("TARGET NCI GDC"),
            p("Source data:")
            ,tags$li("Heath, A.P., Ferretti, V., Agrawal, S., An, M., Angelakos, J.C., Arya, R., Bajari, R., Baqar, B., Barnowski, J.H., Burt, J. and Catton, A., 2021. The NCI Genomic Data Commons. Nature Genetics, 53(3), pp.257-262.")
          )
          ,tabPanel(
            strong("DepMap cell line"),
            p("Source data:")
          )
          ,tabPanel(
            strong("eVITTA gene set libraries"),
            p("Source data:")
          )
        )
      )
      ,size = "l",
      easyClose = TRUE
      ,footer = modalButton("Dismiss")
    )
  )
})