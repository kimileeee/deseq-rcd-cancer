# sourced by 'server.R'
# save as 'GTEx_exp_server.R'
# server elements 'GTEx_exp' sub tab of 'GTEx' tab

source(file.path(config$wd, "functions", "data_function.R"))
source(file.path(config$wd, "functions", "gtex_gexp_function.R"))





# welcome information -----------------------------------------------------

output$ui_gtex_exp_welcome <- shiny::renderUI({fn_gtex_exp_welcome()})

output$ui_gtex_exp_help <- shiny::renderUI({fn_gtex_exp_help()})

# generate eqtl result out ui -------------------------------------------------------

output$ui_gexp_result <- shiny::renderUI({fn_gexp_result(selected_analysis$gtex_exp)})

#  get tissue type --------------------------------------------------------
# 
# observeEvent(input$gtex_expr_submit, {
#   status$gtex_expr_submit <- TRUE
#   shinyjs::disable(id = "GTEx_tissue_submit")
#   shinyjs::enable(id = "analysis_stop")
# })

# stop analysis -----------------------------------------------------------
# observeEvent(input$analysis_stop, {
#   status$GTEx_tissue_submit <- FALSE
#   shinyjs::enable(id = "GTEx_tissue_submit")
#   GTEx_hide <- c("GTEx_exp", "GTEx_gsva")
#   hidePic(GTEx_hide) # hide pic when stop clicked!
# })


# get gene set-----------------------
# GTEx_expr_gene_list <- eventReactive(
#   eventExpr = status$analysis,
#   ignoreNULL = TRUE,
#   valueExpr = {
#     # be sure the following code run after start analysis
#     if (status$analysis == TRUE) {
#       status$gtex_expr_submit <- TRUE
#       shinyjs::disable(id = "GTEx_tissue_submit")
#       shinyjs::enable(id = "analysis_stop")
#       as.character(gene_set$match.gtex)
#     }
#   }
# )

# Selected tissue types ---------------------------------------------------
# GTEx_tissue_type <- callModule(GTEx_normal_Tissue, "GTEx_exp")
# output$GTEx_selected_tissues <- renderText(GTEx_tissue_type())



####### reset tissue type selection when click-----------
# observeEvent(input$GTEx_tissue_reset, {
#   GTEx_tissue_type <- callModule(resetGTExTissueType, "GTEx_exp")
# })
# observeEvent(input$GTEx_tissue_submit, heatmap_gsva_4_geneset(gene_set = gene_set, gtex_expr = gtex_expr, tissue_set = tissue_set))
# analysis core ----------------------------

##### get gene expression profiles in GTEx dataset######

get_gene_exp_profile <- function(gene_set = gene_set, gtex_expr = gtex_expr, tissue_set = tissue_set, filter_gene = 0) {
  gtex_expr <- gtex_expr[gtex_expr$SMTS %in% tissue_set, ]
  if (filter_gene) {
    gtex_expr %>%
      dplyr::mutate(
        expr = purrr::map(
          .x = expr,
          .f = function(.x) {
            .x %>%
              dplyr::filter(symbol %in% gene_set)
          }
        )
      ) -> gtex_gene_list_expr
  } else {
    gtex_gene_list_expr <- gtex_expr
  }
  return(gtex_gene_list_expr)
}


get_gene_mean_profile <- function(gene_set = gene_set, gtex_expr_mean = gtex_expr_mean, tissue_set = tissue_set, filter_gene = 0) {
  gtex_expr_mean <- gtex_expr_mean[gtex_expr_mean$SMTS %in% tissue_set, ]
  if (filter_gene) {
    gtex_expr_mean %>%
      dplyr::mutate(
        Mean = purrr::map(
          .x = Mean,
          .f = function(.x) {
            .x %>%
              dplyr::filter(symbol %in% gene_set)
          }
        )
      ) -> gtex_gene_list_mean
  } else {
    gtex_gene_list_mean <- gtex_expr_mean
  }
  return(gtex_gene_list_mean)
}



###### do heatmap and gsva for gene set ####
heatmap_gsva_4_geneset <- eventReactive(
  {
    status$analysis == TRUE
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      if (selected_analysis$gtex_exp == TRUE) {
        # check box-------------------------
        callModule(module = tissueTypesSelect, id = "GTEx_exp", .sctps = intersect(selected_ctyps(), gtex_data))
        callModule(module = selectAndAnalysis, id = "GTEx_exp", .id = "GTEx_exp")
        load_data_gexp()
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} start: expression profiles of gene set on GTEx dataset processing@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
        # gene_set <- GTEx_expr_gene_list()
        # tissue_set <- c("Heart", "Ovary", "Lung","Muscle","Blood","Uterus","Vagina","Breast","Skin","Testis","Colon","Stomach","Pancreas")
        # print(tissue_set)
        # print(gene_set)
        
        ##### start: draw heatmap for the gene set in GTEx dataset######
        print("start: draw heatmap for the gene set in GTEx dataset")
        gtex_gene_list_expr.mean <- get_gene_mean_profile(gene_set = gene_set$match.gtex, gtex_expr_mean = gtex_expr_mean, tissue_set = input$select_ctps, filter_gene = 1)
        if(nrow(gtex_gene_list_expr.mean)>0){
          gene_n <- nrow(gtex_gene_list_expr.mean$Mean[[1]])
          display_matrix <- data.frame(round(matrix(unlist(lapply(gtex_gene_list_expr.mean$Mean,function(.x){.x[2]})), nrow = gene_n), 2))
          colnames(display_matrix) <- gtex_gene_list_expr.mean$SMTS
          display_matrix$GeneName <- gtex_gene_list_expr.mean$Mean[[1]]$symbol
          display_matrix %>% tidyr::gather(Tissue, RPKM, -GeneName) %>% dplyr::rename("TPM"="RPKM")-> hm_4_p
          callModule(heatmap_GTEX_Plot, "GTEx_exp", data = hm_4_p, status=status, downloadname = "heatmap_GTEX_Plot")
        } else {
          .msg <- glue::glue("No GTEx expression in gene set for your selected tissue.")
          shinyBS::createAlert(
            session = session, anchorId = "GTEx_expr-no_gene_set", title = "Oops", style = "danger",
            content = .msg, append = FALSE
          )
        }
        
        print("end: draw heatmap for the gene set in GTEx dataset")
        ######## calculate and draw GSVA profiles for gene set in selected tissues in GTEx dataset############
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} start: calculating gene set on GTEx dataset@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
      }
      # print(GTEx_expr_gene_list(), "2")
      # print("2")
      ## can not get the tissue info!!!!!#
      # print(GTEx_tissue_type())

      
      #   gtex_expr <- get_gene_exp_profile(gene_set = gene_set, gtex_expr = gtex_expr, tissue_set = tissue_set,filter_gene = 0)
      #
      #   gene_set.lst <- list(atg_lys = gene_set)
      #   fn_gsva <- function(.y, gene_set.lst = gene_set.lst){
      #     .y %>%
      #       tidyr::drop_na() %>%
      #       dplyr::select( -ensembl_gene_id) -> .d
      #
      #     .d_mat <- as.matrix(.d[,-1])
      #     rownames(.d_mat) <- .d$symbol
      #     .es_dif <- gsva(.d_mat, gset.idx.list = gene_set.lst, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
      #     .es_dif %>%
      #       as.data.frame() %>%
      #       tibble::as_tibble() %>%
      #       tibble::add_column(set = "atg_lys", .before = 1) -> .d_es
      #   }
      #   print("run_gsva")
      #   gtex_expr %>%
      #     dplyr::mutate(
      #       gsva = purrr::map(
      #         .x = expr,
      #         .f = function(.x) {
      #           fn_gsva(.x, gene_set.lst = gene_set.lst)
      #         }
      #       )
      #     ) -> gtex_expr_gsva
      #   print("end_gsva")
      #   gtex_expr_gsva %>%
      #     dplyr::select(SMTS, gsva) %>%
      #     dplyr::mutate(
      #       gsva = purrr::map(
      #         .x = gsva,
      #         .f = function(.x) {
      #           .x %>%
      #             dplyr::select(-set) %>%
      #             tidyr::gather(key = barcode, value = gsva)
      #         }
      #       )
      #     ) %>% tidyr::unnest() -> plot_ready
      #
      # 	print("call box plot for gsva")
      # callModule(box_GTEx_GSVA_Plot, "GTEx_gsva", data=plot_ready)
    }
  }
)

# observe(GTEx_expr_gene_list())
observe(heatmap_gsva_4_geneset())

# source by "server.R"

source(file.path(config$wd, "functions", "data_function.R"))
source(file.path(config$wd, "functions", "tcga_expr_function.R"))

expr_clean <- NULL
survival_clean <- NULL


# toga expr welcome -------------------------------------------------------
output$ui_expr_welcome <- shiny::renderUI({fn_expr_welcome()})

output$ui_expr_help <- shiny::renderUI({fn_expr_help()})

# expr analysis result ----------------------------------------------------

output$ui_expr_result <- shiny::renderUI({fn_expr_result(selected_analysis$expr)})



# Start analysis ----------------------------------------------------------

expr_start_analysis <- function(input, output, session, .expr_clean, .survival_clean, .subtype_clean, .msg, .msg_no_result) {
  output$expr_dt_comparison <- DT::renderDataTable({expr_clean_datatable(.expr_clean)})
  if(nrow(.expr_clean)>0){
    exp_plot_height <- .expr_clean$symbol %>% unique() %>% length()*20
    output$expr_bubble_plot <- renderPlot({
    .expr_clean %>% expr_buble_plot()
    } ,height = function(){ifelse(exp_plot_height<200,200,exp_plot_height)})
    
  output$`de-picdownload` <- downloadHandler(
    filename = function() {
      paste("Differential_Expression", ".", input$`de-pictype`, sep = "")
    },
    content = function(file){
      print(file)
      ggsave(file,expr_buble_plot(.expr_clean),device = input$`de-pictype`,width = input$`de-d_width`,height = input$`de-d_height`)
    }
  )
  .msg_de_massage <- NULL
  } else {
    .msg_de_massage <- .msg
    output$expr_bubble_plot <- renderPlot({NULL})
    output$`de-picdownload` <- downloadHandler(NULL)
  }
  
  
  # survival
  if(nrow(.survival_clean)>0){
    sur_plot_height <- .survival_clean$symbol %>% unique() %>% length()*20
    output$survival <- renderPlot({.survival_clean %>% survival_bubble_plot()},height = function(){ifelse(sur_plot_height<200,200,sur_plot_height)})
    output$`sur-picdownload` <- downloadHandler(
      filename = function() {
        paste("Expression_Survival", ".", input$`sur-pictype`, sep = "")
      },
      content = function(file){
        print(file)
        ggsave(file,survival_bubble_plot(.survival_clean),device = input$`sur-pictype`,width = input$`sur-d_width`,height = input$`sur-d_height`)
      }
    )
    .msg_sur_massage <- NULL
  } else {
    .msg_sur_massage <- .msg_no_result
    output$survival <- renderPlot({NULL})
    output$`sur-picdownload` <- downloadHandler(NULL)
  }
  
  
  # subtype
  if(nrow(.subtype_clean)>0){
    sub_plot_height <- .subtype_clean$symbol %>% unique() %>% length()*20
    output$subtype <- renderPlot({.subtype_clean %>% subtype_bubble_plot()},height = function(){ifelse(sub_plot_height<200,200,sub_plot_height)})
    output$`sub-picdownload` <- downloadHandler(
      filename = function() {
        paste("Subtype", ".", input$`sub-pictype`, sep = "")
      },
      content = function(file){
        print(file)
        ggsave(file,subtype_bubble_plot(.subtype_clean),device = input$`sub-pictype`,width = input$`sub-d_width`,height = input$`sub-d_height`)
      }
    )
    .msg_sub_massage <- NULL
  } else {
    .msg_sub_massage <- .msg_no_result
    output$subtype <- renderPlot({NULL})
    output$`sub-picdownload` <- downloadHandler(NULL)
  }
  
  # message output ----
  output[["de_massage"]] <- renderUI({
    tagList(
      shiny::tags$p(.msg_de_massage, style = "color:#CD3700")
    )
  })
  output[["sur_massage"]] <- renderUI({
    tagList(
      shiny::tags$p(.msg_sur_massage, style = "color:#CD3700")
    )
  })
  output[["sub_massage"]] <- renderUI({
    tagList(
      shiny::tags$p(.msg_sub_massage, style = "color:#CD3700")
    )
  })
}


# From start analysis -----------------------------------------------------
expr_analysis <- eventReactive(
  eventExpr = status$analysis,
  ignoreNULL = TRUE,
  valueExpr = {
    # be sure the following code run after start analysis
    if (status$analysis == TRUE) {
        # load data expr ----
        processing$start_loading_start <- TRUE
        
        # Cancer types value box selection ----------------------------------------
        
        callModule(module = cancerTypesSelect, id = "expr", .sctps = intersect(selected_ctyps(), tcga_data))
        
        # Check box ---------------------------------------------------------------
        
        callModule(module = selectAndAnalysis, id = "expr", .id = "expr")
        
        updateProgressBar(session = session, id = "progressbar", value = 40, status = "danger")
        session$onFlushed(function() {progress$expr_loading <- TRUE})
        observeEvent(eventExpr = progress$expr_loading, handlerExpr = {
            if (progress$expr_loading == TRUE) {
              # load data
              load_data_expr()
              processing$start_loading_end <- TRUE
            }
          })
        
        observeEvent(processing$start_loading_end, {
          if (processing$start_loading_end == TRUE) {
            updateProgressBar(session = session, id = "progressbar", value = 70, status = "warning")
            session$onFlushed(function() {
              progress$expr_calc <- TRUE
            })
          }
        })
        
        observeEvent(eventExpr = progress$expr_calc, handlerExpr = {
            if (progress$expr_calc == TRUE) {
              
              .valid_ctps <- intersect(paired_cancer_types, selected_ctyps())
              .invalid_ctps <- setdiff(selected_ctyps(), paired_cancer_types)
              .msg_valid_ctps <- ""
              .msg_invalid_ctps <- ""
              if (length(.valid_ctps) > 1) {
                .msg_valid_ctps <- glue::glue("only {paste0(.valid_ctps, collapse = ', ')} have paired samples")
              } else if (length(.valid_ctps) == 1) {
                .msg_valid_ctps <- glue::glue("only {.valid_ctps} has paired samples")
              } else{
                .msg_valid_ctps <- "no cancer types have paired samples"
              }
              if (length(.invalid_ctps) > 1) {
                .msg_invalid_ctps <- glue::glue("The cancer types {paste0(.invalid_ctps, collapse = ', ')} don't have paired samples.")
              } else if (length(.invalid_ctps) == 1) {
                .msg_invalid_ctps <- glue::glue("The cancer type {.invalid_ctps} doesn't have paired samples.")
              } else {
                .msg_invalid_ctps <- ""
              }
              .msg <- glue::glue("In Tumor vs. Normal module, the analysis is based on paired samples in each cancer types. In your selected cancer types, {.msg_valid_ctps}. {.msg_invalid_ctps}")
              .msg_no_result <- glue::glue("No significant result of gene set in your selected cancer types: {selected_ctyps()}.")
              
              shinyBS::createAlert(
                session = session, anchorId = "expr-no_gene_set", title = "Information", style = "info",
                content = .msg, append = FALSE
              )
              
              expr %>%
                dplyr::filter(cancer_types %in% paired_cancer_types) %>%
                dplyr::filter(cancer_types %in% selected_ctyps()) %>% 
                dplyr::filter(symbol %in% gene_set$match) ->> expr_clean
              expr_survival %>% 
                dplyr::filter(cancer_types %in% selected_ctyps()) %>% 
                dplyr::filter(symbol %in% gene_set$match) ->> survival_clean
              expr_subtype %>% 
                dplyr::filter(cancer_types %in% selected_ctyps()) %>% 
                dplyr::filter(symbol %in% gene_set$match) ->> subtype_clean
              
              print(glue::glue("{paste0(rep('-', 10), collapse = '')} clean data complete @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
              # call module
              callModule(module = expr_start_analysis, id = "expr", .expr_clean = expr_clean, .survival_clean = survival_clean, .subtype_clean = subtype_clean, .msg = .msg, .msg_no_result=.msg_no_result)
              print(glue::glue("{paste0(rep('-', 10), collapse = '')} expr bubble plot complete @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
              
              processing$expr_calc_end <- TRUE
            }
          })
        
        observeEvent(processing$expr_calc_end, {
          if (processing$expr_calc_end == TRUE) {
            updateProgressBar(session = session, id = "progressbar", value = 100, status = "info")
            progress$progress_end <- TRUE
          }
        })
    }
  }
)

observe(expr_analysis())
# source by server.R
# saved as functions_server.R

# Check input gene set ----------------------------------------------------
check_gene_set <- function(.s, status = status, error = error) {
  .s %>%
    stringr::str_split(pattern = "[^[:alnum:]]+", simplify = TRUE) %>%
    .[1, ] %>%
    stringr::str_trim(side = "both") -> .ss
  #print(.ss)

  if (!dplyr::between(length(.ss), 5, 100)) {
    error$gene_set <- "The number of genes should be between 5 and 100."
    status$trigger <- if (status$trigger == TRUE) FALSE else TRUE
    status$gene_set <- FALSE
  }

  .ss
}


# Validate gene with TCGA gene symbol -------------------------------------
# can get the old version from github
validate_gene_set <- function(.v, user_dir = user_dir, user_logs = user_logs, total_gene_symbol = total_gene_symbol, status = status, error = error, gene_set = gene_set) {
  .log_file <- user_logs$gene_set
  total_gene_symbol %>%
    dplyr::mutate(alias = toupper(alias)) %>%
    dplyr::mutate(NCBI_sym = toupper(NCBI_sym))-> total_gene_symbol

  .v_dedup <- tibble::tibble(input = .v[.v != ""]) %>% unique() %>%
    dplyr::mutate(Up = toupper(input))
  .v_dedup %>% dplyr::rename("TCGA_sym" = "Up") %>%
    dplyr::inner_join(total_gene_symbol, by ="TCGA_sym") %>%
    # dplyr::filter(TCGA_sym %in% .v_dedup$Up) %>%
    # .$TCGA_sym  %>%
    unique() -> .v_tcga

  .v_dedup %>% dplyr::rename("GTEX_sym" = "Up") %>%
    dplyr::inner_join(total_gene_symbol, by ="GTEX_sym") %>%
    #dplyr::filter(GTEX_sym %in% .v_dedup$Up)  %>%
    #.$NCBI_sym %>%
    unique() -> .v_gtex

  # gene id not match in gtex data
  setdiff(.v_dedup$input,unique(.v_gtex$input)) -> gtex.diff

  .v_dedup %>% dplyr::rename("alias" = "Up") %>%
    dplyr::filter(input %in% gtex.diff) %>%
    dplyr::inner_join(total_gene_symbol, by ="alias") %>%
    unique() -> .v_alias.gtex.1

  .v_dedup %>% dplyr::rename("NCBI_sym" = "Up") %>%
    dplyr::filter(input %in% gtex.diff) %>%
    dplyr::inner_join(total_gene_symbol, by ="NCBI_sym") %>%
    unique() -> .v_alias.gtex.2

  rbind(.v_alias.gtex.1,.v_alias.gtex.2) %>% unique() -> .v_alias.ncbi

  # gene id not match in TCGA data
  setdiff(.v_dedup$input,unique(.v_tcga$input)) -> tcga.diff

  .v_dedup %>% dplyr::rename("alias" = "Up") %>%
    dplyr::filter(input %in% tcga.diff) %>%
    dplyr::inner_join(total_gene_symbol, by ="alias") %>%
    unique() -> .v_alias.tcga.1

  .v_dedup %>% dplyr::rename("NCBI_sym" = "Up") %>%
    dplyr::filter(input %in% tcga.diff) %>%
    dplyr::inner_join(total_gene_symbol, by ="NCBI_sym") %>%
    unique() -> .v_alias.tcga.2
  rbind(.v_alias.tcga.1,.v_alias.tcga.2) %>% unique() -> .v_alias.tcga

  .v_tcga %>%
    rbind(.v_gtex) %>%
    rbind(.v_alias.tcga) %>%
    rbind(.v_alias.ncbi) %>%
    unique() -> match_all
  gene_set$match <- match_all %>% dplyr::filter(!is.na(TCGA_sym)) %>% .$TCGA_sym %>% unique()
  gene_set$match.gtex <- match_all %>% dplyr::filter(!is.na(GTEX_sym)) %>% .$GTEX_sym %>% unique()

  .non_match <- match_all %>%
    dplyr::select(input) %>%
    unique() %>%
    tidyr::drop_na() %>%
    dplyr::mutate(match = "TRUE") %>%
    dplyr::right_join(.v_dedup,by = "input") %>%
    dplyr::filter(is.na(match)) %>%
    .$input
  gene_set$non_match <- .non_match
  gene_set$n_match <- length(match_all %>%
                               dplyr::select(input) %>%
                               unique() %>% .$input)
  gene_set$n_non_match <- length(.non_match)
  gene_set$n_total <- gene_set$n_match + gene_set$n_non_match

  if (length(gene_set$match) < 5) {
    error$gene_set <- "Please input at least five valid gene symbol."
    status$trigger <- if (status$trigger == TRUE) FALSE else TRUE
    status$gene_set <- FALSE
  }

  .log <- c(
    glue::glue("{paste0(rep('-', 10), collapse = '')} Notice: Input total gene set number is {length(gene_set$n_total)} {paste0(rep('-', 10), collapse = '')}"),
    glue::glue("{paste0(rep('-', 10), collapse = '')} Notice: Unique gene set number is {length(.v_dedup)} {paste0(rep('-', 10), collapse = '')}"),
    glue::glue("#Total input gene set: {paste0(.v, collapse = ', ')}"),
    glue::glue("#Validated genes: {paste0(gene_set$match, collapse = ', ')}"),
    glue::glue("#Invalidated genes: {paste0(gene_set$non_match, collapse = ', ')}")
  )
  write(x = .log, file = .log_file, append = TRUE)
}
# older version of match gene set
# validate_gene_set <- function(.v, user_dir = user_dir, user_logs = user_logs, total_gene_symbol = total_gene_symbol, status = status, error = error, gene_set = gene_set) {
#   .log_file <- user_logs$gene_set
#
#
#   .v_dedup <- .v[.v != ""] %>% unique() %>% sapply(FUN = tolower, USE.NAMES = FALSE)
#   .v_dedup %in% names(total_gene_symbol) -> .inter
#
#
#   gene_set$match <- total_gene_symbol[.v_dedup[.inter]]
#   gene_set$non_match <- .v[!.inter]    #.v_dedup[!.inter]
#   gene_set$n_match <- length(total_gene_symbol[.v_dedup[.inter]])
#   gene_set$n_non_match <- length(.v_dedup[!.inter])
#   gene_set$n_total <- length(total_gene_symbol[.v_dedup[.inter]]) + length(.v_dedup[!.inter])
#
#   if (length(gene_set$match) < 5) {
#     error$gene_set <- "Please input at least five valid gene symbol."
#     status$trigger <- if (status$trigger == TRUE) FALSE else TRUE
#     status$gene_set <- FALSE
#   }
#
#   .log <- c(
#     glue::glue("{paste0(rep('-', 10), collapse = '')} Notice: Input total gene set number is {length(gene_set$n_total)} {paste0(rep('-', 10), collapse = '')}"),
#     glue::glue("{paste0(rep('-', 10), collapse = '')} Notice: Unique gene set number is {length(.v_dedup)} {paste0(rep('-', 10), collapse = '')}"),
#     glue::glue("#Total input gene set: {paste0(.v, collapse = ', ')}"),
#     glue::glue("#Validated genes: {paste0(gene_set$match, collapse = ', ')}"),
#     glue::glue("#Invalidated genes: {paste0(gene_set$non_match, collapse = ', ')}")
#   )
#   write(x = .log, file = .log_file, append = TRUE)
# }


# cnv bar data prepare ----------------------------------------------------


# threshold cnv -----------------------------------------------------------

fn_get_amplitue_threshold <- function(.x) {
  tibble::tibble(
    a_total = sum(.x > 0) / length(.x),
    d_total = sum(.x < 0) / length(.x),
    a_homo = sum(.x == 2) / length(.x),
    d_homo = sum(.x == -2) / length(.x),
    a_hete = sum(.x == 1) / length(.x),
    d_hete = sum(.x == -1) / length(.x)
  )
}

# get cnv percent ---------------------------------------------------------

fn_get_ad <- function(.d) {
  .d %>%
    unlist(use.name = F) %>%
    fn_get_amplitue_threshold()
}
fn_get_percent <- function(cancer_types, filter_cnv) {
  filter_cnv %>%
    tidyr::nest(-symbol) %>%
    dplyr::mutate(ad = purrr::map(data, .f = fn_get_ad)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(ad) %>%
    tibble::add_column(cancer_types = cancer_types, .before = 1)
}
fn_cnv_percecnt <- function(data) {
  data %>%
    dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>%
    dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    dplyr::select(-cancer_types, -filter_cnv) %>%
    tidyr::unnest(rs) -> gene_list_cnv_per
}


fn_gen_combined_core_atg <- function(cancer_types, filter_cnv, g_list, n) {
  # cancer_types <- "KIRC"
  # filter_cnv <- gene_list_cancer_cnv$filter_cnv[[1]]
  # g_list <-gene_list
  # n=1
  filter_cnv %>%
    dplyr::filter(symbol %in% g_list) %>%
    tidyr::drop_na() %>%
    tidyr::gather(key = barcode, value = gistic, -symbol) %>%
    tidyr::spread(key = symbol, value = gistic) %>%
    dplyr::select(-barcode) -> .d


  n_sample <- nrow(.d)

  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == -n)) %>%
    nrow() -> .del
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == -n)) %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == n)) %>%
    nrow() -> .sub_d

  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == n)) %>%
    nrow() -> .amp
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == n)) %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == -n)) %>%
    nrow() -> .sub_a

  tibble::tibble(del_a = .del / n_sample, del_s = (.del - .sub_d) / n_sample, amp_a = .amp / n_sample, amp_s = (.amp - .sub_a) / n_sample)
}
# extract gene set from TCGA data -----------------------------------------


filter_gene_list <- function(.x, gene_list) {
  .x %>%
    dplyr::filter(symbol %in% gene_list)
  # gene_list %>%
  #   dplyr::select(symbol) %>%
  #   dplyr::left_join(.x, by = "symbol")
}


# get exclusive cnv -------------------------------------------------------
fn_cnv_exclusive <- function(V1, V2, .data, cancer_types) {
  # .data <- filter_cnv
  # V1 <- 'TP53'
  # V2 <- 'EZH2'
  .data %>%
    dplyr::filter(symbol %in% c(V1, V2)) %>%
    tidyr::gather(key = barcode, value = gistic, -symbol) %>%
    tidyr::spread(key = symbol, value = gistic) %>%
    dplyr::select(-barcode) -> .d
  .g_name <- colnames(.d)
  # colnames(.d) <- c("A", "B")
  name <- paste(c(cancer_types, .g_name), collapse = "_")
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. == 0)) %>%
    nrow() -> nn
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. != 0)) %>%
    nrow() -> aa
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. != 0)) %>%
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 0)) -> .d_an

  sum(.d_an %>% dplyr::pull(1) != 0) -> an
  sum(.d_an %>% dplyr::pull(2) != 0) -> na
  c(nn = nn, an = an, na = na, aa = aa) %>%
    cometExactTest::comet_exact_test(mutmatplot = F) -> p_val

  tibble::tibble(te = name, nn = nn, an = an, na = na, aa = aa, p_val = p_val)
}
fn_cnv_mutal_exclusive <- function(cancer_types, filter_cnv, cluster) {
  # cancer_types <- te$cancer_types
  # filter_cnv <- te$filter_cnv[[1]]
  filter_cnv %>%
    dplyr::pull(symbol) %>%
    combn(m = 2) %>%
    t() %>%
    dplyr::as_data_frame() -> .gene_pairs

  .gene_pairs %>%
    multidplyr::partition(cluster = cluster) %>%
    multidplyr::cluster_library("magrittr") %>%
    multidplyr::cluster_assign_value("fn_cnv_exclusive", fn_cnv_exclusive) %>%
    multidplyr::cluster_assign_value("filter_cnv", filter_cnv) %>%
    multidplyr::cluster_assign_value("cancer_types", cancer_types) %>%
    dplyr::mutate(rs = purrr::map2(V1, V2, .f = fn_cnv_exclusive, .data = filter_cnv, cancer_types = cancer_types)) %>%
    dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    dplyr::select(-PARTITION_ID) %>%
    dplyr::select(rs) %>%
    tidyr::unnest() %>%
    tidyr::separate(col = te, into = c("cancer_types", "g1", "g2")) -> .gene_pairs_pval

  .gene_pairs_pval %>%
    dplyr::mutate(fdr = p.adjust(p_val, method = "fdr"))
}


# rppa line contact faction -----------------------------------------------

get_rppa_text <- function(data) {
  data %>%
    dplyr::pull(symbol) %>%
    unique() -> gene.text
  data %>%
    dplyr::pull(cancer_types) %>%
    unique() -> cancer.text
  data %>%
    dplyr::pull(pathway) %>%
    unique() -> pathway.text

  c.text <- data.frame(x = 0.5, y = 1, text = "test", type = "test")
  g.l <- data$symbol %>% unique() %>% length()
  c.l <- data$cancer_types %>% unique() %>% length()
  p.l <- data$pathway %>% unique() %>% length()

  # condition 1: gene is more -----------------------------------------------

  if (g.l >= c.l & g.l >= p.l) {
    for (i in 1:length(gene.text)) {
      data.frame(x = 4, y = 2 * i - 1, text = gene.text[i], type = "gene") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }
    c.text %>%
      dplyr::filter(type == "gene") %>%
      dplyr::select(y) %>%
      max() -> g.m

    g.m / c.l -> c.i
    for (i in 1:length(cancer.text)) {
      data.frame(x = 1, y = i * c.i - 1, text = cancer.text[i], type = "cancer") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }

    g.m / 10 -> p.i
    for (i in 1:length(pathway.text)) {
      data.frame(x = 7, y = i * p.i - 1, text = pathway.text[i], type = "pathway") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }
  }

  # condition 2: cancer is more ---------------------------------------------

  if (c.l >= g.l & c.l >= p.l) {
    for (i in 1:length(cancer.text)) {
      data.frame(x = 1, y = 2 * i - 1, text = cancer.text[i], type = "cancer") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }
    c.text %>%
      dplyr::filter(type == "cancer") %>%
      dplyr::select(y) %>%
      max() -> c.m

    c.m / g.l -> g.i
    for (i in 1:length(gene.text)) {
      data.frame(x = 4, y = i * g.i - 1, text = gene.text[i], type = "gene") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }

    c.m / 10 -> p.i
    for (i in 1:length(pathway.text)) {
      data.frame(x = 7, y = i * p.i - 1, text = pathway.text[i], type = "pathway") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }
  }

  # condition 3: pathway is more --------------------------------------------

  if (p.l >= c.l & p.l >= g.l) {
    for (i in 1:length(pathway.text)) {
      data.frame(x = 7, y = i * 2 - 1, text = pathway.text[i], type = "pathway") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }

    c.text %>%
      dplyr::filter(type == "pathway") %>%
      dplyr::select(y) %>%
      max() -> p.m

    p.m / c.l -> c.i
    for (i in 1:length(cancer.text)) {
      data.frame(x = 1, y = i * c.i - 1, text = cancer.text[i], type = "cancer") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }

    p.m / g.l -> g.i
    for (i in 1:length(gene.text)) {
      data.frame(x = 4, y = i * g.i - 1, text = gene.text[i], type = "gene") -> tmp.text
      rbind(c.text, tmp.text) -> c.text
    }
  }
  return(c.text[-1, ])
}

get_rppa_seg <- function(data,cancer_text) {
  # name <- c("x1","y1","x2","y2","Cancer","Regulation")
  # print(n)
  data[1,1] %>% as.character() -> cancer
  data[1,2] %>% as.character() -> gene
  data[1,3] %>% as.character() -> pathway
  data[1,4] %>% as.numeric() -> diff
  if (diff > 0) {
    line_type <- "Activate"
  } else {
    line_type <- "Inhibit"
  }
  cancer_text %>%
    dplyr::filter(text %in% gene) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(x = x - 0.5) -> g1.pos
  cancer_text %>%
    dplyr::filter(text %in% gene) %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(x = x + 0.5) -> g2.pos

  cancer_text %>%
    dplyr::filter(text %in% cancer) %>%
    dplyr::select(x, y) -> c.pos
  cancer_text %>%
    dplyr::filter(text %in% pathway) %>%
    dplyr::select(x, y) -> p.pos
  .d_seq_tmp1 <- data.frame(x1 = c.pos$x, y1 = c.pos$y, x2 = g1.pos$x, y2 = g1.pos$y, Cancer = cancer, Regulation = "Activate")
  .d_seq_tmp2 <- data.frame(x1 = g2.pos$x, y1 = g2.pos$y, x2 = p.pos$x, y2 = p.pos$y, Cancer = cancer, Regulation = line_type)
  rbind(.d_seq_tmp1,.d_seq_tmp2) -> .d_seg
  .d_seg$Cancer <- .d_seg$Cancer %>% as.character()
  .d_seg$Regulation <- .d_seg$Regulation %>% as.character()
  tibble::as_tibble(.d_seg)
}

get_rppa_seg1 <- function(cancer_text, data) {
  .d_seg <- data.frame(x1 = 0, y1 = 0, x2 = 0, y2 = 0, Cancer = "test", Regulation = "test")
  nrow(data) -> n

  for (i in 1:n) {
    data[i, 1] -> cancer
    data[i, 2] -> gene
    data[i, 3] -> pathway
    data[i, 4] -> diff
    if (diff > 0) {
      line_type <- "Activate"
    } else {
      line_type <- "Inhibit"
    }
    cancer_text %>%
      dplyr::filter(text %in% gene) %>%
      dplyr::select(x, y) %>%
      dplyr::mutate(x = x - 0.5) -> g1.pos
    cancer_text %>%
      dplyr::filter(text %in% gene) %>%
      dplyr::select(x, y) %>%
      dplyr::mutate(x = x + 0.5) -> g2.pos

    cancer_text %>%
      dplyr::filter(text %in% cancer) %>%
      dplyr::select(x, y) -> c.pos
    cancer_text %>%
      dplyr::filter(text %in% pathway) %>%
      dplyr::select(x, y) -> p.pos
    .d_seq_tmp1 <- data.frame(x1 = c.pos$x, y1 = c.pos$y, x2 = g1.pos$x, y2 = g1.pos$y, Cancer = cancer$cancer_types, Regulation = "Activate")
    .d_seq_tmp2 <- data.frame(x1 = g2.pos$x, y1 = g2.pos$y, x2 = p.pos$x, y2 = p.pos$y, Cancer = cancer$cancer_types, Regulation = line_type)
    rbind(.d_seg, .d_seq_tmp1) %>% rbind(.d_seq_tmp2) -> .d_seg
  }
  .d_seg[-1, ] %>%
    unique() -> .d_seg
  return(.d_seg)
}


# maftools subsetMaf edit to suit our addition. ---------------------------
my_subsetMaf <- function (maf, tsb = NULL, genes = NULL, fields = NULL, cancer = NULL,
                          mafObj = FALSE, includeSyn = TRUE, isTCGA = FALSE)
{
  maf.silent <- maf@maf.silent
  maf.dat <- maf@data
  maf.anno <- maf@clinical.data
  if (!is.null(tsb)) {
    if (isTCGA) {
      tsb = substr(x = tsb, start = 1, stop = 12)
    }
    maf.dat = maf.dat[Tumor_Sample_Barcode %in% tsb, ]
    maf.silent = maf.silent[Tumor_Sample_Barcode %in% tsb,
                            ]
  }
  if (!is.null(genes)) {
    maf.dat = maf.dat[Hugo_Symbol %in% genes, ]
    maf.silent = maf.silent[Hugo_Symbol %in% genes, ]
  }
  if (!is.null(cancer)) {
    # maf.dat = maf.dat[eval(parse(text = query))]
    maf.dat = maf.dat[Cancer_Types %in% cancer,]
    # maf.silent = maf.silent[eval(parse(text = query))]
    maf.silent = maf.silent[Cancer_Types %in% cancer,]
  }
  default.fields = c("Hugo_Symbol", "Chromosome", "Start_Position",
                     "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                     "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")
  if (!is.null(fields)) {
    default.fields = unique(c(default.fields, fields))
    if (length(default.fields[!default.fields %in% colnames(maf.dat)]) >
        0) {
      message("Missing fields. Ignoring them.. ")
      print(default.fields[!default.fields %in% colnames(maf.dat)])
      default.fields = default.fields[default.fields %in%
                                        colnames(maf.dat)]
    }
    maf.dat = maf.dat[, default.fields, with = FALSE]
    maf.silent = maf.silent[, default.fields, with = FALSE]
  }
  if (mafObj) {
    maf.silent = droplevels.data.frame(maf.silent)
    maf.dat = droplevels.data.frame(maf.dat)
    maf.anno = droplevels.data.frame(maf.anno)
    mafSummary = my_summarizeMaf(maf.dat, chatty = FALSE, anno = maf.anno)
    m = my_MAF(data = maf.dat, variants.per.sample = mafSummary$variants.per.sample,
            variant.type.summary = mafSummary$variant.type.summary,
            variant.classification.summary = mafSummary$variant.classification.summary,
            gene.summary = mafSummary$gene.summary, summary = mafSummary$summary,
            maf.silent = maf.silent, clinical.data = mafSummary$sample.anno)
    return(m)
  }
  else {
    if (includeSyn) {
      return(rbind(maf.dat, maf.silent, use.names = TRUE,
                   fill = TRUE))
    }
    else {
      return(maf.dat)
    }
  }
}

my_summarizeMaf = function(maf, anno = NULL, chatty = TRUE){

  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[!Variant_Type %in% 'CNV', NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]

    if(chatty){
      if(length(NCBI_Build) > 1){
        message('NOTE: Mutiple reference builds found!')
        NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
        message(NCBI_Build)
      }
    }
  }else{
    NCBI_Build = NA
  }

  if('Center' %in% colnames(maf)){
    Center = unique(maf[!Variant_Type %in% 'CNV', Center])
    #Center = Center[is.na(Center)]
    if(length(Center) > 1){
      Center = do.call(paste, c(as.list(Center), sep=";"))
      if(chatty){
        message('Mutiple centers found.')
        print(Center)
      }
    }
  }else{
    Center = NA
  }

  #nGenes
  nGenes = length(unique(maf[,Hugo_Symbol]))

  #Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
  flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
            "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
            "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]

  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')

  if(any(colnames(vc.cast) %in% c('Amp', 'Del'))){
    vc.cast.cnv = vc.cast[,c('Tumor_Sample_Barcode', colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')]), with =FALSE]
    vc.cast.cnv$CNV_total = rowSums(vc.cast.cnv[,2:ncol(vc.cast.cnv)], na.rm = TRUE)

    vc.cast = vc.cast[,!colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')], with =FALSE]
    vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]

    vc.cast = merge(vc.cast, vc.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]

    vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
    vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))

  }else{
    vc.cast = vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])][order(total, decreasing = TRUE)]

    vc.mean = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean)))), 3)
    vc.median = round(as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median)))), 3)
  }

  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)

  if(any(colnames(vt.cast) %in% c('CNV'))){
    vt.cast.cnv = vt.cast[,c('Tumor_Sample_Barcode', colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')]), with =FALSE]

    vt.cast = vt.cast[,!colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')], with =FALSE]
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]

    vt.cast = merge(vt.cast, vt.cast.cnv, by = 'Tumor_Sample_Barcode', all = TRUE)[order(total, CNV, decreasing = TRUE)]
  }else{
    vt.cast = vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])][order(total, decreasing = TRUE)]
  }

  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  #----
  if(any(colnames(hs.cast) %in% c('Amp', 'Del'))){
    hs.cast.cnv = hs.cast[,c('Hugo_Symbol', colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')]), with = FALSE]
    hs.cast.cnv$CNV_total = rowSums(x = hs.cast.cnv[,2:ncol(hs.cast.cnv), with = FALSE], na.rm = TRUE)

    hs.cast = hs.cast[,!colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')], with = FALSE]
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE], na.rm = TRUE)]

    hs.cast = merge(hs.cast, hs.cast.cnv, by = 'Hugo_Symbol', all = TRUE)[order(total, CNV_total, decreasing = TRUE)]
  }else{
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
    hs.cast = hs.cast[order(total, decreasing = TRUE)]
  }
  #----

  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = maf[, .(AlteredSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  numAlteredSamples = merge(numMutatedSamples, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)
  #Merge and sort
  hs.cast = merge(hs.cast, numAlteredSamples, by = 'Hugo_Symbol', all = TRUE)[order(MutatedSamples, total, decreasing = TRUE)]
  #Replace NAs with 0
  hs.cast$AlteredSamples = ifelse(test = is.na(x = hs.cast$AlteredSamples), yes = 0, no = hs.cast$AlteredSamples)
  hs.cast$MutatedSamples = ifelse(test = is.na(x = hs.cast$MutatedSamples), yes = 0, no = hs.cast$MutatedSamples)
  #Make a summarized table
  summary = data.table::data.table(ID = c('NCBI_Build', 'Center','Samples', 'nGenes',colnames(vc.cast)[2:ncol(vc.cast)]),
                                   summary = c(NCBI_Build, Center, nrow(vc.cast), nGenes, colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))
  summary[,Mean := vc.mean]
  summary[,Median := vc.median]

  if(chatty){
    print(summary)

    message("Gene Summary..")
    print(hs.cast)
  }

  #Check for flags.
  if(nrow(hs.cast) > 10){
    topten = hs.cast[1:10, Hugo_Symbol]
    topten = topten[topten %in% flags]
    if(chatty){
      if(length(topten) > 0){
        message('NOTE: Possible FLAGS among top ten genes:')
        print(topten)
      }
    }
  }


  if(chatty){
    message("Checking clinical data..")
  }

  if(is.null(anno)){
    if(chatty){
      message("NOTE: Missing clinical data! It is strongly recommended to provide clinical data associated with samples if available.")
    }
    sample.anno = tsb[,.(Tumor_Sample_Barcode)]
  }else if(is.data.frame(x = anno)){
    sample.anno  = data.table::setDT(anno)
    if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
      message(paste0('Available fields in provided annotations..'))
      print(colnames(sample.anno))
      stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.'))
    }
  }else{
    if(file.exists(anno)){
      sample.anno = data.table::fread(anno, stringsAsFactors = FALSE)
      if(!'Tumor_Sample_Barcode' %in% colnames(sample.anno)){
        message(paste0('Available fields in ', basename(anno), '..'))
        print(colnames(sample.anno))
        stop(paste0('Tumor_Sample_Barcode column not found in provided clinical data. Rename column name containing sample names to Tumor_Sample_Barcode if necessary.'))
      }
    }
  }

  #clean up annotation data
  colnames(sample.anno) = gsub(pattern = ' ', replacement = '_', x = colnames(sample.anno), fixed = TRUE) #replace spaces in column names for annotation data
  sample.anno = as.data.frame(apply(sample.anno, 2, function(y) trimws(y))) #remove trailing whitespaces
  sample.anno[sample.anno == ""] = NA #Replace blanks with NA
  sample.anno = as.data.frame(apply(sample.anno, 2, function(y) gsub(pattern = " ", replacement = "_", x = y))) #replace spaces with _
  data.table::setDT(x = sample.anno)
  colnames(sample.anno)[1] = c("Tumor_Sample_Barcode")

  maf.tsbs = levels(tsb[,Tumor_Sample_Barcode])
  sample.anno = sample.anno[Tumor_Sample_Barcode %in% maf.tsbs][!duplicated(Tumor_Sample_Barcode)]
  anno.tsbs = sample.anno[,Tumor_Sample_Barcode]

  if(!length(maf.tsbs[!maf.tsbs %in% anno.tsbs]) == 0){
    if(chatty){
      message('Annotation missing for below samples in MAF')
      print(maf.tsbs[!maf.tsbs %in% anno.tsbs])
    }
  }

  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary, sample.anno = sample.anno))
}

## MAF object
my_MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                          variant.classification.summary = 'data.table', gene.summary = 'data.table',
                                          summary = 'data.table', maf.silent = 'data.table', clinical.data = 'data.table'))



# Loading screen ----------------------------------------------------------

loading_screen <- function(){
  shinyjs::hide(id = "loading-content", anim = TRUE, animType = "fade")
  shinyjs::show("app-content")
}
# generate cnv result out ui -------------------------------------------------------

output$ui_cnv_result <- shiny::renderUI({
  fn_cnv_result(selected_analysis$cnv)
})



# analysis start ----------------------------------------------------------
cnv_analysis <- eventReactive(
  {
    status$analysis
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      if (selected_analysis$cnv == TRUE) {
        # Cancer types value box selection ----------------------------------------

        callModule(module = cancerTypesSelect, id = "cnv", .sctps = intersect(selected_ctyps(), tcga_data))

        # Check box ---------------------------------------------------------------

        callModule(module = selectAndAnalysis, id = "cnv", .id = "cnv")
        .msg <- c("NOTICE: ")
        
        # load data----
        load_data_cnv()
        cnv_gene_new <- gene_set$match
        cnv_cancer_new <- selected_ctyps()

        # cancer overlap ----
        cancer_in_tcga_data_cnv <- intersect(selected_ctyps(),tcga_data)
        # print(cnv_cancer_diff)
        # print(cnv_gene_diff)

        if (length(gene_set$match) != 0) {
          if (!setequal(cnv_gene_new, cnv_gene_old) | !setequal(cnv_cancer_new, cnv_cancer_old)) {
            cnv_gene_old <- gene_set$match
            cnv_cancer_old <- selected_ctyps()

            # cnv percent plot ------------------------------------------------------------

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} start cnv pie percent data processing@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # get cancer type cnv ----

            cnv %>%
              dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_set$match)) %>%
              dplyr::select(-cnv) %>%
              dplyr::filter(cancer_types %in% selected_ctyps()) -> gene_list_cancer_cnv

            # get data for plot ----
            gene_list_cancer_cnv %>%
              tidyr::unnest() %>%
              tidyr::drop_na() -> cnv_plot_ready
            if (nrow(cnv_plot_ready) > 0) {

              # cancer rank ----
              cnv_plot_ready %>%
                dplyr::group_by(cancer_types) %>%
                dplyr::summarise(v = sum(a_total - d_total)) %>%
                dplyr::arrange(dplyr::desc(v)) -> cnv_cancer_rank

              # gene rank ----
              cnv_plot_ready %>%
                dplyr::group_by(symbol) %>%
                dplyr::summarise(v = sum(a_total - d_total)) %>%
                dplyr::arrange(v) -> cnv_gene_rank

              # plot generate ----

              # pie plot ----
              cnv_plot_ready %>%
                dplyr::select(-a_total, -d_total) %>%
                tidyr::gather(key = type, value = per, -c(cancer_types, symbol)) %>%
                dplyr::mutate(
                  symbol = factor(x = symbol, levels = cnv_gene_rank$symbol),
                  cancer_types = factor(x = cancer_types, levels = cnv_cancer_rank$cancer_types)
                ) -> pie_plot_ready

              # cnvpie_getheight <- function(cn) {
              #   if (cn <= 5) {
              #     return(0.27)
              #   }
              #   if (cn > 5 && cn <= 20) {
              #     return(0.25 - (cn - 5) * 0.01)
              #   } else {
              #     return(0.15)
              #   }
              # }
              cnv_pie_gn <- pie_plot_ready$symbol %>% unique() %>% length()
              cnv_pie_cn <- pie_plot_ready$cancer_types %>% unique() %>% length()
              if(cnv_pie_cn < 7){
                cnv_pie_width <- 2
                cnv_pie_height <- 0.2*cnv_pie_gn} 
              if(cnv_pie_cn >= 7 && cnv_pie_cn<15){
                cnv_pie_width <- cnv_pie_cn * 0.27
                cnv_pie_height <- cnv_pie_gn * 0.27
              }
              if(cnv_pie_cn >= 15){
                cnv_pie_width <- cnv_pie_cn * 0.27
                cnv_pie_height <- cnv_pie_gn * 0.27
              }
              print(paste0("cnv_pie_width:",cnv_pie_width))
              print(paste0("cnv_pie_height:",cnv_pie_height))
              # cnv_pie_h <- cnvpie_getheight(cn = pie_plot_ready$cancer_types %>% unique() %>% length())
              # cnv_pie_height <- pie_plot_ready$symbol %>% unique() %>% length() * 0.27 #cnv_pie_h
              
              # if (cnv_pie_height > 15) {
              #   cnv_pie_height <- 15
              # }
              # if (cnv_pie_height < 3) {
              #   cnv_pie_height <- 3
              # }
              print(glue::glue("{paste0(rep('-', 10), collapse = '')} End cnv percent data processing@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

              print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start gernerate cnv pie profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

              callModule(
                piePlot, "cnv_pie", data = pie_plot_ready, y = "per",
                fill = "type", facet_grid = "symbol ~ cancer_types",
                outfile = file.path(user_dir, "pngs", paste(user_id, "-CNV_pie_profile.png", sep = "")), height = cnv_pie_height,
                width = cnv_pie_width,
                status_monitor = "analysis", status, downloadname = "cnv_percent_profile_figure"
              )
              .msg_cnv_pie <- NULL
            } else {
              .msg_cnv_pie <- paste(glue::glue("No significant [CNV Pie distribution] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types {paste0(cancer_in_tcga_data_cnv,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
              output[["cnv_pie-plot"]] <- callModule(white_plot, "cnv_pie", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_1.png", sep = "")))
            }

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} End gernerate cnv pie profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # homo cnv plot ----
            print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start gernerate homo cnv profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            cnv_homo_plot_ready <- cnv_plot_ready %>%
              dplyr::select(cancer_types, symbol, a_homo, d_homo) %>%
              tidyr::gather(key = type, value = per, -cancer_types, -symbol) %>%
              dplyr::mutate(effect = plyr::revalue(type, replace = c("a_homo" = "Homozygous Amplification", "d_homo" = "Homozygous Deletion"))) %>%
              dplyr::mutate(color = plyr::revalue(type, replace = c("a_homo" = "brown4", "d_homo" = "aquamarine4")))
            if (nrow(cnv_homo_plot_ready) > 0) {
              callModule(
                cnv_pointPlot, "cnv_homo", data = cnv_homo_plot_ready, cancer = "cancer_types",
                gene = "symbol", size = "per", color = "color", sizename = "Homo CNV%",
                colorname = "SCNA Type", wrap = "~ effect", status_monitor = "analysis", status,
                downloadname = "cnv_homo_figure"
              )
              .msg_cnv_homo <- NULL
            } else {
              .msg_cnv_homo <- paste(glue::glue("No significant [Homo CNV profile] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(cancer_in_tcga_data_cnv,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
              output[["cnv_homo-plot"]] <- renderPlot({
                NULL
              })
            }


            print(glue::glue("{paste0(rep('-', 10), collapse = '')} End gernerate homo cnv profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # hete cnv plot ----
            print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start gernerate hete cnv profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            cnv_hete_plot_ready <- cnv_plot_ready %>%
              dplyr::select(cancer_types, symbol, a_hete, d_hete) %>%
              tidyr::gather(key = type, value = per, -cancer_types, -symbol) %>%
              dplyr::mutate(effect = plyr::revalue(type, replace = c("a_hete" = "Heterozygous Amplification", "d_hete" = "Heterozygous Deletion"))) %>%
              dplyr::mutate(color = plyr::revalue(type, replace = c("a_hete" = "brown1", "d_hete" = "aquamarine3")))

            if (nrow(cnv_hete_plot_ready) > 0) {
              callModule(
                cnv_pointPlot, "cnv_hete", data = cnv_hete_plot_ready, cancer = "cancer_types",
                gene = "symbol", size = "per", color = "color", sizename = "Hete CNV%",
                colorname = "SCNA Type", wrap = "~ effect", status_monitor = "analysis", status,
                downloadname = "cnv_hete_figure"
              )
              .msg_cnv_hete <- NULL
            } else {
              .msg_cnv_hete <- paste(glue::glue("No significant [Hete CNV profile] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(cancer_in_tcga_data_cnv,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
              output[["cnv_hete-plot"]] <- renderPlot({
                NULL
              })
            }

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} End gernerate hete cnv profile plot@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # cnv bar plot ------------------------------------------------------------

            # print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start processing cnv overall percent data@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
            # cnv_raw %>%
            #   dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_set$match)) %>%
            #   dplyr::select(-cnv) %>%
            #   dplyr::filter(cancer_types %in% selected_ctyps()) -> gene_list_cancer_cnv_raw

            # bar stack plot ----
            # gene_list_cancer_cnv_raw %>%
            #   dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_gen_combined_core_atg, g_list = gene_set$match, n = 1)) %>%
            #   dplyr::select(-filter_cnv) %>%
            #   tidyr::unnest(rs) %>%
            #   dplyr::mutate(del_a = -del_a) %>%
            #   dplyr::mutate(del_s = -del_s) %>%
            #   tidyr::gather(key = type, value = per, -cancer_types) %>%
            #   dplyr::mutate(cnv_type = "Hete CNV") -> cnv_hete_bar_plot_ready
            # 
            # gene_list_cancer_cnv_raw %>%
            #   dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_gen_combined_core_atg, g_list = gene_set$match, n = 2)) %>%
            #   dplyr::select(-filter_cnv) %>%
            #   tidyr::unnest(rs) %>%
            #   dplyr::mutate(del_a = -del_a) %>%
            #   dplyr::mutate(del_s = -del_s) %>%
            #   tidyr::gather(key = type, value = per, -cancer_types) %>%
            #   dplyr::mutate(cnv_type = "Homo CNV") -> cnv_homo_bar_plot_ready
            # 
            # rbind(cnv_hete_bar_plot_ready, cnv_homo_bar_plot_ready) -> cnv_bar_plot_ready
            # 
            # print(glue::glue("{paste0(rep('-', 10), collapse = '')} End processing cnv overall percent data@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
            # if (nrow(cnv_bar_plot_ready) > 0) {
            #   callModule(cnvbarPlot, "cnv_bar", data = cnv_bar_plot_ready, x = "cancer_types", y = "per", fill = "type", status_monitor = "analysis", status, downloadname = "cnv_hete_figure")
            #   .msg_cnv_bar <- NULL
            # } else {
            #   .msg_cnv_bar <- paste(glue::glue("No significant [Overall CNV frenquency] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(cancer_in_tcga_data_cnv,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
            #   output[["cnv_bar-plot"]] <- renderPlot({
            #     NULL
            #   })
            # }

            # cnv cor to expressin ----------------------------------------------------

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start cnv cor to expression ploting@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
            cnv_cor %>%
              dplyr::mutate(spm = purrr::map(spm, filter_gene_list, gene_list = gene_set$match)) %>%
              dplyr::filter(cancer_types %in% selected_ctyps()) %>%
              tidyr::unnest() -> gene_list_cancer_cnv_cor

            if (nrow(gene_list_cancer_cnv_cor) > 0) {
              # cnv to expression plot  ----
              gene_list_cancer_cnv_cor %>%
                dplyr::group_by(symbol) %>%
                dplyr::summarise(rank = sum(spm)) %>%
                dplyr::arrange(rank) -> gene_rank.cnvcor

              gene_list_cancer_cnv_cor %>%
                dplyr::group_by(cancer_types) %>%
                dplyr::summarise(rank = sum(spm)) %>%
                dplyr::arrange(rank) -> cancer_rank.cnvcor

              callModule(methy_diff_pointPlot, "cnv_exp", data = gene_list_cancer_cnv_cor, cancer = "cancer_types", gene = "symbol", size = "logfdr", color = "spm", cancer_rank = cancer_rank.cnvcor, gene_rank = gene_rank.cnvcor, sizename = "-Log10(FDR)", colorname = "Pearson Correlation", title = "Pearson Correlation between CNV and mRNA RSEM.", status_monitor = "analysis", status, downloadname = "cnv_correlate_to_expr")
              .msg_cnv_exp <- NULL
            } else {
              .msg_cnv_exp <- paste(glue::glue("No significant [CNV to Expression] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(cancer_in_tcga_data_cnv,collapse=', ')}. Please try more cancers or more genes."), sep = " ")

              output[["cnv_exp-plot"]] <- renderPlot({
                NULL
              })
            }

            # infomation UI for each part
            output[["cnv_pie-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_cnv_pie, style = "color:#CD3700")
              )
            })
            output[["cnv_hete-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_cnv_hete, style = "color:#CD3700")
              )
            })
            output[["cnv_homo-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_cnv_homo, style = "color:#CD3700")
              )
            })
            output[["cnv_bar-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_cnv_bar, style = "color:#CD3700")
              )
            })
            output[["cnv_exp-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_cnv_exp, style = "color:#CD3700")
              )
            })

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} End cnv cor to expression ploting@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            .msg <- paste(.msg, glue::glue("Since we just show significant results, so a small size of gene and cancer set may cause no significant result in some plots. If it happens, try more genes and cancer types."), sep = " ")

            # alert for information
            shinyBS::createAlert(
              session = session, anchorId = "cnv-no_gene_set", title = "Information", style = "info",
              content = .msg, append = FALSE
            )
            status$cnv_submit <- FALSE
            shinyjs::enable("cnv-submit")
          } else {
            print("not do")
            cnv_gene_old <- gene_set$match
            cnv_cancer_old <- selected_ctyps()
          }
          print(cnv_gene_old)
          print(cnv_cancer_old)
        } else {
          shinyBS::createAlert(
            session = session, anchorId = "snv-no_gene_set", title = "Oops",
            content = "No input gene set! Please go to Welcome page to input gene set.", style = "danger", append = FALSE
          )
        }
      }
    }
  }
)
# monitor ---------------------------------------------------------------------
observe(cnv_analysis())
# sourced by 'server.R'
# save as 'tcga_meth_server.R'
# server elements 'tcga_meth' sub tab of 'tcga' tab

source(file.path(config$wd, "functions", "tcga_meth_function.R"))

#  get cancer type --------------------------------------------------------

# meth_cancer_type <- callModule(cancerType, "meth")




# analysis core -----------------------------------------------------------
# generate meth result out ui -------------------------------------------------------

output$ui_meth_result <- shiny::renderUI({
  fn_meth_result(selected_analysis$meth)
})

# analysis core -----------------------------------------------------------
meth_analysis <- eventReactive(
  {
    status$analysis == TRUE
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      if (selected_analysis$meth == TRUE) {
        # Cancer types value box selection ----------------------------------------

        callModule(module = cancerTypesSelect, id = "meth", .sctps = intersect(selected_ctyps(), tcga_data))
        # Check box ---------------------------------------------------------------

        callModule(module = selectAndAnalysis, id = "meth", .id = "meth")
        load_data_meth()
        if (length(gene_set$match) != 0) {
          .msg <- c("NOTICE: ")
          # load data----
          load_data_meth()

          # cancer overlap
          cancer_in_tcga_data_meth <- intersect(selected_ctyps(),tcga_data)

          print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start methy part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

          # get gene set meth ----------------------------------------
          meth_diff %>%
            dplyr::mutate(filter_methyDiff = purrr::map(methy_comparison, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-methy_comparison) %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) -> gene_list_meth_diff

          meth_survival %>%
            dplyr::mutate(filter_SurDiff = purrr::map(diff_pval, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-diff_pval) -> gene_list_meth_sur

          meth_cor %>%
            dplyr::mutate(filter_cor = purrr::map(spm, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-spm) -> gene_list_meth_cor

          # get cancer type meth ----------------------------------------
          if (nrow(gene_list_meth_diff) > 0) {
            gene_list_meth_diff %>%
              tidyr::unnest() %>%
              tidyr::drop_na() -> gene_list_cancer_methdiff


            # ploting -----------------------------------------------------------------
            # meth diff point ----
            gene_list_cancer_methdiff %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(diff)) %>%
              dplyr::arrange(rank) -> gene_rank.methdiff
            gene_list_cancer_methdiff %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(diff)) %>%
              dplyr::arrange(rank) -> cancer_rank.methdiff

            callModule(methy_diff_pointPlot, "meth_diff", data = gene_list_cancer_methdiff, cancer = "cancer_types", gene = "symbol", size = "fdr", color = "diff", cancer_rank = cancer_rank.methdiff, gene_rank = gene_rank.methdiff, sizename = "-Log10(FDR)", colorname = "Methylation Diff (T - N)", title = "Methylation difference between tumor and normal samples.", status_monitor = "analysis", status, downloadname="Differential_methylation")
            .msg_meth_diff <- NULL
          } else {
            .msg_meth_diff <- paste(glue::glue("The [Differential Methylation] analysis based on paired sample in each cancer types.
In this analysis, only {nrow(meth_diff)} cancer types have paired samples. They are {paste0(meth_diff$cancer_types, collapse = ', ')}."), sep = " ")
            output[["meth_diff-plot"]] <- renderPlot({NULL})
          }

          # meth survival point ----
          gene_list_meth_sur %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) %>%
            tidyr::unnest() %>%
            tidyr::drop_na() -> gene_list_cancer_methsur
          if (nrow(gene_list_cancer_methsur) > 0) {
            gene_list_cancer_methsur %>%
              dplyr::mutate(a = ifelse(Hyper_worse == "Low", -1, 1)) %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(a)) %>%
              dplyr::arrange(rank) -> gene_rank.methsur

            gene_list_cancer_methsur %>%
              dplyr::mutate(a = ifelse(Hyper_worse == "Low", -1, 1)) %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(a)) %>%
              dplyr::arrange(rank) -> cancer_rank.methsur

            callModule(snv_sur_pointPlot, "meth_survival", data = gene_list_cancer_methsur, cancer = "cancer_types", gene = "symbol", size = "log10logrankP", color = "Hyper_worse", cancer_rank = cancer_rank.methsur, gene_rank = gene_rank.methsur, sizename = "logRank Pvalue", colorname = "Effect of HyperMethy on survival risk", title = "Overall survival difference between hypermethylation and hypomethylation.", status_monitor = "analysis", status, downloadname="Methylation_survival")
            .msg_meth_survival <- NULL
          } else {
            .msg_meth_survival <- paste(.msg, glue::glue("No significant [Methylation Survival] result of gene: {paste0(gene_set$match, collapse = ', ')} in your selected cancer type: {paste0(cancer_in_tcga_data_meth,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
            output[["meth_survival-plot"]] <- renderPlot({
              NULL
            })
          }

          # meth correlate to expression point ----
          gene_list_meth_cor %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) %>%
            tidyr::unnest() %>%
            tidyr::drop_na() -> gene_list_cancer_methcor

          if (nrow(gene_list_cancer_methcor) > 0) {
            gene_list_cancer_methcor %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(spm)) %>%
              dplyr::arrange(rank) -> gene_rank.methcor

            gene_list_cancer_methcor %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(spm)) %>%
              dplyr::arrange(rank) -> cancer_rank.methcor

            callModule(methy_diff_pointPlot, "meth_exp", data = gene_list_cancer_methcor, cancer = "cancer_types", gene = "symbol", size = "logfdr", color = "spm", cancer_rank = cancer_rank.methcor, gene_rank = gene_rank.methcor, sizename = "-Log10(P.value)", colorname = "Spearman Correlation Coefficient", title = "Spearman Correlation Coefficient of methylation and gene expression.", status_monitor = "analysis", status, downloadname="Methylation_affect_exp")
            .msg_meth_exp <- NULL
          } else {
            .msg_meth_exp <- paste(.msg, glue::glue("No significant [Methylation to Expression] result of gene: {paste0(gene_set$match, collapse = ', ')} in your selected cancer type: {paste0(cancer_in_tcga_data_meth,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
            output[["meth_exp-plot"]] <- renderPlot({
              NULL
            })
          }
          print(glue::glue("{paste0(rep('-', 10), collapse = '')} End methy part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

          # infomation UI for each part --------------------------
          output[["meth_diff-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_diff,style= "color:#CD3700")
            )
          })
          output[["meth_survival-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_survival,style= "color:#CD3700")
            )
          })
          output[["meth_exp-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_exp,style= "color:#CD3700")
            )
          })

          .msg <- paste(.msg, glue::glue("Since we just show significant results, so a small size of gene and cancer set may cause no significant result in some plots, if it happens, try more genes and cancer types."), sep = " ")
          # alert for information
          shinyBS::createAlert(
            session = session, anchorId = "meth-no_gene_set", title = "Information", style = "info",
            content = .msg, append = FALSE
          )
        } else {
          shinyBS::createAlert(
            session = session, anchorId = "meth-no_gene_set", title = "Oops",
            content = "No input gene set! Please go to Welcome page to input gene set.", style = "danger", append = FALSE
          )
        }
      }
    }
  }
)


# monitors -------------------------------------------------------
observe(meth_analysis())
# sourced by 'server.R'
# save as 'tcga_meth_server.R'
# server elements 'tcga_meth' sub tab of 'tcga' tab

source(file.path(config$wd, "functions", "tcga_meth_function.R"))

#  get cancer type --------------------------------------------------------

# meth_cancer_type <- callModule(cancerType, "meth")




# analysis core -----------------------------------------------------------
# generate meth result out ui -------------------------------------------------------

output$ui_meth_result <- shiny::renderUI({
  fn_meth_result(selected_analysis$meth)
})

# analysis core -----------------------------------------------------------
meth_analysis <- eventReactive(
  {
    status$analysis == TRUE
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      if (selected_analysis$meth == TRUE) {
        # Cancer types value box selection ----------------------------------------

        callModule(module = cancerTypesSelect, id = "meth", .sctps = intersect(selected_ctyps(), tcga_data))
        # Check box ---------------------------------------------------------------

        callModule(module = selectAndAnalysis, id = "meth", .id = "meth")
        load_data_meth()
        if (length(gene_set$match) != 0) {
          .msg <- c("NOTICE: ")
          # load data----
          load_data_meth()

          # cancer overlap
          cancer_in_tcga_data_meth <- intersect(selected_ctyps(),tcga_data)

          print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start methy part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

          # get gene set meth ----------------------------------------
          meth_diff %>%
            dplyr::mutate(filter_methyDiff = purrr::map(methy_comparison, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-methy_comparison) %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) -> gene_list_meth_diff

          meth_survival %>%
            dplyr::mutate(filter_SurDiff = purrr::map(diff_pval, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-diff_pval) -> gene_list_meth_sur

          meth_cor %>%
            dplyr::mutate(filter_cor = purrr::map(spm, filter_gene_list, gene_list = gene_set$match)) %>%
            dplyr::select(-spm) -> gene_list_meth_cor

          # get cancer type meth ----------------------------------------
          if (nrow(gene_list_meth_diff) > 0) {
            gene_list_meth_diff %>%
              tidyr::unnest() %>%
              tidyr::drop_na() -> gene_list_cancer_methdiff


            # ploting -----------------------------------------------------------------
            # meth diff point ----
            gene_list_cancer_methdiff %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(diff)) %>%
              dplyr::arrange(rank) -> gene_rank.methdiff
            gene_list_cancer_methdiff %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(diff)) %>%
              dplyr::arrange(rank) -> cancer_rank.methdiff

            callModule(methy_diff_pointPlot, "meth_diff", data = gene_list_cancer_methdiff, cancer = "cancer_types", gene = "symbol", size = "fdr", color = "diff", cancer_rank = cancer_rank.methdiff, gene_rank = gene_rank.methdiff, sizename = "-Log10(FDR)", colorname = "Methylation Diff (T - N)", title = "Methylation difference between tumor and normal samples.", status_monitor = "analysis", status, downloadname="Differential_methylation")
            .msg_meth_diff <- NULL
          } else {
            .msg_meth_diff <- paste(glue::glue("The [Differential Methylation] analysis based on paired sample in each cancer types.
In this analysis, only {nrow(meth_diff)} cancer types have paired samples. They are {paste0(meth_diff$cancer_types, collapse = ', ')}."), sep = " ")
            output[["meth_diff-plot"]] <- renderPlot({NULL})
          }

          # meth survival point ----
          gene_list_meth_sur %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) %>%
            tidyr::unnest() %>%
            tidyr::drop_na() -> gene_list_cancer_methsur
          if (nrow(gene_list_cancer_methsur) > 0) {
            gene_list_cancer_methsur %>%
              dplyr::mutate(a = ifelse(Hyper_worse == "Low", -1, 1)) %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(a)) %>%
              dplyr::arrange(rank) -> gene_rank.methsur

            gene_list_cancer_methsur %>%
              dplyr::mutate(a = ifelse(Hyper_worse == "Low", -1, 1)) %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(a)) %>%
              dplyr::arrange(rank) -> cancer_rank.methsur

            callModule(snv_sur_pointPlot, "meth_survival", data = gene_list_cancer_methsur, cancer = "cancer_types", gene = "symbol", size = "log10logrankP", color = "Hyper_worse", cancer_rank = cancer_rank.methsur, gene_rank = gene_rank.methsur, sizename = "logRank Pvalue", colorname = "Effect of HyperMethy on survival risk", title = "Overall survival difference between hypermethylation and hypomethylation.", status_monitor = "analysis", status, downloadname="Methylation_survival")
            .msg_meth_survival <- NULL
          } else {
            .msg_meth_survival <- paste(.msg, glue::glue("No significant [Methylation Survival] result of gene: {paste0(gene_set$match, collapse = ', ')} in your selected cancer type: {paste0(cancer_in_tcga_data_meth,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
            output[["meth_survival-plot"]] <- renderPlot({
              NULL
            })
          }

          # meth correlate to expression point ----
          gene_list_meth_cor %>%
            dplyr::filter(cancer_types %in% selected_ctyps()) %>%
            tidyr::unnest() %>%
            tidyr::drop_na() -> gene_list_cancer_methcor

          if (nrow(gene_list_cancer_methcor) > 0) {
            gene_list_cancer_methcor %>%
              dplyr::group_by(symbol) %>%
              dplyr::summarise(rank = sum(spm)) %>%
              dplyr::arrange(rank) -> gene_rank.methcor

            gene_list_cancer_methcor %>%
              dplyr::group_by(cancer_types) %>%
              dplyr::summarise(rank = sum(spm)) %>%
              dplyr::arrange(rank) -> cancer_rank.methcor

            callModule(methy_diff_pointPlot, "meth_exp", data = gene_list_cancer_methcor, cancer = "cancer_types", gene = "symbol", size = "logfdr", color = "spm", cancer_rank = cancer_rank.methcor, gene_rank = gene_rank.methcor, sizename = "-Log10(P.value)", colorname = "Spearman Correlation Coefficient", title = "Spearman Correlation Coefficient of methylation and gene expression.", status_monitor = "analysis", status, downloadname="Methylation_affect_exp")
            .msg_meth_exp <- NULL
          } else {
            .msg_meth_exp <- paste(.msg, glue::glue("No significant [Methylation to Expression] result of gene: {paste0(gene_set$match, collapse = ', ')} in your selected cancer type: {paste0(cancer_in_tcga_data_meth,collapse=', ')}. Please try more cancers or more genes."), sep = " ")
            output[["meth_exp-plot"]] <- renderPlot({
              NULL
            })
          }
          print(glue::glue("{paste0(rep('-', 10), collapse = '')} End methy part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

          # infomation UI for each part --------------------------
          output[["meth_diff-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_diff,style= "color:#CD3700")
            )
          })
          output[["meth_survival-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_survival,style= "color:#CD3700")
            )
          })
          output[["meth_exp-massage"]] <- renderUI({
            tagList(
              shiny::tags$p(.msg_meth_exp,style= "color:#CD3700")
            )
          })

          .msg <- paste(.msg, glue::glue("Since we just show significant results, so a small size of gene and cancer set may cause no significant result in some plots, if it happens, try more genes and cancer types."), sep = " ")
          # alert for information
          shinyBS::createAlert(
            session = session, anchorId = "meth-no_gene_set", title = "Information", style = "info",
            content = .msg, append = FALSE
          )
        } else {
          shinyBS::createAlert(
            session = session, anchorId = "meth-no_gene_set", title = "Oops",
            content = "No input gene set! Please go to Welcome page to input gene set.", style = "danger", append = FALSE
          )
        }
      }
    }
  }
)


# monitors -------------------------------------------------------
observe(meth_analysis())
# sourced by 'server.R'
# save as 'tcga_rppa_server.R'
# server elements 'tcga_rppa' sub tab of 'tcga' tab


source(file.path(config$wd, "functions", "tcga_rppa_function.R"))


# generate rppa result out ui -------------------------------------------------------

output$ui_rppa_result <- shiny::renderUI({
  fn_rppa_result(selected_analysis$rppa)
})

# analysis core -----------------------------------------------------------

rppa_analysis <- eventReactive(
  {
    status$analysis == TRUE
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      if (selected_analysis$rppa == TRUE) {
        # Cancer types value box selection ----------------------------------------
        
        callModule(module = cancerTypesSelect, id = "rppa", .sctps = intersect(selected_ctyps(), tcga_data))
        # Check box ---------------------------------------------------------------
        
        callModule(module = selectAndAnalysis, id = "rppa", .id = "rppa")
        load_data_rppa()
        
        # cancer overlap
        cancer_in_tcga_data_rppa <- intersect(selected_ctyps(),tcga_data)
        
      if (length(gene_set$match) != 0) {
        shinyBS::createAlert(
          session = session, anchorId = "rppa-no_gene_set", title = "Information", style = "info",
          content = "Need some time to draw picture, please wait.", append = FALSE
        )
        .msg <- c("NOTICE: Too much cancers and genes will make [Relation network] complicated, it will hard to see, try less genes or less cancers. [Global percentage] and [Heatmap percentage] will not change when cancer selection changes, cause it's a global percentage included all cancer types (see help page).")

        # remove pic result generate before ----
        # output$rppa_rela_plot <- renderImage({})
        
        # ploting -----------------------------------------------------------------
        # global plot -----------------------------------------------
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start rppa global analysis part@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
        
        # get gene set /cancer type data ----
        rppa_per %>%
          dplyr::filter(symbol %in% gene_set$match) -> gene_list_rppa_per
        if(nrow(gene_list_rppa_per)>0){
          gene_list_rppa_per %>%
            tidyr::unnest() %>%
            tidyr::gather(-symbol, -pathway, key = "class", value = "per") %>%
            dplyr::mutate(class = plyr::revalue(class, replace = c("a" = "Activation", "i" = "Inhibition", "n" = "None"))) -> rppa_pie_plot_ready
          
          # arugument for plot
          rppa_pie_height <- gene_set$match %>% length() * 0.25
          # if (rppa_pie_height > 15) {
          #   rppa_pie_height <- 15
          # }
          if (rppa_pie_height < 3) {
            rppa_pie_height <- 3
          }
          rppa_pie_outfile <- file.path(user_dir, "pngs", paste(user_id, "-", "TCGA_rppa_pie_profile.png", sep = ""))
          
          # draw ----
          callModule(rppaPiePlot, "rppa_pie", data = rppa_pie_plot_ready, y = "per", fill = "class", facet_grid = " symbol~pathway", height = rppa_pie_height, outfile = rppa_pie_outfile, status, downloadname = "Pathway_activity_pie_percentage")
          
          # rppa global percentage ----
          # data process
          gene_list_rppa_per %>%
            tidyr::unnest() %>%
            dplyr::filter(a + i > 5 / 32) %>%
            dplyr::select(-n) %>%
            tidyr::gather(-symbol, -pathway, key = "class", value = "per") %>%
            dplyr::mutate(per = ifelse(class == "i", -per * 100, per * 100)) %>%
            dplyr::mutate(class = plyr::revalue(class, replace = c("a" = "A", "i" = "I"))) %>%
            tidyr::unite(pathway, c(pathway, class)) -> rppa_per_ready
          
          # pic draw
          rppa_heat_height <- gene_set$match %>% length() * 0.1
          if (rppa_heat_height > 15) {
            rppa_heat_height <- 15
          }
          if (rppa_heat_height < 3) {
            rppa_heat_height <- 3
          }
          rppa_heat_outfile <- file.path(user_dir, "pngs", paste(user_id, "-", "TCGA_rppa_heatmap_percentage.png", sep = ""))
          callModule(rppa_heat_per, "rppa_per", rppa_per_ready = rppa_per_ready, pathway = "pathway", symbol = "symbol", per = "per", height = rppa_heat_height, outfile = rppa_heat_outfile, status, downloadname = "Pathway_activity_global_percentage")
          .msg_rppa_global <- NULL
        } else {
          .msg_rppa_global <- c("No significant result for your input genes in this part of analysis.")
          # callModule(white_plot, "rppa_per", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_1.png", sep = "")))
          output[["rppa_per-plot"]] <- renderPlot({NULL})
          callModule(white_plot, "rppa_pie", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_2.png", sep = "")))
        }
        # rppa global pie plot----
        
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start rppa global analysis part@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
        
        # rppa relation plot ---------------------------------------------------
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} start rppa analysis part@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
        
        # data processing
                rppa_relation %>%
          dplyr::filter(cancer_types %in% selected_ctyps()) %>%
          dplyr::mutate(data = purrr::map(data, filter_gene_list, gene_list = gene_set$match)) -> gene_list_cancer_rppa_rela
                
                if(nrow(gene_list_cancer_rppa_rela)>0){
                  gene_list_cancer_rppa_rela%>%
                    tidyr::unnest() -> gene_list_cancer_rppa_rela
                  if(nrow(gene_list_cancer_rppa_rela)>0) {
                    # rppa line contact ----
                    # get data
                    cancer_text <- get_rppa_text(gene_list_cancer_rppa_rela)
                    gene_list_cancer_rppa_rela %>%
                      dplyr::mutate(n=1:nrow(gene_list_cancer_rppa_rela)) %>%
                      tidyr::nest(-n) %>%
                      dplyr::group_by(n) %>%
                      dplyr::mutate(seg=purrr::map(data,.f=get_rppa_seg,cancer_text=cancer_text)) %>%
                      dplyr::ungroup() %>%
                      dplyr::select(-n,-data) %>%
                      tidyr::unnest() ->plot_seg
                    # get_rppa_seg1(gene_list_cancer_rppa_rela,cancer_text = cancer_text) -> plot_seg
                    cancer_text %>%
                      dplyr::filter(type == "cancer") -> cancer.text
                    cancer_text %>%
                      dplyr::filter(type == "gene") -> gene.text
                    cancer_text %>%
                      dplyr::filter(type == "pathway") -> path.text
                    rppa_line_height <- gene_set$match %>% length() * 0.15
                    if (rppa_heat_height > 15) {
                      rppa_heat_height <- 15
                    }
                    if (rppa_heat_height < 3) {
                      rppa_heat_height <- 3
                    }
                    # plot draw
                    output$`rppa_line-plot` <- renderImage({
                      status[["analysis"]]

                      rppa_line_outfile <- file.path(user_dir, "pngs", paste(user_id, "-", "TCGA_rppa_network_profile.png", sep = ""))

                      ggsave(rppa_line_outfile, rppa_line_contact(plot_seg, cancer.text, gene.text, path.text), device = "png", width = 3, height = rppa_line_height)
                      list(
                        src = rppa_line_outfile,
                        contentType = "image/png",
                        # width = "100%" ,
                        # height = 900,
                        alt = "This is alternate text"
                      )
                    }, deleteFile = FALSE)
                    
                    # output$`rppa_line-plot` <- renderPlot({
                    #   status[["analysis"]]
                    #   rppa_line_contact(plot_seg, cancer.text, gene.text, path.text)
                    # },height = function(){ifelse(rppa_line_height<200,200,rppa_line_height)})
                    
                    output$`rppa_line-picdownload` <- downloadHandler(
                      filename = function() { paste("Pathway_relation_network", '.',input$`rppa_line-pictype`, sep='') },
                      content = function(file) {
                        ggsave(file,rppa_line_contact(plot_seg, cancer.text, gene.text, path.text),device = input$`rppa_line-pictype`,width = input$`rppa_line-d_width`,height = input$`rppa_line-d_height`)
                      }
                    )
                    .msg_rppa_line <- NULL
                  } else{
                    .msg_rppa_line <- paste(glue::glue("No regulation relationship between {paste0(gene_set$match, collapse = ', ')} and pathway activity in {paste0(cancer_in_tcga_data_rppa,collapse=', ')}. Please try more cancers or more genes."))
                    # output[["rppa_line-plot"]] <- renderPlot({NULL})
                    callModule(white_plot, "rppa_line", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_3.png", sep = "")))
                  }
                } else {
                  .msg_rppa_line <- c("No significant result in this [Regulation network], try more genes or cancers.")
                  # output[["rppa_line-plot"]] <- renderPlot({NULL})
                  callModule(white_plot, "rppa_line", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_3.png", sep = "")))
                }
                
                # create error message -------
                output[["rppa_line-massage"]] <- renderUI({
                  tagList(
                    shiny::tags$p(.msg_rppa_line, style = "color:#CD3700")
                  )
                })
                output[["rppa_pie-massage"]] <- renderUI({
                  tagList(
                    shiny::tags$p(.msg_rppa_global, style = "color:#CD3700")
                  )
                })
                output[["rppa_per-massage"]] <- renderUI({
                  tagList(
                    shiny::tags$p(.msg_rppa_global, style = "color:#CD3700")
                  )
                })
        
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} End rppa relation analysis part@ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
        # alert for information
        shinyBS::createAlert(
          session = session, anchorId = "rppa-no_gene_set", title = "Information", style = "info",
          content = .msg, append = FALSE
        )
      } else {
        shinyBS::createAlert(
          session = session, anchorId = "rppa-no_gene_set", title = "Oops",
          content = "No input gene set! Please go to Welcome page to input gene set.", style = "danger", append = FALSE
        )
      }
      }
    }
  }
)



# monitor -----------------------------------------------------------------

observe(rppa_analysis())
# sourced by 'server.R'
# save as 'tcga_snv_server.R'
# server elements 'tcga_snv' sub tab of 'tcga' tab

source(file.path(config$wd, "functions", "tcga_snv_function.R"))


#  get cancer type --------------------------------------------------------
# selected_ctyps <- callModule(cancerType, "snv")



# button control --------------------------------------------------------

# submit cancer type -------------------------------------------------------

# snv_submit_analysis <- function(input, output, session){
#   observeEvent(input$submit, {
#     status$snv_submit <- TRUE
#     print(status$snv_submit)
#   })
# }
#
# callModule(snv_submit_analysis,"snv")

# analysis core -----------------------------------------------------------
# monitor for gene list change-----------------------------------
# snv_gene_list <- eventReactive(
#   eventExpr = status$analysis,
#   ignoreNULL = TRUE,
#   valueExpr = {
#     # be sure the following code run after start analysis
#     if (status$analysis == TRUE) {
#       status$snv_submit <- TRUE
#       shinyjs::disable(id = "snv-submit")
#       shinyjs::disable(id = "snv-switch")
#       as.character(gene_set$match)
#     }
#   }
# )

# snv result out ui -------------------------------------------------------

output$ui_snv_result <- shiny::renderUI({
  fn_snv_result(selected_analysis$snv)
})

# get gene set snv --------------------------------------------------------
snv_analysis <- eventReactive(
  {
    eventExpr <- status$analysis
  },
  ignoreNULL = TRUE,
  valueExpr = {
    if (status$analysis == TRUE) {
      # laod snv data ----
      if (selected_analysis$snv == TRUE) {
        load_data_snv()
        # Cancer types value box selection ----------------------------------------

        callModule(module = cancerTypesSelect, id = "snv", .sctps = intersect(selected_ctyps(), tcga_data))
        # Check box ---------------------------------------------------------------

        callModule(module = selectAndAnalysis, id = "snv", .id = "snv")
        if (length(gene_set$match) != 0) {
          if (length(selected_ctyps() != 0)) {
            .msg <- c("NOTICE: ")

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} Start snv part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # cancer overlap
            .cancer_in_tcga_data_snv <- intersect(selected_ctyps(),tcga_data)
            print(.cancer_in_tcga_data_snv)

            # snv percent -------------------------------------------------------------
            snv %>%
              dplyr::mutate(filter_snv = purrr::map(mut_count, filter_gene_list, gene_list = gene_set$match)) %>%
              dplyr::select(-mut_count) %>%
              dplyr::filter(cancer_types %in% selected_ctyps()) -> gene_list_cancer_snv

            # plot out ----------------------------------------------------------------

            # snv percentage ----------------------------------------------------------
            if (nrow(gene_list_cancer_snv) > 0) {
              gene_list_cancer_snv %>%
                tidyr::unnest(filter_snv) ->gene_list_cancer_snv
              if(nrow(gene_list_cancer_snv)>0){
                gene_list_cancer_snv %>%
                  tidyr::drop_na() %>%
                  dplyr::mutate(x_label = paste(cancer_types, " (n=", n, ")", sep = "")) %>%
                  #dplyr::mutate(sm_count = ifelse(sm_count > 0, sm_count, NA)) %>%
                  dplyr::mutate(per = ifelse(per > 0.02, per, 0)) -> snv_per_plot_ready
                snv_per_plot_ready %>%
                  dplyr::group_by(x_label) %>%
                  dplyr::summarise(s = sum(per)) %>%
                  dplyr::arrange(dplyr::desc(s)) -> snv_per_cancer_rank
                snv_per_plot_ready %>%
                  dplyr::group_by(symbol) %>%
                  dplyr::summarise(s = sum(sm_count)) %>%
                  dplyr::arrange(s) -> snv_per_gene_rank

                callModule(
                  snv_per_heatmap, "snv_percentage", data = snv_per_plot_ready,
                  cancer = "x_label", gene = "symbol", fill = "per", label = "sm_count",
                  cancer_rank = snv_per_cancer_rank, gene_rank = snv_per_gene_rank, status_monitor = "analysis", status,
                  downloadname = "SNV_percentage_profile"
                )
                .msg_snv_percentage <- NULL
              } else {
                .msg_snv_percentage <- paste(glue::glue("No significant [SNV percentage profile] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(.cancer_in_tcga_data_snv, collapse=', ')}. Please try other cancers or genes."), sep = " ")
                output[["snv_percentage-plot"]] <- renderPlot({
                  NULL
                })
              }
            } else {
              .msg_snv_percentage <- paste(glue::glue("No significant [SNV percentage profile] result of gene: {paste0(gene_set$match, collapse = ',')} in your selected cancer types: {paste0(.cancer_in_tcga_data_snv, collapse=', ')}. Please try other cancers or genes."), sep = " ")
              output[["snv_percentage-plot"]] <- renderPlot({
                NULL
              })
            }
            # snv survival ------------------------------------------------------------
            snv_survival %>%
              dplyr::mutate(filter_survival = purrr::map(diff_pval, filter_gene_list, gene_list = gene_set$match)) %>%
              dplyr::select(-diff_pval) %>%
              dplyr::filter(cancer_types %in% selected_ctyps()) %>%
              tidyr::unnest() %>%
              dplyr::mutate(logP = -log10(logRankP)) %>%
              dplyr::mutate(logP = ifelse(logP > 15, 15, logP)) %>%
              dplyr::mutate(logP = ifelse(logP < -log10(0.05), NA, logP)) %>%
              tidyr::drop_na() -> snv_sur_plot_ready -> snv_sur_plot_ready
            # survival ----------------------------------------------------------------
            if (nrow(snv_sur_plot_ready) > 0) {
              snv_sur_plot_ready %>%
                dplyr::mutate(s = ifelse(estimate > 0, 1, -1)) %>%
                dplyr::mutate(s = ifelse(logRankP > 0.05, 0, s)) %>%
                dplyr::group_by(cancer_types) %>%
                dplyr::summarise(r = sum(s)) %>%
                dplyr::arrange(dplyr::desc(r)) -> snv_sur_cancer_rank

              snv_sur_plot_ready %>%
                dplyr::mutate(s = ifelse(estimate > 0, 1, -1)) %>%
                dplyr::mutate(s = ifelse(logRankP > 0.05, 0, s)) %>%
                dplyr::group_by(symbol) %>%
                dplyr::summarise(r = sum(s)) %>%
                dplyr::arrange(dplyr::desc(r)) -> snv_sur_gene_rank
              callModule(
                snv_sur_pointPlot, "snv_survival", data = snv_sur_plot_ready, cancer = "cancer_types",
                gene = "symbol", size = "logP", color = "worse", cancer_rank = snv_sur_cancer_rank,
                gene_rank = snv_sur_gene_rank, sizename = "logRank P", colorname = "Effect of mutation on survival risk", title = "OS between mut and non-mut genes.", status_monitor = "analysis", status, downloadname = "SNV_affect_survival"
              )
              .msg_snv_survival <- NULL
            } else {
              .msg_snv_survival <- paste(glue::glue("No significant [SNV survival] result of gene: {paste0(gene_set$match, collapse = ', ')} in your selected cancer type: {paste0(.cancer_in_tcga_data_snv, collapse=', ')}. Please try other cancers or genes."), sep = " ")
              output[["snv_survival-plot"]] <- renderPlot({
                NULL
              })
            }


            .msg <- paste(.msg, glue::glue("Please be patient, need some time to draw pictrue. Since we just show significant results, so a small size of gene and cancer set may cause no significant result in some plots. If it happens, try more genes and cancer types."), sep = " ")
            # alert for information
            shinyBS::createAlert(
              session = session, anchorId = "snv-no_gene_set", title = "Information", style = "info",
              content = .msg, append = FALSE
            )

            # maf ---------------------------------------------------------------------
            if (length(gene_set$match) >= 2) {
              snv_InpSel <- paste0(selected_ctyps(), collapse = "','")
              query <- as.expression(paste0("Cancer_Types %in% c('", snv_InpSel, "')"))
              # my_subsetMaf(mc3_pass, genes = gene_set$match, mafObj = T,query = query) -> gene_list_maf #
              gene_set_in_maf <- intersect(maf_gene_all,gene_set$match) %>% length()
              cancer_in_maf <- intersect(maf_cancer_all,selected_ctyps())
              cancer_noin_maf <- setdiff(.cancer_in_tcga_data_snv,maf_cancer_all)

              if(length(cancer_noin_maf)>0){
                cancer_no_data.msg <- glue::glue(" {paste0(cancer_noin_maf, collapse=',')} don't have data in this analysis.")
              } else {cancer_no_data.msg <- NULL}

              if(gene_set_in_maf>=2 && length(cancer_in_maf)>0){
                maftools::subsetMaf(mc3_pass, genes = gene_set$match, mafObj = T, query = query) -> gene_list_maf
                # 1. snv summary
                snv_su_out <- file.path(user_dir, "pngs", paste(user_id, "-SNV_summary_profile.png", sep = ""))
                callModule(snv_maf_summaryPlot, "snv_summary", gene_list_maf = gene_list_maf, outfile = snv_su_out, status_monitor = "analysis", status, downloadname = "SNV_summary")

                # 2. oncoplot
                snv_onco_out <- file.path(user_dir, "pngs", paste(user_id, "-SNV_oncoplot_profile.png", sep = ""))
                callModule(snv_maf_oncoPlot, "snv_oncoplot", gene_list_maf = gene_list_maf, pancan_color = pancan_color, outfile = snv_onco_out, status_monitor = "analysis", status, downloadname = "SNV_oncoplot")
                .msg_snv_oncoplot <- cancer_no_data.msg
                .msg_snv_summary <- cancer_no_data.msg
              } else {
                .msg_snv_oncoplot <- paste(glue::glue("Your selected genes: {paste0(gene_set$match, collapse = ', ')} are not mutate in your selected cancer type: {paste0(.cancer_in_tcga_data_snv, collapse = ', ')}.{cancer_no_data.msg} Please try other cancers or genes."), sep = " ")
                .msg_snv_summary <- paste(glue::glue("Your selected genes: {paste0(gene_set$match, collapse = ', ')} are not mutate in your selected cancer type: {paste0(.cancer_in_tcga_data_snv, collapse = ', ')}.{cancer_no_data.msg} Please try other cancers or genes."), sep = " ")

                callModule(white_plot, "snv_oncoplot", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_2.png", sep = "")))
                callModule(white_plot, "snv_summary", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_1.png", sep = "")))
              }
            } else {
              .msg_snv_oncoplot <- "Cannot create SNV oncoplot for single gene. Minimum two genes required ! "
              .msg_snv_summary <- "Cannot create SNV summary plot for single gene. Minimum two genes required ! "
              callModule(white_plot, "snv_oncoplot", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_2.png", sep = "")))
              callModule(white_plot, "snv_summary", status_monitor = "analysis", status = status, outfile = file.path(user_dir, "pngs", paste(user_id, "-white_1.png", sep = "")))
            }

            print(glue::glue("{paste0(rep('-', 10), collapse = '')} End maf part analysis @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))

            # infomation UI for each part
            output[["snv_percentage-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_snv_percentage, style = "color:#CD3700")
              )
            })
            output[["snv_summary-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_snv_summary, style = "color:#CD3700")
              )
            })
            output[["snv_oncoplot-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_snv_oncoplot, style = "color:#CD3700")
              )
            })
            output[["snv_survival-massage"]] <- renderUI({
              tagList(
                shiny::tags$p(.msg_snv_survival, style = "color:#CD3700")
              )
            })
          } else {
            shinyBS::createAlert(
              session = session, anchorId = "snv-no_cancer_set", title = "Oops",
              content = "No cancer selected! Please select at least one cancer type.", style = "danger", append = FALSE
            )
          }
        } else {
          shinyBS::createAlert(
            session = session, anchorId = "snv-no_gene_set", title = "Oops",
            content = "No input gene set! Please go to Welcome page to input gene set.", style = "danger", append = FALSE
          )
        }
      }
    }
  }
)

# monitors -------------------------------------------------------
# observe(snv_global_analysis())
observe(snv_analysis())
# source by server.R
# source by ui.R
# saved as functions_server.R

##########################################
# Cancer selection data get and comfim####
## 1. cancer select for each part ui ######
## 2. cancer type selection confirm  ######
## cancerTypeInput & cancerType############
##########################################

# GTEx normal tissue choice ##############

GTEx_expr_Brain_choice <- list("Brain" = "Brain")
GTEx_expr_Liver_choice <- list("Liver" = "Liver")
GTEx_expr_Heart_choice <- list("Heart" = "Heart")
GTEx_expr_Ovary_choice <- list("Ovary" = "Ovary")
GTEx_expr_Lung_choice <- list("Lung" = "Lung")
GTEx_expr_Breast_choice <- list("Breast" = "Breast")
GTEx_expr_Skin_choice <- list("Skin" = "Skin")
GTEx_expr_Blood_choice <- list("Blood" = "Blood")
GTEx_expr_Testis_choice <- list("Testis" = "Testis")
GTEx_expr_Colon_choice <- list("Colon" = "Colon")
GTEx_expr_other_choice <- list(
  "Adipose Tissue" = "Adipose Tissue", "Muscle" = "Muscle",
  "Blood Vessel" = "Blood Vessel", "Salivary Gland" = "Salivary Gland", "Adrenal Gland" = "Adrenal Gland",
  "Thyroid" = "Thyroid", "Spleen" = "Spleen", "Small Intestine" = "Small Intestine",
  "Cervix Uteri" = "Cervix Uteri", "Bladder" = "Bladder", "Fallopian Tube" = "Fallopian Tube",
  "Uterus" = "Uterus", "Pituitary" = "Pituitary", "Esophagus" = "Esophagus",
  "Nerve" = "Nerve", "Vagina" = "Vagina", "Pancreas" = "Pancreas", "Prostate" = "Prostate",
  "Stomach" = "Stomach", "Kidney" = "Kidney"
)
GTEx_expr_input_selection <- paste("input$", grep("GTEx_expr_.*_choice", ls(), value = T), sep = "")


##### GTEx tissue 4 UI#####
GTExTissueType <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # GTEx tissues selection----
      column(
        width = 10,
        offset = 1,
        shiny::tags$br(),
        shiny::tags$h3("GTEx tissues selection", class = "text-success"),
        shiny::tags$br(),

        shinydashboard::tabBox(
          width = 12, title = "Tissue",
          tabPanel(
            "Brain",
            shiny::tags$h4("Brain", class = "text-success"),
            checkboxGroupButtons(
              inputId = ns("Brain"), label = NULL,
              choices = GTEx_expr_Brain_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Liver",
            checkboxGroupButtons(
              inputId = ns("Liver"), label = NULL,
              choices = GTEx_expr_Liver_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Heart",
            checkboxGroupButtons(
              inputId = ns("Heart"), label = NULL,
              choices = GTEx_expr_Heart_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Ovary",
            checkboxGroupButtons(
              inputId = ns("Ovary"), label = NULL,
              choices = GTEx_expr_Ovary_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Lung",
            checkboxGroupButtons(
              inputId = ns("Lung"), label = NULL,
              choices = GTEx_expr_Lung_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Breast",
            checkboxGroupButtons(
              inputId = ns("Breast"), label = NULL,
              choices = GTEx_expr_Breast_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Skin",
            checkboxGroupButtons(
              inputId = ns("Skin"), label = NULL,
              choices = GTEx_expr_Skin_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Blood",
            checkboxGroupButtons(
              inputId = ns("Blood"), label = NULL,
              choices = GTEx_expr_Blood_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Testis",
            checkboxGroupButtons(
              inputId = ns("Testis"), label = NULL,
              choices = GTEx_expr_Testis_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Colon",
            checkboxGroupButtons(
              inputId = ns("Colon"), label = NULL,
              choices = GTEx_expr_Colon_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Other tissues",
            checkboxGroupButtons(
              inputId = ns("gtex_expr_other_tissue"), label = NULL,
              choices = GTEx_expr_other_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          )
        )
      ),
      shiny::tags$hr(width = "85%")
    )
  )
}



GTEx_normal_Tissue <- function(input, output, session) {
  GTEx_normal_tissue <- reactive({
    eval(parse(text = GTEx_expr_input_selection)) -> GTEx_normal_tissue
  })

  return(GTEx_normal_tissue)
}

GTEx_eqtl_Tissue <- function(input, output, session) {
  GTEx_eqtl_tissue <- reactive({
    eval(parse(text = GTEx_eqtl_input_selection)) -> GTEx_eqtl_tissue
  })

  return(GTEx_eqtl_tissue)
}

# sub normal tissue----------

sub_GTEx_tissue <- list(
  "Brain" = "Brain",
  "Liver" = "Liver",
  "Heart" = "Heart",
  "Ovary" = "Ovary",
  "Lung" = "Lung",
  "Breast" = "Breast",
  "Skin" = "Skin",
  "Blood" = "Blood",
  "Testis" = "Testis",
  "Colon" = "Colon",
  "gtex_expr_other_tissue" = c(
    "Adipose Tissue", "Muscle", "Blood Vessel",
    "Salivary Gland", "Adrenal Gland", "Thyroid",
    "Spleen", "Small Intestine", "Cervix Uteri",
    "Bladder", "Fallopian Tube", "Uterus",
    "Pituitary", "Esophagus", "Nerve",
    "Vagina", "Pancreas", "Prostate", "Stomach", "Kidney"
  )
)

resetGTExTissueType <- function(input, output, session) {
  for (i in c(tabPannel_element_ten, "other_tissue")) {
    shinyjs::reset(i)
  }
  GTEx_normal_tissue <- reactive({
    c("") -> GTEx_normal_tissue
  })
  return(GTEx_normal_tissue)
}



##### GTEx eqtl tissue ####
GTEx_eqtl_Brain_choice <- list("Brain" = "Brain")
GTEx_eqtl_Liver_choice <- list("Liver" = "Liver")
GTEx_eqtl_Heart_choice <- list("Heart" = "Heart")
GTEx_eqtl_Ovary_choice <- list("Ovary" = "Ovary")
GTEx_eqtl_Lung_choice <- list("Lung" = "Lung")
GTEx_eqtl_Breast_choice <- list("Breast" = "Breast")
GTEx_eqtl_Skin_choice <- list("Skin" = "Skin")
GTEx_eqtl_Blood_choice <- list("Blood" = "Blood")
GTEx_eqtl_Testis_choice <- list("Testis" = "Testis")
GTEx_eqtl_Colon_choice <- list("Colon" = "Colon")
GTEx_eqtl_other_choice <- list(
  "Adipose Tissue" = "Adipose Tissue", "Muscle" = "Muscle",
  "Artery" = "Artery", "Salivary Gland" = "Salivary Gland", "Adrenal Gland" = "Adrenal Gland",
  "Thyroid" = "Thyroid", "Spleen" = "Spleen", "Small Intestine" = "Small Intestine",
  "Uterus" = "Uterus", "Pituitary" = "Pituitary", "Esophagus" = "Esophagus",
  "Nerve" = "Nerve", "Vagina" = "Vagina", "Pancreas" = "Pancreas", "Prostate" = "Prostate",
  "Stomach" = "Stomach", "Artery" = "Artery", "Cells" = "Cells"
)
GTEx_eqtl_input_selection <- paste("input$", grep("GTEx_eqtl_.*_choice", ls(), value = T), sep = "")

GTExTissueeqtl <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # GTEx tissues selection----
      column(
        width = 10,
        offset = 1,
        shiny::tags$br(),
        shiny::tags$h3("GTEx tissues selection", class = "text-success"),
        shiny::tags$br(),

        shinydashboard::tabBox(
          width = 12, title = "Tissue",
          tabPanel(
            "Brain",
            shiny::tags$h4("Brain", class = "text-success"),
            checkboxGroupButtons(
              inputId = ns("Brain"), label = NULL,
              choices = GTEx_eqtl_Brain_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Liver",
            checkboxGroupButtons(
              inputId = ns("Liver"), label = NULL,
              choices = GTEx_eqtl_Liver_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Heart",
            checkboxGroupButtons(
              inputId = ns("Heart"), label = NULL,
              choices = GTEx_eqtl_Heart_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Ovary",
            checkboxGroupButtons(
              inputId = ns("Ovary"), label = NULL,
              choices = GTEx_eqtl_Ovary_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Lung",
            checkboxGroupButtons(
              inputId = ns("Lung"), label = NULL,
              choices = GTEx_eqtl_Lung_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Breast",
            checkboxGroupButtons(
              inputId = ns("Breast"), label = NULL,
              choices = GTEx_eqtl_Breast_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Skin",
            checkboxGroupButtons(
              inputId = ns("Skin"), label = NULL,
              choices = GTEx_eqtl_Skin_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Blood",
            checkboxGroupButtons(
              inputId = ns("Blood"), label = NULL,
              choices = GTEx_eqtl_Blood_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Testis",
            checkboxGroupButtons(
              inputId = ns("Testis"), label = NULL,
              choices = GTEx_eqtl_Testis_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Colon",
            checkboxGroupButtons(
              inputId = ns("Colon"), label = NULL,
              choices = GTEx_eqtl_Colon_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          ),
          tabPanel(
            "Other tissues",
            checkboxGroupButtons(
              inputId = ns("gtex_eqtl_other_tissue"), label = NULL,
              choices = GTEx_eqtl_other_choice,
              justified = TRUE,
              checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon")),
              direction = "vertical",
              individual = TRUE
            )
          )
        )
      ),
      shiny::tags$hr(width = "85%")
    )
  )
}

sub_GTEx_eqtl_tissue <- list(
  "Brain" = "Brain",
  "Liver" = "Liver",
  "Heart" = "Heart",
  "Ovary" = "Ovary",
  "Lung" = "Lung",
  "Breast" = "Breast",
  "Skin" = "Skin",
  "Blood" = "Blood",
  "Testis" = "Testis",
  "Colon" = "Colon",
  "gtex_eqtl_other_tissue" = c(
    "Adipose Tissue", "Muscle", "Artery",
    "Salivary Gland", "Adrenal Gland",
    "Thyroid", "Spleen", "Small Intestine",
    "Uterus", "Pituitary", "Esophagus",
    "Nerve", "Vagina", "Pancreas", "Prostate",
    "Stomach", "Artery", "Cells"
  )
)
# cancer type choice ------------------------------------------------------

Kidney_choice <- list(
  "Kidney Chromophobe(KICH)" = "KICH",
  "Kidney Renal Clear Cell Carcinoma(KIRC)" = "KIRC",
  "Kidney Renal Papillary Cell Carcinoma(KIRP)" = "KIRP"
)
Adrenal_Gland_choice <- list(
  "Adrenocortical Carcinoma(ACC)" = "ACC",
  "Pheochromocytoma and Paraganglioma(PCPG)" = "PCPG"
)
Brain_choice <- list(
  "Glioblastoma Multiforme(GBM)" = "GBM",
  "Brain Lower Grade Glioma(LGG)" = "LGG"
)
Colorectal_choice <- list(
  "Colon Adenocarcinoma(COAD)" = "COAD",
  "Rectum Adenocarcinoma(READ)" = "READ"
)
Lung_choice <- list(
  "Lung Adenocarcinoma(LUAD)" = "LUAD",
  "Lung Squamous Cell Carcinoma(LUSC)" = "LUSC"
)
Uterus_choice <- list(
  "Uterine Corpus Endometrial Carcinoma(UCEC)" = "UCEC",
  "Uterine Carcinosarcoma(UCS)" = "UCS"
)
Bile_Duct_choice <- list("Bladder Urothelial Carcinoma(BLCA)" = "BLCA")
Bone_Marrow_choice <- list("Acute Myeloid Leukemia(LAML)" = "LAML")
Breast_choice <- list("Breast Invasive Carcinoma(BRCA)" = "BRCA")
Cervix_choice <- list("Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma(CESC)" = "CESC")
other_tissue_choice <- list(
  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma(DLBC)" = "DLBC",
  "Esophageal Carcinoma(ESCA)" = "ESCA",
  "Stomach Adenocarcinoma(STAD)" = "STAD",
  "Head and Neck Squamous Cell Carcinoma(HNSC)" = "HNSC",
  "Liver Hepatocellular Carcinoma(LIHC)" = "LIHC",
  "Mesothelioma(MESO)" = "MESO",
  "Ovarian Serous Cystadenocarcinoma(OV)" = "OV",
  "Pancreatic Adenocarcinoma(PAAD)" = "PAAD",
  "Prostate Adenocarcinoma(PRAD)" = "PRAD",
  "Sarcoma(SARC)" = "SARC",
  "Skin Cutaneous Melanoma(SKCM)" = "SKCM",
  "Testicular Germ Cell Tumors(TGCT)" = "TGCT",
  "Thyroid Carcinoma(THCA)" = "THCA",
  "Thymoma(THYM)" = "THYM",
  "Uveal Melanoma(UVM)" = "UVM",
  "Cholangiocarcinoma(CHOL)" = "CHOL"
)

# cancer type selection ---------------------------------------------------
cancerTypeInput <- function(id) {
  ns <- NS(id)

  tagList(
    # value box for selected cancer types ----
    fluidRow(shiny::uiOutput(outputId = ns("cancer_types_select")))
    # shiny::tags$hr(width = "85%")
  )
}

tissueTypeInput <- function(id) {
  ns <- NS(id)

  tagList(
    # value box for selected cancer types ----
    fluidRow(shiny::uiOutput(outputId = ns("tissue_types_select")))
    # shiny::tags$hr(width = "85%")
  )
}

# Value box for selection cancer types ------------------------------------


cancerTypesSelect <- function(input, output, session, .sctps) {
  output$cancer_types_select <- renderUI({
    shiny::tagList(
      column(
        width = 4, offset = 2,
        tags$head(tags$style(HTML('.info-box {min-height: 60px;} .info-box-icon {height: 60px; line-height: 60px;} .info-box-content {padding-top: 2px; padding-bottom: 2px;}'))),
        infoBox(
          title = "Number of selected cancers", value = length(.sctps),
          width = 12, color = "aqua", fill = TRUE
        ) # icon = icon("users"),
      ),
      column(
        width = 4,
        infoBox(
          title = "Number of unselected cancers", value = 33 - length(.sctps),
          width = 12, color = "red", fill = TRUE
        ) # icon = icon("credit-card"),
      ),
      column(
        width = 8, offset = 2,
        box(
          solidHeader = TRUE, status = "primary",
          title = "Selected Cancer Types", width = 12,
          paste0(.sctps, collapse = ", ")
        )
      )
    )
  })
}

tissueTypesSelect <- function(input, output, session, .sctps) {
  output$tissue_types_select <- renderUI({
    print(.sctps)
    shiny::tagList(
      column(
        width = 4, offset = 2,
        infoBox(
          title = "Number of selected tissues", value = length(.sctps),
          width = 12, color = "aqua", fill = TRUE
        ) # icon = icon("users"),
      ),
      column(
        width = 4,
        infoBox(
          title = "Number of unselected tissues", value = 30 - length(.sctps),
          width = 12, color = "red", fill = TRUE
        ) # icon = icon("credit-card"),
      ),
      column(
        width = 8, offset = 2,
        box(
          solidHeader = TRUE, status = "primary",
          title = "Selected Tissue Types", width = 12,
          paste0(.sctps, collapse = ", ")
        )
      )
    )
  })
}
# select and submit for UI----

selectAndAnalysisInput <- function(id) {
  ns <- NS(id)
  shiny::tagList(
    fluidRow(
      column(
        width = 8, offset = 2,
        shinyBS::bsAlert(anchorId = ns("no_gene_set")),
        shinyBS::bsAlert(anchorId = ns("no_paired_sample")),
        shinyBS::bsAlert(anchorId = ns("note"))
      )
    )
  )
}


# Simplified cancer types -------------------------------------------------

sub_cancer_types <- list(
  Kidney = c("KICH", "KIRC", "KIRP"),
  Adrenal_Gland = c("ACC", "PCPG"),
  Brain = c("GBM", "LGG"),
  Colorectal = c("COAD", "READ"),
  Lung = c("LUAD", "LUSC"),
  Uterus = c("UCEC", "UCS"),
  Bile_Duct = c("BLCA"),
  Bone_Marrow = c("LAML"),
  Breast = c("BRCA"),
  Cervix = c("CESC"),
  other_tissue = c("DLBC", "ESCA", "STAD", "HNSC", "LIHC", "MESO", "OV", "PAAD", "PRAD", "SARC", "SKCM", "TGCT", "THCA", "THYM", "UVM", "CHOL")
)


# Check and uncheck submit ------------------------------------------------

check_sub_cancer_types <- function(input, output, session, .cts, .check) {
  names(.cts) %>%
    purrr::walk(
      .f = function(.x) {
        .selected <- if (.check) .cts[[.x]] else character(0)
        updateCheckboxGroupButtons(
          session = session,
          inputId = .x, selected = .selected,
          checkIcon = list(yes = icon("ok", lib = "glyphicon"), no = icon("remove", lib = "glyphicon"))
        )
      }
    )
}


# Call by server to update check ------------------------------------------

selectAndAnalysis <- function(input, output, session, .id) {
  observeEvent(
    eventExpr = input$switch,
    handlerExpr = {
      if (input$switch) {
        check_sub_cancer_types(input, output, session, sub_cancer_types, TRUE)
        check_sub_cancer_types(input, output, session, sub_GTEx_tissue, TRUE)
        check_sub_cancer_types(input, output, session, sub_GTEx_eqtl_tissue, TRUE)
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} Select all {.id} @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
      } else {
        check_sub_cancer_types(input, output, session, sub_cancer_types, FALSE)
        check_sub_cancer_types(input, output, session, sub_GTEx_tissue, FALSE)
        check_sub_cancer_types(input, output, session, sub_GTEx_eqtl_tissue, FALSE)
        print(glue::glue("{paste0(rep('-', 10), collapse = '')} Deselect all {.id} @ {Sys.time()} {paste0(rep('-', 10), collapse = '')}"))
      }
    }
  )
}


# cancerType server function ----------------------------------------------
# pair with cancerTypeInput in functions_ui.R##
# Call by *_*_server.R by callModule(cancerType,"id pair with UI part")
cancerType <- function(input, output, session) {
  reactive({
    c(
      input$Kidney, input$Adrenal_Gland, input$Brain, input$Colorectal,
      input$Lung, input$Uterus, input$Bile_Duct, input$Bone_Marrow, input$Breast,
      input$Cervix, input$other_tissue
    )
  })
}

resetcancerType <- function(input, output, session) {
  shinyjs::reset("Kidney")
  shinyjs::reset("Adrenal_Gland")
  shinyjs::reset("Brain")
  shinyjs::reset("Colorectal")
  shinyjs::reset("Lung")
  shinyjs::reset("Uterus")
  shinyjs::reset("Bile_Duct")
  shinyjs::reset("Bone_Marrow")
  shinyjs::reset("Breast")
  shinyjs::reset("Cervix")
  shinyjs::reset("other_tissue")
  cancer_type <- c("")
  # cancer_type <- reactive({
  #   c("") -> cancer_type
  # })
  return(cancer_type)
}



# remove pic when stop clicked ----------------------------------------------


removePic <- function(input, output, session, outtype) {
  if (outtype == "image") {
    output$plot <- renderImage({
      NULL
    })
  }
  if (outtype == "plot") {
    output$plot <- renderPlot({
      NULL
    })
  }
}


###############################################################
# Plot function to generate plot in ui#########################
## 1. plotoutout in ui and in server ###########################
## PlotInput & Plot#################################
###############################################################
# call in ui by PlotInput("cnv_pie",..) OR  PlotInput("cnv_bar",..)
# call in ser by callModule(Plot,"cnv_pie",...) OR ...

download_bt <- function(id){
  ns <- NS(id)
  tagList(
    shinyWidgets::dropdownButton(
      tags$h3("Download Options"),
      prettyRadioButtons(
        inputId = ns("pictype"),
        label = "Selcet format for your pictur",
        choices = list("PDF" = "pdf", "PNG" = "png","EPS"="eps"),
        inline = TRUE,
        icon = icon("check"),
        bigger = TRUE, status = "info",
        animation = "jelly"
      ),
      numericInput(
        inputId = ns("d_width"),
        label = "Width",
        value = 4,
        min = 1,
        max = 10
      ),
      
      numericInput(
        inputId = ns("d_height"),
        label = "Height",
        value = 6,
        min = 3,
        max = 20
      ),
      downloadButton(
        outputId = ns("picdownload"),
        label = "Download"
      ),
      circle = TRUE, status = "default",
      right = TRUE,
      icon = icon("download"), width = "300px",
      tooltip = shinyWidgets::tooltipOptions(title = "Click to download")
    )
  )
}
PlotInput <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("massage")),
    column(
      width = 2, offset = 0,
      download_bt(id)
    ),
    column(
      width = 12, offset = 0,
      plotOutput(ns("plot"), height = "100%") %>% withSpinner(color = "#0dc5c1",size = 0.5, proxy.height = "200px")
    )
  )
}

##################### GTEx expression heatmap plot by zhangq#########################

heatmap_GTEX_Plot <- function(input, output, session, data, status, downloadname) {
  # plot function
  plotinput <- function(){
    ggplot(data, aes(Tissue, GeneName)) +
      geom_tile(aes(fill = TPM)) +
      geom_text(aes(label = TPM)) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(title = "Expression value of query genes in GTEx dataset") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) -> p
  }
  
  # get output
  plot_height <- data$GeneName %>% unique() %>% length()*20
  output$plot <- renderPlot({
    status$analysis
    print(plotinput())
  } ,height = function(){ifelse(plot_height<200,200,plot_height)})
  
  # get download output
  output$picdownload <-downloadHandler(
    filename = function() {
      paste(downloadname, ".", input$pictype, sep = "")
    },
    content = function(file){
      ggsave(file,plotinput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}

box_GTEx_GSVA_Plot <- function(input, output, session, data) {
  output$plot <- renderPlot({
    data %>%
      dplyr::group_by(SMTS) %>%
      dplyr::summarise(m = median(gsva)) %>%
      dplyr::arrange(m) %>%
      dplyr::pull(SMTS) -> lev
    tcc <- tibble(SMTS = lev, color = rainbow(length(lev)))
    data %>%
      dplyr::mutate(SMTS = factor(SMTS, levels = lev)) %>%
      ggplot(aes(x = SMTS, y = gsva)) +
      stat_boxplot(geom = "errorbar", width = 0.3) +
      geom_boxplot(outlier.colour = NA) +
      geom_point(aes(color = SMTS), position = position_jitter(width = 0.05), alpha = 0.4, size = 0.8) +
      scale_color_manual(name = "Tissues", values = dplyr::slice(tcc, match(lev, SMTS)) %>% dplyr::pull(color)) +
      theme(
        axis.line = element_line(color = "black"), axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        panel.grid = element_blank(), panel.background = element_rect(fill = "white", color = NA), panel.spacing.x = unit(0, "lines")
      ) +
      guides(color = F) + labs(x = "Tissue", y = "GSVA Score", title = "") -> p
    return(p)
  })
}




Plot <- function(input, output, session) { # data(raw data, gene set, cancer types),plot type(decide function type)
  # size <- reactive(as.numeric(input$num))
  # x <- reactive({
  #   input$num %>% as.numeric %>% rnorm()
  #   })
  # y <- reactive({
  #   input$num %>% as.numeric %>% rnorm()
  # })
  x <- rnorm(50)
  y <- nrow(50)

  output$plot <- renderPlot({
    plot(x, y)
  }) # fun argument decide what function will be called.
}

# cnv Point plot --------------------------------------------------------------


cnv_pointPlot <- function(input, output, session, data, cancer, gene, size, color, sizename, colorname, wrap, status_monitor, status,downloadname) {

  # Example: callModule(pointPlot,"cnv_pie",data=cnv_plot_ready_1,cancer="cancer_types",
  #                     gene="symbol",size="per",color="color",sizename="CNV%",
  #                     colorname="SCNA Type",wrap="~ effect")
  # data should include x/y, point size and point color.
  plotinput <- reactive({
    data %>%
      ggplot(aes_string(y = gene, x = cancer)) +
      geom_point(aes_string(size = size, color = color)) +
      xlab("Cancer type") +
      ylab("Symbol") +
      scale_size_continuous(
        name = sizename,
        breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 1),
        limits = c(0.05, 1),
        labels = c("5", "10", "20", "40", "60", "100")
      ) +
      ggthemes::scale_color_gdocs(
        name = colorname,
        labels = c("Deletion", "Amplification")
      ) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      facet_wrap(as.formula(wrap)) +
      theme(strip.text.x = element_text(size = 15)) -> p
  })
  
  plot_height <- data$symbol %>% unique() %>% length()*20
  
  output$plot <- renderPlot({
    status[[status_monitor]]
    print(plotinput())
  } ,height = function(){ifelse(plot_height<200,200,plot_height)})
  output$picdownload <-downloadHandler(
    filename = function() {
      paste(downloadname, ".", input$pictype, sep = "")
    },
    content = function(file){
      ggsave(file,plotinput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}


# cnv pie plot ----------------------------------------------------------------

piePlot <- function(input, output, session, data, y, fill, facet_grid, outfile, height, width, status_monitor, status,downloadname) {
  # Example:
  # callModule(piePlot,"cnv_pie",data=pie_plot_ready,y="per",
  #            fill="type",facet_grid="cancer_types ~ symbol")
  # data should include ...

  imgplotInput <- reactive({
    data %>%
      ggplot(aes_string(x = factor(1), y = y, fill = fill)) +
      geom_bar(stat = "identity", position = "stack", color = NA) +
      # scale_y_continuous(limits = c(0,1))
      coord_polar("y") +
      facet_grid(as.formula(facet_grid)) + # cancer_types ~ symbol
      # scale_x_discrete(limits = cnv_gene_rank$symbol) +
      # scale_x_discrete(expand=c(0,0)) +
      # scale_y_discrete(expand=c(0,0)) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),

        strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
        strip.text.x = element_text(size = 4, angle = 90, vjust = 0),
        strip.background = element_blank(),

        legend.title = element_blank(),
        legend.text = element_text(size = 4),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),

        panel.background = element_blank(),
        panel.spacing = unit(0, "null"), # unit(0.01, "lines"),
        panel.spacing.x = unit(0, "null"),

        plot.margin = rep(unit(0, "null"), 4),
        axis.ticks.length = unit(0, "cm")
      ) +
      scale_fill_manual(
        limits = c("a_hete", "a_homo", "d_hete", "d_homo", "other"),
        label = c("Hete Amp", "Homo Amp", "Hete Del", "Homo Del", "None"),
        # Amp RColorBrewer name = "Spectral"
        # Del RColorBrewer name = "BrBG"
        values = c("brown1", "brown4", "aquamarine3", "aquamarine4", "grey")
      ) -> p
  })
  
  output$plot <- renderImage({
    status[[status_monitor]]
    # outfile <- paste("/project/huff/huff/github/GSCALite/userdata","/","TCGA_cnv_pie_rellation_network",'.png',sep="")
    ggsave(outfile, imgplotInput(), device = "png", width = width, height = height)
    list(
      src = outfile,
      contentType = "image/png",
      # width = 400,
      # height = "900px",
      alt = "This is alternate text"
    )
  }, deleteFile = FALSE)
  
  output$picdownload <- downloadHandler(
    filename = function() {
      paste(downloadname, ".", input$pictype, sep = "")
    },
    content = function(file) {
      ggsave(file, imgplotInput(), device = input$pictype, width = input$d_width, height = input$d_height)
    }
  )
}


# gene set CNV frenquencey in each cancer ---------------------------------
# bar stak plot
cnvbarPlot <- function(input, output, session, data, x, y, fill, status_monitor, status, downloadname) {
  # Example:
  # callModule(piePlot,"cnv_pie",data=pie_plot_ready,y="per",
  #            fill="type",facet_grid="cancer_types ~ symbol")
  # data should include ...
  plotinput <- reactive({
    data %>%
      ggplot(aes_string(x = x, y = y, fill = fill)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~cnv_type, ncol = 2) +
      theme(strip.text.x = element_text(size = 15)) +
      ggsci::scale_fill_npg(
        name = "Type",
        limits = c("amp_a", "amp_s", "del_a", "del_s"),
        labels = c("Amp", "Amp Only", "Del", "Del Only")
      ) +
      ggthemes::theme_gdocs() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)
      ) +
      labs(x = "Cancer Types", y = "CNV Frequency") -> p
  })
  
  # plot_height <- data$symbol %>% unique() %>% length()*20
  
  output$plot <- renderPlot({
    status[[status_monitor]]
    print(plotinput())
  } )
  
  output$picdownload <- downloadHandler(
    filename = function(){
      paste(downloadname,".",input$pictype, sep='')
    },
    content = function(file) {
      ggsave(file,plotinput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}


# snv percentage plot -----------------------------------------------------

snv_per_heatmap <- function(input, output, session, data, cancer, gene, fill, label, cancer_rank, gene_rank, status_monitor, status, downloadname) {
  #per_max <- data$per %>% max()
  
  plotInput <- reactive({
    data %>%
      ggplot(aes_string(x = cancer, y = gene, fill = fill)) +
      geom_tile() +
      geom_text(aes_string(label = label)) +
      scale_x_discrete(position = "top", limits = cancer_rank$x_label) +
      scale_y_discrete(limits = gene_rank$symbol) +
      scale_fill_gradient2(
        name = "Mutation Frequency (%)",
        limit = c(0, 1),
        breaks = c(seq(0, 0.2, 0.05), seq(0.3, 1, 0.1)),
        label = c("0", "5", "10", "15", "20", "30", "40", "50", "60", "70", "80", "90", "100"),
        high = "red",
        na.value = "white"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = -0.05, size = "10"),
        axis.title.y = element_text(size = "15"),
        panel.grid = element_line(colour = "grey", linetype = "dashed")
      ) +
      guides(fill = guide_legend(
        title = "Mutation Frequency (%)",
        title.position = "left",
        title.theme = element_text(angle = 90, vjust = 2),
        reverse = T,
        keywidth = 0.6,
        keyheight = 0.8
      )) +
      labs(x = "", y = "") -> p
  })
  plot_height <- data$symbol %>% unique() %>% length()*20
  
  output$plot <- renderPlot({
    # data$per %>% max() ->max.limit
    # max.limit/10 -> inter.limit
    status[[status_monitor]]
    print(plotInput())
  } ,height = function(){ifelse(plot_height<200,200,plot_height)})
  
  output$picdownload <- downloadHandler(
    filename = function() { paste(downloadname, '.',input$pictype, sep='') },
    content = function(file) {
      ggsave(file,plotInput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}


# snv survival point plot -------------------------------------------------

snv_sur_pointPlot <- function(input, output, session, data, cancer, gene, size, color, cancer_rank, gene_rank, sizename, colorname, title, status_monitor, status, downloadname) {
  # Example: callModule(pointPlot,"cnv_pie",data=cnv_plot_ready_1,cancer="cancer_types",
  #                     gene="symbol",size="per",color="color",sizename="CNV%",
  #                     colorname="SCNA Type",wrap="~ effect")
  # data should include x/y, point size and point color.
  plotInput <- reactive({
    data %>%
      ggplot(aes_string(y = gene, x = cancer)) +
      geom_point(aes_string(size = size, color = color)) +
      labs(title = title) +
      xlab("Cancer type") +
      ylab("Symbol") +
      scale_x_discrete(limit = cancer_rank$cancer_types) +
      scale_y_discrete(limit = gene_rank$symbol) +
      scale_size_continuous(
        name = sizename,
        breaks = c(-log10(0.05), 5, 10, 15),
        limits = c(-log10(0.05), 15),
        labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$"))
      ) +
      #ggthemes::scale_color_gdocs(
      scale_color_manual(
        name = colorname,
        labels = c("High", "Low"),
        values = c("#e31a1c", "#1f78b4")
      ) +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        ),
        plot.title = element_text(size = 20),
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
      ) -> p
  })
  
  plot_height <- data$symbol %>% unique() %>% length()*20
  
  output$plot <- renderPlot({
    status[[status_monitor]]
    print(plotInput())
  } ,height = function(){ifelse(plot_height<200,200,plot_height)})
  
  output$picdownload <- downloadHandler(
    filename = function() { paste(downloadname, '.',input$pictype, sep='') },
    content = function(file) {
      ggsave(file,plotInput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}



# snv maf summary ---------------------------------------------------------

# 1. ui part -----------------------------------------------------------------

imagePlotInput <- function(id, width="100%", height=300) {
  ns <- NS(id)
  shiny::tagList(
    uiOutput(ns("massage")),
    column(
      width = 2, offset = 0,
      download_bt(id)
    ),
    br(),
    br(),
    column(
      width = 12, offset = 0,
      imageOutput(ns("plot"), width = width, height = height) %>% withSpinner(color = "#0dc5c1",size = 0.5, proxy.height = "200px")
    )
  )
}

# 2. server part ----------------------------------------------------------

snv_maf_summaryPlot <- function(input, output, session, gene_list_maf, outfile, status_monitor, status, downloadname) {
  plotInput <- reactive({
    maftools::plotmafSummary(gene_list_maf, fs = 3) 
    #p
  })
  output$plot <- renderImage({
    status[[status_monitor]]
    # png(outfile, width = 1000, height= 700)
    
    # dev.off()
    ggsave(plotInput(), filename = outfile, device = "png", width = 3, height = 2)
    list(
      src = outfile,
      contentType = "image/png",
      # width = 1000,
      # height = 700,
      alt = "This is alternate text"
    )
  }, deleteFile = FALSE)
  
  output$picdownload <- downloadHandler(
    filename = function() { paste(downloadname, '.',input$pictype, sep='') },
    content = function(file) {
      ggsave(file,plotInput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}

snv_maf_oncoPlot <- function(input, output, session, gene_list_maf, pancan_color, outfile, status_monitor, status, downloadname) {
  plotInput <- function(){
    col <- RColorBrewer::brewer.pal(n = 8, name = "Paired")
    names(col) <- c(
      "Frame_Shift_Del", "Missense_Mutation", "Nonsense_Mutation", "Multi_Hit", "Frame_Shift_Ins",
      "In_Frame_Ins", "Splice_Site", "In_Frame_Del"
    )
    gene_list_maf %>% maftools::getClinicalData() %>% dplyr::select(Cancer_Types) %>% unique() %>% t() %>% as.character() -> snv_maf_cancer_type
    pancan_color %>%
      dplyr::filter(cancer_types %in% snv_maf_cancer_type) %>%
      dplyr::select(color, cancer_types) -> snv_maf_cancer_type_color
    
    fabcolors <- snv_maf_cancer_type_color$color
    names(fabcolors) <- snv_maf_cancer_type_color$cancer_types
    
    fabcolors <- list(Cancer_Types = fabcolors)
    if (length(snv_maf_cancer_type) > 1) {
      maftools::oncoplot(
        # my_oncoplot(
        maf = gene_list_maf, removeNonMutated = T, colors = col,
        clinicalFeatures = "Cancer_Types", sortByMutation = TRUE, sortByAnnotation = TRUE,
        annotationColor = fabcolors, top = 10
      )
    } else {
      maftools::oncoplot(
        # my_oncoplot(
        maf = gene_list_maf, removeNonMutated = T, colors = col,
        clinicalFeatures = "Cancer_Types", sortByMutation = TRUE, # sortByAnnotation = TRUE,
        annotationColor = fabcolors, top = 10
      )
    }
  }
  output$plot <- renderImage({
    status[[status_monitor]]
    png(outfile, width = 800, height = 600)
    plotInput()
    # maftools::oncoplot(maf = gene_list_maf, top = 10)#, fontSize = 12
    dev.off()
    list(
      src = outfile,
      contentType = "image/png",
      width = 800,
      height = 600,
      alt = "This is alternate text"
    )
  }, deleteFile = FALSE)
  
  fn_save <- function(.file,.pictyp,width,height){
    print(.file)
    if(.pictyp == "png"){
      png(.file, width, height,units ="in",res=1200)
      plotInput()
      dev.off()
    } else{
      pdf(.file, width, height)
      plotInput()
      dev.off()
    }
  }
  output$picdownload <- downloadHandler(
    filename = function() {paste(downloadname, '.',input$pictype, sep='') },
    content = function(filename) {
      fn_save(.file=filename, .pictyp=input$pictype, width = input$d_width,height = input$d_height)
    }
  )
}


# methylation plot --------------------------------------------------------


# 1. methy diff -----------------------------------------------------------
methy_diff_pointPlot <- function(input, output, session, data, cancer, gene, size, color, cancer_rank, gene_rank, sizename, colorname, title, status_monitor, status, downloadname) {
  plotinput <- reactive({
    CPCOLS <- c("red", "white", "blue")
    data %>%
      ggplot(aes_string(y = gene, x = cancer)) +
      geom_point(aes_string(size = size, color = color)) +
      scale_y_discrete(limit = gene_rank$symbol) +
      scale_x_discrete(limit = cancer_rank$cancer_types) +
      labs(title = title) +
      ylab("Symbol") +
      xlab("Cancer types") +
      scale_size_continuous(
        name = sizename # "-Log10(FDR)"
      ) +
      scale_color_gradient2(
        name = colorname, # "Methylation diff (T - N)",
        low = CPCOLS[3],
        mid = CPCOLS[2],
        high = CPCOLS[1]
      ) +
      theme(
        legend.position = "bottom",
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        ),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(vjust = 1, hjust = 1, angle = 40, size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(size = 20)
      ) -> p
  })
  
  plot_height <- data$symbol %>% unique() %>% length()*20
  
  output$plot <- renderPlot({
    status[[status_monitor]]
    print(plotinput())
  },height = function(){ifelse(plot_height<200,200,plot_height)})
  
  output$picdownload <- downloadHandler(
    filename = function() { paste(downloadname, '.',input$pictype, sep='') },
    content = function(file) {
      ggsave(file,plotinput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}




# rppa --------------------------------------------------------------------
# line contact ----
rppa_line_contact <- function(plot_seg, cancer.text, gene.text, path.text) {
  
  ggplot() -> p
  for (cancers in plot_seg$Cancer %>% unique()) {
    # cancers="LUSC"
    plot_seg %>%
      dplyr::filter(Cancer == cancers) -> data
    curvature <- runif(1, 0.1, 0.3)
    p +
      geom_curve(
        data = data, mapping = aes(
          x = x1,
          y = y1,
          xend = x2,
          yend = y2,
          colour = Cancer,
          linetype = Regulation
        ),
        # colour = "red",
        curvature = curvature
      ) -> p
  }
  
  p +
    guides(color = FALSE) +
    geom_text(
      data = cancer.text,
      mapping = aes(x = x, y = y, label = text, color = text),
      hjust = 1,
      size = 2
    ) +
    geom_text(
      data = gene.text,
      mapping = aes(x = x - 0.4, y = y, label = text),
      hjust = 0,
      size = 2
    ) +
    geom_text(
      data = path.text,
      mapping = aes(x = x, y = y, label = text),
      hjust = 0,
      size = 2
    ) +
    expand_limits(x = c(-1, 10)) +
    theme(
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      # text = element_text(size=5),
      plot.title = element_text(hjust = 0.5, size = 7),
      plot.margin = rep(unit(0, "null"), 4),
      legend.position = "bottom",
      legend.text = element_text(size = 3),
      legend.key.size = unit(0.25, "cm"),
      legend.title = element_text(size = 4)
    ) +
    xlab("") +
    ylab("") +
    labs(title = "Interaction map of gene and pathway.") -> p
}

# rppa pie ----
rppaPiePlot <- function(input, output, session, data, y, fill, facet_grid, height, outfile, status, downloadname) {
  # Example:
  # callModule(piePlot,"cnv_pie",data=pie_plot_ready,y="per",
  #            fill="type",facet_grid="cancer_types ~ symbol")
  # data should include ...

  imgInput <- function(){
    data %>%
      ggplot(aes_string(x = factor(1), y = y, fill = fill)) +
      geom_bar(stat = "identity", position = "stack", color = NA) +
      # scale_y_continuous(limits = c(0,1))
      coord_polar("y") +
      facet_grid(as.formula(facet_grid)) + #  symbol~ cancer_types
      # scale_x_discrete(limits = cnv_gene_rank$symbol) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),

        strip.text.y = element_text(angle = 0, hjust = 0, size = 5),
        strip.text.x = element_text(size = 5, angle = 90, vjust = 0),
        strip.background = element_blank(),

        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom",

        panel.background = element_blank(),
        panel.spacing = unit(0.02, "lines"),
        plot.margin = rep(unit(0, "null"), 4),
        axis.ticks.length = unit(0, "cm")
      ) +
      scale_fill_manual(
        limits = c("Activation", "Inhibition", "None"),
        label = c("Activation", "Inhibition", "None"),
        # Amp RColorBrewer name = "Spectral"
        # Del RColorBrewer name = "BrBG"
        values = c("brown1", "aquamarine3", "grey")
      ) -> p
  }
  
  output$plot <- renderImage({
    status$analysis
    ggsave(outfile, imgInput(), device = "png", width = 3, height = height)
    list(
      src = outfile,
      contentType = "image/png",
      # width = "100%" ,
      # height = 900,
      alt = "This is alternate text"
    )
  }, deleteFile = TRUE)
  
  output$picdownload <- downloadHandler(
    filename = function() {
      paste(downloadname, ".", input$pictype, sep = "")
    },
    content = function(file) {
      ggsave(file, imgInput(), device = input$pictype, width = input$d_width, height = input$d_height)
    }
  )
}

# rppa heatmap percent ----
rppa_heat_per <- function(input, output, session, rppa_per_ready, pathway, symbol, per, height, outfile, status, downloadname) {
  plotInput <- function(){
    rppa_per_ready %>%
      ggplot(aes(x = pathway, y = symbol)) +
      xlab("Pathway") + ylab("Symbol") +
      guides(fill = guide_colorbar("Percent")) +
      geom_tile(aes(fill = per), col = "white") +
      geom_text(
        label = ceiling(rppa_per_ready$per)
        # size = 1
      ) +
      scale_fill_gradient2(
        high = "red",
        mid = "white",
        low = "blue"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black",size = 10),
        axis.title = element_text(size = 13),
        # legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom",
        plot.margin = rep(unit(0, "null"), 4),
        axis.ticks.length = unit(0, "cm"),
        # legend.text = element_text(size = 5),
        # axis.title.x = element_text(size = 6),
        # axis.title.y = element_text(size = 6),
        # legend.title = element_text(size = 6),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(
          colour = "grey",
          linetype = "dashed",
          size = 0.2
        )
      ) +
      xlab("Pathway (A:Activate; I:Inhibit)") -> p
  }
  
  plot_height <- rppa_per_ready$symbol %>% unique() %>% length()*20
  output$plot <- renderPlot({
    status$analysis
    print(plotInput())
  },height = function(){ifelse(plot_height<200,200,plot_height)})
  
  # output$plot <- renderImage({
  #   status$analysis
  #   ggsave(outfile, plotInput(), device = "png", width = 4, height = height)
  #   list(
  #     src = outfile,
  #     contentType = "image/png",
  #     # width = "100%" ,
  #     # height = 900,
  #     alt = "This is alternate text"
  #   )
  # }, deleteFile = FALSE)
  
  output$picdownload <- downloadHandler(
    filename = function() { paste(downloadname, '.',input$pictype, sep='') },
    content = function(file) {
      ggsave(file,plotInput(),device = input$pictype,width = input$d_width,height = input$d_height)
    }
  )
}





#### GTEx eqtl table output-------------

GTEx_eqtl_Output <- function(id) {
  ns <- NS(id)
  column(
    width = 10, offset = 1,
    shinydashboard::tabBox(
      id = "gtex_eqtl_table", title = "TABLE", width = 12,
      # datatable
      tabPanel(
        title = "Table of eQTL in GTEX dataset",
        DT::dataTableOutput(outputId = ns("gtex_eqtl_dt"))
      )
    )
  )
}


# white plot generate -----------------------------------------------------

white_plot <- function(input, output, session, status_monitor, status, outfile){
  data <- data.frame(x=c(1:100),y=c(1:100))
  plotInput<- function(){
      data %>%
      ggplot(aes(x=x,y=y)) +
      theme(
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )}
  output$plot <- renderImage({
    status$analysis
    ggsave(plotInput(), filename = outfile, device = "png", width = 3, height = 2)
    list(
      src = outfile,
      contentType = "image/png",
      # width = 1000,
      # height = 700,
      alt = "This is alternate text"
    )
  }, deleteFile = FALSE)
}
