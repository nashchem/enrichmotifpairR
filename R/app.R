#' Opens interactive shiny to run enrichmotifpairR.
#' @examples
#' \dontrun{
#' runEnrichMotifPair()
#' }
#' @import shiny
#' @import DT
#' @export
#' 

runEnrichMotifPair <- function(...){
    data("example_peaks_data")
    # Define UI 
    ui <- fluidPage(
        
        # Application title
        titlePanel("EnrichMotifPair"),
        h5("This is a tool, built on the enrichmotifpairR library, for discovering enriched TF motifs and their binding partner motifs
       in a set of genomic regions."),
        h5("A optional set of background regions can be provided as control. If none are provided,
       EnrichMotifPair will generate background regions based on the target inputs."),
        h5("Heatmaps are generated to visualize the degree of enrichment of select motifs, and networks display the range of TF-motif interactions."),
        h5("A set of example DHS peaks from ESCs are provided. Use these as inputs by selecting 'Use Example Dataset'."),
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Inputs
            checkboxInput("example", "Use example dataset", FALSE),
            fileInput("target_data", 
                      "Upload target peaks file (BED)",
                      multiple = FALSE,
                      accept = c("text/tsv", ".tsv")),
            checkboxInput("background", "Use own background peaks?", FALSE),
            selectInput('genome', 'Select genome', choices = c('hg38', 'hg19')),
            selectInput('motif_db', 'Select motif database', choices = c("CISBP", "JASPAR", "ENCODE", "HOMER", "HOCOMOCCO", "JASPAR_UNVALIDATED")),
            selectInput('distribution', 'Select distribution test', choices = c('binom', 'hyper')),
            selectInput('adjustment', 'Select adjustment method', choices = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none")),
            numericInput('pval_thresh', 'Set p-value threshold', value = 0.01),
            actionButton(
                inputId = "submit_loc",
                label = "Submit"
            ),
            
        ),
        
        
        
        #downloadButton("downloadData", "Download"),
        uiOutput("get_the_item1"),
        uiOutput("get_the_item2"),
        
        mainPanel(
            tabsetPanel(
                tabPanel("Tables", 
                         h4("Enriched Motifs"),
                         DT::dataTableOutput("contents1"),
                         h4("Enriched Motif Pairs"),
                         DT::dataTableOutput("contents2")
                ),
                tabPanel("Heatmap", 
                         h5("Visualize top enriched binding partners for a set of TFs"),
                         uiOutput("get_heatmap"),
                         plotOutput("enrichmentpair"),
                         textAreaInput("tfs", "Input TFs (one per line)"), 
                         actionButton(
                             inputId = "submit_heatmap",
                             label = "Go"
                         ),
                         textOutput("enrichment")
                ),
                
                tabPanel("Network", 
                         h5("Visualize a network of all binding partners for a TF"),
                         uiOutput("get_network"),
                         plotOutput("enrichmentnetwork"),
                         textAreaInput("tf", "Input TF"),
                         actionButton(
                             inputId = "submit_network",
                             label = "Go"
                         )
                )
            ),
        ))
    
    
    
    
    # Define server logic
    server <- function(input, output){
        observeEvent({input$example},
                     if(input$example == T){
                         removeUI(
                             selector = "div:has(> #background)",
                             immediate = T,
                             multiple = T
                         )
                         
                     }
        )
        
        observeEvent({input$background},
                     if(input$background == T){
                         insertUI(
                             selector = "#background",
                             where = "afterEnd",
                             ui = fileInput("background_data",
                                            "Upload background BED peaks file",
                                            multiple = FALSE,
                                            accept = c("text/tsv", ".tsv")
                             )
                         )
                     }
        )
        
        run_motifs = reactiveValues(tmp=NULL)
        enrich_motifs <- eventReactive(input[["submit_loc"]], {
            
            # Create a Progress object
            progress <- shiny::Progress$new()
            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress$close())
            prog = 0
            progress$set(message = "Initializing", value = prog)
            
            if (input$example == F){
                df_target <- data.table::fread(input$target_data$datapath)
                colnames(df_target) = c("chr", "start", "end")
                
                
                
                if (input$background == T){
                    df_background <- data.table::fread(input$background_data$datapath)
                    colnames(df_background) = c("chr", "start", "end")
                } else {
                    progress$set(message = "Generating background peaks", value = 0)
                    df_background = NULL
                    prog <- prog + 0.5
                }
            }
            
            if (input$example == T){
                df_target = example_peaks_data$`H1-ESC_DHS_peaks`
                df_background = example_peaks_data$`H1-ESC_DHS_peaks_matched_background`
                prog <- prog + 0.5
            }
            
            Sys.sleep(5.0)
            
            progress$set(message = "Finding motif pairs", value = prog)
            motifs = list()
            results <- findEnrichMotifPair(
                target_data = df_target,
                background_data = df_background,
                genome_ver = input$genome,
                scramble_data = FALSE,
                motif_database = input$motif_db,
                Pvalue_computation = input$distribution,
                Pvalue_threshold = input$pval_thresh,
                Pvalue_adjust_method = input$adjustment
            )
            motifs[[1]] <- results$motif_enrich
            motifs[[2]] <- results$motif_pair_enrich
            run_motifs$tmp = T
            return(motifs)
        })
        
        #results = reactive({enrich_motifs()})
        
        output$contents1 = DT::renderDataTable({
            enrich_motifs()[[1]]
        })
        
        output$contents2 = DT::renderDataTable({
            enrich_motifs()[[2]]
        })
        
        
        
        plotHeatmapInput <- eventReactive(input[["submit_heatmap"]], {
            
            input_tfs <- input$tfs
            input_tfs_split <- unlist(strsplit(input_tfs, "\n"))
            plotEnrichPair(enrich_motifs()[[2]], input_tfs_split)
        })
        
        
        output$enrichmentpair = renderPlot({
            
            print(plotHeatmapInput())
        }
        )
        
        plotNetworkInput = eventReactive(input[["submit_network"]], {
            
            plotNetwork(enrich_motifs()[[2]], TF_name = input$tf, 
                        color_TF = "#70d9e0", color_bind_TF = "#e841da")
        }
        )
        
        output$enrichmentnetwork = renderPlot({
            
            print(plotNetworkInput())
        }
        )
        
        
        
        output$get_the_item1 <- renderUI({
            req(run_motifs$tmp)
            downloadButton('downloadData1', label = 'Enriched Motifs')})
        
        output$get_the_item2 <- renderUI({
            req(run_motifs$tmp)
            downloadButton('downloadData2', label = 'Enriched Pairs')})
        
        output$get_heatmap <- renderUI({
            req(run_motifs$tmp)
            downloadButton('downloadHeatmap', label = 'Download Image')})
        
        output$get_network <- renderUI({
            req(run_motifs$tmp)
            downloadButton('downloadNetwork', label = 'Download Image')})
        
        output$downloadData1 <- downloadHandler(
            filename = function() {
                paste(input$target_data, "_enriched_motifs.tsv", sep = "")
            },
            content = function(file) {
                write.table(enrich_motifs()[[1]], file, row.names = F, col.names = T, quote = F, sep = "\t")
            }
        )
        
        output$downloadData2 <- downloadHandler(
            filename = function() {
                paste(input$target_data, "_enriched_motif_pairs.tsv", sep = "")
            },
            content = function(file) {
                write.table(enrich_motifs()[[2]], file, row.names = F, col.names = T, quote = F, sep = "\t")
            }
        )
        
        output$downloadHeatmap <- downloadHandler(
            filename = function() {
                paste(input$target_data, "enrichMotifPair_heatmap.pdf", sep = "")
            },
            content = function(file) {
                pdf(file)
                print(plotHeatmapInput())
                dev.off()
            }
        )
        
        output$downloadNetwork <- downloadHandler(
            filename = "enrichMotifPair_network.pdf"
            ,
            content = function(file) {
                
                # device <- function(..., width, height) {
                #     grDevices::png(..., width = width, height = height,
                #                    res = 300, units = "in")
                # }
                
                pdf(file)
                print(plotNetworkInput())
                dev.off()
                
            }
        )
    }
    
    # Run the application 
    shinyApp(ui = ui, server = server) 
}



