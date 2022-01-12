library(shiny)
library(enrichmotifpairR)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("EnrichMotifPair"),
    h5("This is a tool, built on the enrichmotifpairR library, for discovering enriched TF motifs and their binding partner motifs
       in a set of genomic regions."),
    h5("A optional set of background regions can be provided as control. If none are provided,
       EnrichMotifPair will generate background regions based on the target inputs."),
    h5("Heatmaps are generated to visualize the degree of enrichment of select motifs, and networks display the range of TF-motif interactions"),
    
    # Sidebar panel for inputs ----
    sidebarPanel(
        
        # Inputs
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
                     plotOutput("enrichmentpair"),
                     textAreaInput("tfs", "Input TFs (one per line)"), 
                     textOutput("enrichment")
                     ),
                    
            tabPanel("Network", 
                     h5("Visualize a network of all binding partners for a TF"),
                     plotOutput("enrichmentnetwork"),
                     textAreaInput("tf", "Input TF")
                     )
            ),
    ))
    



# Define server logic required to draw a histogram
server <- function(input, output){
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
        df_target <- data.table::fread(input$target_data$datapath)
        colnames(df_target) = c("chr", "start", "end")
        
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        prog = 0
        progress$set(message = "Initializing", value = prog)
        
        if (input$background == T){
            df_background <- data.table::fread(input$background_data$datapath)
            colnames(df_background) = c("chr", "start", "end")
        } else {
            progress$set(message = "Generating background peaks", value = 0)
            df_background = NULL
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
    
    # output$enrichment = renderPlot({
    #     plotEnrichment(enrich_motifs()[[1]], input$tfs)
    #     }
    # )
    
    output$enrichmentpair = renderPlot({
        req(input$tfs)
        input_tfs <- input$tfs
        input_tfs_split <- unlist(strsplit(input_tfs, "\n"))
        plotEnrichPair(enrich_pairs, input_tfs_split)
    }
    )
    
    output$enrichmentnetwork = renderPlot({
        req(input$tf)
        plotNetwork(enrich_motifs()[[2]], TF_name = input$tf, 
                    color_TF = "#70d9e0", color_bind_TF = "#e841da")
    }
    )
    
    
    #output$enrichment = renderText({input$tfs})
   
   
    output$get_the_item1 <- renderUI({
        req(run_motifs$tmp)
        downloadButton('downloadData1', label = 'Enriched Motifs')})
    
    output$get_the_item2 <- renderUI({
        req(run_motifs$tmp)
        downloadButton('downloadData2', label = 'Enriched Pairs')})
    
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
    
}

# Run the application 
shinyApp(ui = ui, server = server)

