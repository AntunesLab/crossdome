library(shiny)

hla_valid <- sort(unique(crossdome::hla_database$hla_allele))

ui <- navbarPage("Crossdome App",
                tabPanel("Home", fluid = TRUE,
                   sidebarLayout(
                       sidebarPanel(
                           textInput("peptide", h3("Query peptide:"), value = "EVDPIGHLY"),
                           selectInput("allele", h3("HLA I allele:"),
                                       choices = hla_valid,
                                       selected = 1),
                           radioButtons("weights",
                                              h3("Allele input:"),
                                              choices = list(
                                                  "Biochemical properties only " = 1,
                                                  "Crystal-derived TCR Weight" = 2,
                                                  "Predicted TCR Weight" = 3
                                                  ),
                                              selected = 1),
                           actionButton("submit", "Submit", class = "btn-primary"),
                           width = 3
                       ),
                       mainPanel(
                           textOutput('input_text'),
                           hr(),
                           DT::dataTableOutput('table'),
                           #downloadButton('download',"Download the data"),
                           width = 8
                       )
                    )
                ),
                tabPanel("Analyze", verbatimTextOutput('selected')),
                tabPanel("Help", includeHTML("help.html")),
                tabPanel("About", includeHTML("about.html"))
)

server <- function(input, output) {

    crossdomeParameters <- eventReactive(input$submit, {
        paste(input$peptide, input$allele, input$weights)
    })

    crossdomeInput <- eventReactive(input$submit, {
        if(TRUE) {
            position_weight <- list(
                `1` = rep(1, 9),
                `2` = c(0.33, 0.31, 0.33, 0.66, 0.21, 0.25, 0.12, 0.18, 0.33),
                `3` = rep(1, 9)
            )

            position_weight <- position_weight[[input$weights]]
            cross_background <- crossdome::cross_universe(subject = c(), allele = input$allele)

            cross_result <- data.frame()

            withProgress(
                message = "Calculating...", value = 0, {
                    cross_result <- crossdome::cross_compose(
                        query = input$peptide,
                        subject = cross_background$subject,
                        allele = input$allele,
                        position_weight = position_weight
                    )
                }
            )
        }
        dplyr::mutate_if(cross_result, is.numeric, round, 2)
    })

    output$input_text <- renderText({
        crossdomeParameters()
    })

    output$table <- DT::renderDataTable(
        server = FALSE,
        crossdomeInput(),
        rownames = FALSE,
        extensions = c('Buttons'),
        filter = list(
            position = 'top',
            clear = FALSE
        ),
        options = list(
            dom = 'Bfrtip',
            pageLength = 50,
            buttons = list(
                list(
                    extend = 'csv', text = "Download Current Page",
                    filename = paste('crossdome.', Sys.Date(), '.filtered', sep=''),
                    exportOptions = list(
                        modifier = list(page = "current")
                    )
                ),
                list(
                    extend = 'csv', text = "Download All Results",
                    filename = paste('crossdome.', Sys.Date(), '.all', sep=''),
                    exportOptions = list(
                        modifier = list(page = "all")
                    )
                )
            )
        )
    )

    output$selected <- renderPrint({
        selected <- input$table_rows_selected
        if(selected) {
            cat(selected, sep = ",")
        }
    })

    output$download <- downloadHandler(
        filename = function() {
            paste('crossdome.', Sys.Date(), '.tsv', sep='')
        },
        content = function(file) {
            write.table(crossdomeInput(), file, sep = "\t", row.names = FALSE)
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
