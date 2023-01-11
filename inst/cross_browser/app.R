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
                           radioButtons("protocol",
                                              h3("Protocol/Mode:"),
                                              choices = list(
                                                  "Biochemical properties" = 1,
                                                  "Structure-derived TCR Weight" = 2
                                                  ),
                                        selected = 1),
                           fluidRow(
                               column(width = 3,
                                      numericInput("p1", 'P1', value = "3"),
                               ),
                               column(width = 3,
                                      numericInput("p2", 'P1', value = "0.5"),
                               ),
                               column(width = 3,
                                      numericInput("p3", 'P2', value = "0.5"),
                               ),
                           ),
                           fluidRow(
                             column(width = 3,
                                    numericInput("p4", 'P4', value = "4"),
                             ),
                             column(width = 3,
                                    numericInput("p5", 'P5', value = "2"),
                             ),
                             column(width = 3,
                                    numericInput("p6", 'P6', value = "0.5"),
                             )
                           ),
                           fluidRow(
                             column(width = 3,
                                    numericInput("p7", 'P7', value = "1"),
                             ),
                             column(width = 3,
                                    numericInput("p8", 'P8', value = "1"),
                             ),
                             column(width = 3,
                                    numericInput("p9", 'P9', value = "0.5"),
                             )
                           ),
                           actionButton("submit", "Submit", class = "btn-primary"),
                           width = 3
                       ),
                       mainPanel(
                           textOutput('input_text'),
                           hr(),
                           DT::dataTableOutput('table'),
                           # downloadButton('download',"Download the data"),
                           width = 8
                       )
                    )
                ),
                tabPanel("Analysis", includeHTML("analysis.html")),
                tabPanel("Help", includeHTML("help.html")),
                tabPanel("About", includeHTML("about.html"))
)

server <- function(input, output) {

    crossdomeParameters <- eventReactive(input$submit, {

      if(input$protocol == 1) {

        paste(
          input$peptide,
          input$allele,
          "Biochemical properities",
          sep = " - ")

        } else {

        paste(
          input$peptide,
          input$allele,
          paste0("Structural weights: ",
                paste0(c(input$p1, input$p2, input$p3, input$p4, input$p5, input$p6, input$p7, input$p8, input$p9), collapse = ", ")),
          sep = " - ")

      }

    })

    crossdomeInput <- eventReactive(input$submit, {

        if(TRUE) {

            position_weight <- list(
                `1` = rep(1, 9),
                `2` = c(
                        input$p1,
                        input$p2,
                        input$p3,
                        input$p4,
                        input$p5,
                        input$p6,
                        input$p7,
                        input$p8,
                        input$p9
                    )
            )

            position_weight <- position_weight[[input$protocol]]
            database <- crossdome::cross_background(off_targets = NULL, allele = input$allele)

            cross_result <- data.frame()

            withProgress(
                message = "Calculating...", value = 0, {
                    cross_result <- crossdome::cross_compose(
                        query = input$peptide,
                        background = database,
                        position_weight = position_weight
                    )
                }
            )
        }

        cross_result <- cross_result@result
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

    # To finish
    output$selected <- renderPrint({
        selected <- input$table_rows_selected
        if(selected) {
            cat(selected, sep = ",")
        } else {
            cat(1:50)
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

shinyApp(ui = ui, server = server)
