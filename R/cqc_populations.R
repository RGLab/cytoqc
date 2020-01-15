#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom shiny shinyApp fluidPage fluidRow column wellPanel selectInput plotOutput renderPlot nearPoints
#' @importFrom ggcyto autoplot
#' @export
cqc_gs_check_populations <- function(gs, nodes = NULL, plot = FALSE) {
  counts <- gs_pop_get_stats(gs, nodes = nodes)
  counts_split <- counts %>%
    mutate(idx = as.numeric(row.names(counts))) %>%
    group_split(pop)
  names(counts_split) <- sapply(counts_split, function(subcounts) subcounts[1, "pop"])

  if (plot) {
    plots <- lapply(counts_split, qoutlier_plot)
    shinyApp(
      ui = fluidPage(
        fluidRow(
          column(
            12,
            wellPanel(
              selectInput(
                inputId = "pop",
                label = "Choose a gated subpopulation:",
                choices = gs_get_pop_paths(gs)
              )
            )
          )
        ),
        fluidRow(
          column(
            6,
            "Click a sample point to view gate.\n
             Red points are beyond 1.5*IQR:"
          )
        ),
        fluidRow(
          column(
            6,
            wellPanel(
              plotOutput("boxplot", click = "boxplot_click")
            )
          ),
          column(
            6,
            wellPanel(
              # For diagnostics
              # verbatimTextOutput("info")
              plotOutput("gate")
            )
          )
        )
      ),
      server = function(input, output) {
        output$boxplot <- renderPlot({
          plots[[input$pop]]
        })
        # For diagnostics
        # output$info <- renderPrint({
        #   nearPoints(counts_split[[input$pop]], input$boxplot_click, threshold = 10, maxpoints = 1, addDist = TRUE)
        # })
        observeEvent(input$boxplot_click, {
          clickedpoint <- nearPoints(counts_split[[input$pop]],
            input$boxplot_click,
            threshold = 10, maxpoints = 1, addDist = TRUE
          )
          output$gate <- renderPlot({
            autoplot(gs[[clickedpoint$sample]], input$pop)
          })
        })
      }
    )
  } else {
    wide <- pivot_wider(counts, id_cols = sample, names_from = pop, values_from = count)
    outlier_idx <- as.data.frame(lapply(wide[-1], function(col) {
      qoutlier(col, plot = FALSE)
    }))
    colnames(outlier_idx) <- colnames(wide)[-1]
    rownames(outlier_idx) <- wide[[1]]
    outlier_idx
  }
}
