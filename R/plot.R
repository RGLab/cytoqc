
#' @export
plot.cqc_match_result_and_solution <- function(res){
  # Decide which one we call
}

#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_brewer scale_fill_manual 
#' labs coord_flip theme element_text
#' @importFrom tidyr pivot_longer pivot_wider crossing
#' @importFrom tibble as_tibble
#' @import dplyr dplyr
#' @importFrom shiny shinyApp fluidPage fluidRow column 
#' wellPanel selectInput plotOutput renderPlot renderText 
#' nearPoints verbatimTextOutput hoverOpts runGadget
shiny_plot <- function(res){
  trimmed <- format(res)
  tbl <- trimmed[trimmed$Ref != "",-1] %>%
    mutate(Ref = factor(Ref, levels=rev(levels(Ref)))) %>%
    as_tibble() %>% 
    pivot_longer(-Ref, names_to = "Group", values_to = "Match") %>%
    mutate(Match_Case = factor(case_when(
      Match == "\u2713" ~ "Exact Match",
      is.na(Match)      ~ "Missing",
      TRUE              ~ "Approximate Match"
    ))) %>%
    mutate(Fill_Text = factor(case_when(
      Match == "\u2713" | is.na(Match) ~ "",
      TRUE                             ~ Match
    )
    ))
  
  tileplot <- ggplot(tbl, aes(Group, Ref, fill=Match_Case)) + 
            geom_tile(color="black") + 
            scale_fill_manual(
              values=c("#d1495b", "#edae49", "#66a172"), #red, yellow, green
              limits=c("Missing", "Approximate Match", "Exact Match"),
              name="Match Status") +
            labs(x="Reference", y="Group") +
            theme(axis.title=element_text(size=14, face="bold"),
                  axis.text=element_text(size=12, face="bold"),
                  legend.text=element_text(size=12, face="bold"),
                  legend.title=element_text(size=14, face="bold"))
  
  unmatched_info <- trimmed[trimmed$Ref=="", -c(1,2)]
  unmatched_info <- sapply(colnames(unmatched_info), function(cn){
    paste0("Group ", cn, ": ", unmatched_info[,cn])
    })
  unmatched_info <- paste("Unmatched Values", paste(unmatched_info, collapse="\n"), sep="\n")
  
  ui = fluidPage(
    fluidRow(
        plotOutput("tileplot", hover = hoverOpts("tileplot_hover", delay=100, delayType = "throttle"))
    ),
    fluidRow(
      verbatimTextOutput("hover_info")
    ),
    fluidRow(
      verbatimTextOutput("unmatched_info")
    )
  )
  server = function(input, output){
    output$tileplot <- renderPlot(tileplot)
    output$unmatched_info <- renderText(unmatched_info)
    observeEvent(input$tileplot_hover, {
      hoverpoint <- nearPoints(tbl,
                               input$tileplot_hover,
                               threshold = 100, maxpoints = 1
      )
      output$hover_info <- renderText(
        paste0(
          case_when(
            hoverpoint$Match_Case == "Exact Match" ~ paste0("Exact match to reference: ", hoverpoint$Ref),
            hoverpoint$Match_Case == "Missing"     ~ paste0("Missing marker: ", hoverpoint$Ref),
            TRUE                                   ~ paste0("Approximate match:  ", hoverpoint$Match, " => ", hoverpoint$Ref)
          )
        )
      )
    }, ignoreNULL = FALSE)
  }
  runGadget(shinyApp(ui,server))
}

#' @importFrom dplyr mutate_all
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling cell_spec
kable_plot <- function(res){
  trimmed <- format(res)
  unmatched <- trimmed[trimmed$Ref == "",-1]
  unmatched$Ref <- "Unmatched"
  trimmed <- trimmed[trimmed$Ref != "",-1]
  trimmed <- trimmed %>%
    mutate_at(-1,function(col){
      cell_spec(col, "html", color = "black", background = case_when(
        col == "\u2713" ~ "#66a172",
        is.na(col)      ~ "#d1495b", #red
        TRUE            ~ "#edae49" #yellow
      ))
    })
  rbind(trimmed, unmatched) %>%
    kable(format="html", escape=F) %>%
    kable_styling()
}

