# app.R
# Shiny app: PR curve with clickable plotly points that report "how many investigated"

library(shiny)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  titlePanel("Precision vs Recall — Click to see how many investigated"),
  fluidRow(
    column(
      width = 8,
      plotlyOutput("prPlot", height = 600)
    ),
    column(
      width = 4,
      h4("Clicked point"),
      verbatimTextOutput("clickInfo"),
      tags$hr(),
      h5("Notes"),
      tags$ul(
        tags$li("Click a point on the curve to see the index (number investigated)."),
        tags$li("Points are clickable; the line is for visual continuity and has hover disabled.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  # ---- Generate data (vectorized; same distributional idea as your script) ----
  ps_df <- reactive({
    set.seed(1)  # for reproducibility; remove or expose as input if you like
    
    N       <- 1e4
    propH1  <- 0.05
    nH1     <- round(propH1 * N)
    nH0     <- N - nH1
    
    # H1 p-values (scrunched near 0)
    tmp     <- rexp(2 * nH1)
    psH1    <- (tmp / (2 * max(tmp)))
    psH1    <- sort(psH1)[1:nH1]^1.2
    
    # H0 p-values Uniform(0,1)
    psH0    <- runif(nH0)
    
    ps      <- c(psH1, psH0)
    H1_ind  <- c(rep(1L, nH1), rep(0L, nH0))
    
    df <- data.frame(pVal = ps, H1.ind = H1_ind)
    df <- df[order(df$pVal), , drop = FALSE]           # investigate from smallest p upwards
    df$k <- seq_len(nrow(df))                           # row index = how many investigated
    
    # Fast cumulative counts
    tp <- cumsum(df$H1.ind)                             # true positives up to k
    fp <- df$k - tp                                     # false positives up to k
    tot_H1 <- sum(df$H1.ind)
    fn <- tot_H1 - tp                                   # false negatives remaining
    
    precision <- tp / (tp + fp)
    recall    <- tp / (tp + fn)
    
    # protect against 0/0 cases (only at k=0 which we don't have)
    precision[is.nan(precision)] <- NA_real_
    
    transform(df,
              tp = tp, fp = fp, fn = fn,
              precision = precision, recall = recall)
  })
  
  # ---- Plot (ggplot -> plotly); points carry a "key" so clicks return the index ----
  output$prPlot <- renderPlotly({
    df <- ps_df()
    
    g <- ggplot(df, aes(x = recall, y = precision)) +
      geom_path(linewidth = 1, alpha = 0.7) +
      geom_point(
        aes(
          customdata = k,   # for click events
          text = paste0(
            "Investigated (k): ", k,
            "<br>p-value threshold: ", signif(pVal, 4),
            "<br>TP: ", tp, " | FP: ", fp, " | FN: ", fn,
            "<br>Precision: ", signif(precision, 4),
            "<br>Recall: ", signif(recall, 4),
            "<br>H1 indicator: ", H1.ind
          )
        ),
        size = 1.0, alpha = 0.5
      ) +
      scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
      coord_equal() +
      labs(
        title = "Screen Simulation Summary: Precision vs Recall",
        x = "Recall (TP / (TP + FN))",
        y = "Precision (TP / (TP + FP))"
      ) +
      theme_minimal(base_size = 12)
    
    p <- ggplotly(g, tooltip = "text")
    
    # Disable hover for line trace
    if (length(p$x$data) >= 1) {
      path_ids <- which(vapply(p$x$data, function(tr) tr$mode == "lines", logical(1)))
      if (length(path_ids)) {
        p <- style(p, hoverinfo = "skip", traces = path_ids)
      }
    }
    
    p
  })
  
  
  # ---- Click info ----
  output$clickInfo <- renderPrint({
    ed <- event_data("plotly_click")
    if (is.null(ed) || is.null(ed$customdata)) {
      cat("Click a point on the chart.\n")
      return()
    }
    
    k <- as.integer(ed$customdata)  # <- use customdata here
    df <- ps_df()
    row <- df[k, ]
    
    out <- list(
      investigated_k = row$k,
      p_value_threshold = signif(row$pVal, 6),
      TP = row$tp,
      FP = row$fp,
      FN = row$fn,
      precision = signif(row$precision, 6),
      recall = signif(row$recall, 6),
      H1_indicator = row$H1.ind
    )
    print(out)
  })
  
}

shinyApp(ui, server)
