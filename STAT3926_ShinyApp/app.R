

library(shiny)
library(ggplot2)
library(scales)
library(dplyr)

# ── Core calculation functions ───────────────────────────────
calc_power <- function(N_total, delta, sigma, alpha = 0.05, k = 2) {
  z_a   <- qnorm(1 - alpha / 2)
  n_ctrl <- N_total / (1 + k)
  n_trt  <- k * n_ctrl
  se     <- sigma * sqrt(1/n_ctrl + 1/n_trt)
  ncp    <- delta / se
  pnorm(ncp - z_a) + pnorm(-ncp - z_a)
}

calc_enrolled_n <- function(sigma, delta, alpha = 0.05, power = 0.90,
                             k = 2, loss, withdraw, noncompliance) {
  z_a      <- qnorm(1 - alpha / 2)
  z_b      <- qnorm(power)
  n_ctrl   <- (z_a + z_b)^2 * sigma^2 * (1 + 1/k) / delta^2
  n_trt    <- k * n_ctrl
  base     <- ceiling(n_ctrl + n_trt)
  inflate  <- 1 / ((1 - loss) * (1 - withdraw) * (1 - noncompliance))
  ceiling(base * inflate)
}

calc_mdd <- function(N_total, sigma, alpha = 0.05, power = 0.90, k = 2) {
  z_a    <- qnorm(1 - alpha / 2)
  z_b    <- qnorm(power)
  n_ctrl <- N_total / (1 + k)
  n_trt  <- k * n_ctrl
  (z_a + z_b) * sigma * sqrt(1/n_ctrl + 1/n_trt)
}

# ── Colour palette ───────────────────────────────────────────
PAL <- c("#2166ac", "#4dac26", "#d6604d", "#b2182b",
         "#762a83", "#e07b39")

# ============================================================
# UI
# ============================================================
ui <- fluidPage(

  tags$head(tags$style(HTML("
    body        { font-family: Georgia, serif; background: #f9f9f7; color: #222; }
    h2          { color: #2166ac; font-size: 1.3em; margin-top: 0; }
    h4          { color: #444; font-size: 0.95em; margin-bottom: 4px; }
    .well       { background: #fff; border: 1px solid #ddd;
                  border-radius: 6px; padding: 16px; }
    .metric-box { background: #fff; border: 1px solid #ccc;
                  border-radius: 6px; padding: 14px 18px;
                  margin-bottom: 10px; text-align: center; }
    .metric-val { font-size: 2em; font-weight: bold; color: #2166ac; }
    .metric-lbl { font-size: 0.82em; color: #666; margin-top: 2px; }
    .warn       { color: #b2182b; font-weight: bold; }
    .ok         { color: #4dac26; font-weight: bold; }
    .note       { font-size: 0.82em; color: #666; margin-top: 6px; }
    .tab-header { font-size: 1em; font-weight: bold; }
    hr          { border-color: #ddd; }
    .shiny-input-container { margin-bottom: 10px; }
  "))),

  # Header
  div(style = "background:#2166ac; color:white; padding:14px 20px;
               margin-bottom:18px; border-radius:6px;",
    h2(style = "color:white; margin:0;",
       "Sample Size Sensitivity Tool"),
    div(style = "font-size:0.85em; margin-top:4px; opacity:0.9;",
        "Combination therapy vs placebo | Acute low back pain RCT |
         University of Sydney Statistical Consulting Service")
  ),

  fluidRow(

    # ── Left panel: inputs ──────────────────────────────────
    column(3,
      div(class = "well",
        h4("Study Parameters"),
        hr(),

        sliderInput("delta", "Mean difference to detect (NRS units)",
                    min = 0.25, max = 3.0, value = 1.0, step = 0.05),

        sliderInput("sigma", "Assumed standard deviation (SD)",
                    min = 0.5, max = 4.0, value = 1.0, step = 0.05),

        sliderInput("power_target", "Target power",
                    min = 0.70, max = 0.99, value = 0.90, step = 0.01),

        sliderInput("alpha", "Significance level (alpha)",
                    min = 0.01, max = 0.10, value = 0.05, step = 0.005),

        sliderInput("k_ratio", "Allocation ratio (treatment : control)",
                    min = 1.0, max = 4.0, value = 2.0, step = 0.5),

        hr(),
        h4("Attrition Assumptions"),

        sliderInput("loss", "Loss to follow-up",
                    min = 0, max = 0.40, value = 0.20, step = 0.01,
                    post = "%",
                    animate = FALSE),

        sliderInput("withdraw", "Withdrawal rate",
                    min = 0, max = 0.40, value = 0.15, step = 0.01,
                    post = "%"),

        sliderInput("noncompliance", "Non-compliance rate",
                    min = 0, max = 0.30, value = 0.10, step = 0.01,
                    post = "%"),

        hr(),
        h4("Reference line on plots"),
        numericInput("ref_n", "Proposed N (reference line)",
                     value = 324, min = 10, max = 2000, step = 1)
      )
    ),

    # ── Right panel: outputs ────────────────────────────────
    column(9,

      # Metric cards
      fluidRow(
        column(3, div(class = "metric-box",
          div(class = "metric-val", textOutput("card_n")),
          div(class = "metric-lbl", "Required enrolled N"),
          div(class = "note", "(post-attrition)")
        )),
        column(3, div(class = "metric-box",
          div(class = "metric-val", textOutput("card_power")),
          div(class = "metric-lbl", "Power at proposed N"),
          div(class = "note", uiOutput("card_power_flag"))
        )),
        column(3, div(class = "metric-box",
          div(class = "metric-val", textOutput("card_mdd")),
          div(class = "metric-lbl", "Min. detectable difference"),
          div(class = "note", "at proposed N (NRS units)")
        )),
        column(3, div(class = "metric-box",
          div(class = "metric-val", textOutput("card_d")),
          div(class = "metric-lbl", "Cohen's d"),
          div(class = "note", uiOutput("card_d_note"))
        ))
      ),

      br(),

      # Tabs for figures
      tabsetPanel(
        tabPanel("Power vs N",
          br(),
          plotOutput("plot_power", height = "380px"),
          div(class = "note", style = "padding:8px 4px;",
            "Each curve shows power across sample sizes for a different SD value.
             The vertical line marks the proposed N. The horizontal dashed lines
             mark 80% and 90% power targets.")
        ),
        tabPanel("Required N vs SD",
          br(),
          plotOutput("plot_n_sd", height = "380px"),
          div(class = "note", style = "padding:8px 4px;",
            "Shows how the required enrolled sample size changes as the assumed SD
             increases. The dashed red line is the proposed N; the grey line is
             the feasibility limit of 1000.")
        ),
        tabPanel("Minimum Detectable Difference",
          br(),
          plotOutput("plot_mdd", height = "380px"),
          div(class = "note", style = "padding:8px 4px;",
            "Shows the smallest mean difference detectable at 90% power
             as a function of enrolled N. The horizontal dashed line marks
             the target effect of 1.0 NRS unit.")
        ),
        tabPanel("Scenario Table",
          br(),
          tableOutput("scenario_table"),
          div(class = "note", style = "padding:8px 4px;",
            "Base N = completers before attrition adjustment.
             Alloc-adj N = after 2:1 efficiency correction (x 9/8).
             Enrolled N = after multiplicative attrition inflation.
             Power and MDD are evaluated at the proposed N.")
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # Fix slider display as percentages
  observe({
    updateSliderInput(session, "loss",
                      value = input$loss,
                      label = paste0("Loss to follow-up  (",
                                     round(input$loss * 100), "%)"))
    updateSliderInput(session, "withdraw",
                      value = input$withdraw,
                      label = paste0("Withdrawal rate  (",
                                     round(input$withdraw * 100), "%)"))
    updateSliderInput(session, "noncompliance",
                      value = input$noncompliance,
                      label = paste0("Non-compliance rate  (",
                                     round(input$noncompliance * 100), "%)"))
  })

  # Reactive: required N under current inputs
  req_n <- reactive({
    calc_enrolled_n(input$sigma, input$delta, input$alpha,
                    input$power_target, input$k_ratio,
                    input$loss, input$withdraw, input$noncompliance)
  })

  pwr_at_ref <- reactive({
    calc_power(input$ref_n, input$delta, input$sigma,
               input$alpha, input$k_ratio)
  })

  mdd_at_ref <- reactive({
    calc_mdd(input$ref_n, input$sigma, input$alpha,
             input$power_target, input$k_ratio)
  })

  cohens_d <- reactive({ input$delta / input$sigma })

  # ── Metric cards ─────────────────────────────────────────
  output$card_n <- renderText({
    n <- req_n()
    if (n > 5000) ">5000" else formatC(n, format = "d", big.mark = ",")
  })

  output$card_power <- renderText({
    paste0(round(pwr_at_ref() * 100, 1), "%")
  })

  output$card_power_flag <- renderUI({
    p <- pwr_at_ref()
    tgt <- input$power_target
    if (p < tgt - 0.001)
      span(class = "warn",
           paste0("Below target (", round(tgt*100), "%)"))
    else
      span(class = "ok",
           paste0("Meets target (", round(tgt*100), "%)"))
  })

  output$card_mdd <- renderText({
    round(mdd_at_ref(), 2)
  })

  output$card_d <- renderText({
    round(cohens_d(), 2)
  })

  output$card_d_note <- renderUI({
    d <- cohens_d()
    label <- if (d < 0.3) "small effect"
             else if (d < 0.5) "small-medium"
             else if (d < 0.8) "medium effect"
             else "large effect"
    span(style = "color:#444;", label)
  })

  # ── Plot 1: Power vs N ────────────────────────────────────
  output$plot_power <- renderPlot({
    N_seq  <- seq(10, min(req_n() * 2.5, 1500), length.out = 400)
    sd_scenarios <- c(1.0, 2.0, 2.1, 2.6)
    sd_labels    <- c("SD = 1.0 (client assumption)",
                      "SD = 2.0 (optimistic estimate)",
                      "SD = 2.1 (literature estimate)",
                      "SD = 2.6 (conservative estimate)")

    df <- do.call(rbind, lapply(seq_along(sd_scenarios), function(i) {
      data.frame(
        N     = N_seq,
        Power = calc_power(N_seq, input$delta, sd_scenarios[i],
                           input$alpha, input$k_ratio),
        SD    = sd_labels[i]
      )
    }))
    df$SD <- factor(df$SD, levels = sd_labels)

    ggplot(df, aes(x = N, y = Power, colour = SD)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = input$power_target, linetype = "dashed",
                 colour = "#b2182b", linewidth = 0.8) +
      geom_hline(yintercept = 0.80, linetype = "dashed",
                 colour = "#f4a582", linewidth = 0.7) +
      geom_vline(xintercept = input$ref_n, linetype = "dotted",
                 colour = "#e07b39", linewidth = 1.0) +
      annotate("text", x = input$ref_n + 12, y = 0.12,
               label = paste0("N = ", input$ref_n),
               colour = "#e07b39", hjust = 0, size = 3.5,
               family = "serif") +
      scale_colour_manual(values = PAL) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         limits = c(0, 1.02)) +
      labs(x = "Total Enrolled N", y = "Power (1 - β)",
           title = "Power vs Sample Size Across SD Assumptions",
           subtitle = paste0("δ = ", input$delta,
                             ",  α = ", input$alpha,
                             ",  allocation ratio = ", input$k_ratio, ":1"),
           colour = NULL) +
      theme_minimal(base_family = "serif", base_size = 12) +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(colour = "#2166ac", face = "bold"))
  })

  # ── Plot 2: Required N vs SD ──────────────────────────────
  output$plot_n_sd <- renderPlot({
    sd_seq <- seq(0.5, 4.0, by = 0.05)
    n_seq  <- sapply(sd_seq, function(s)
      calc_enrolled_n(s, input$delta, input$alpha, input$power_target,
                      input$k_ratio, input$loss, input$withdraw,
                      input$noncompliance))

    df <- data.frame(SD = sd_seq, N = pmin(n_seq, 5000))

    ggplot(df, aes(x = SD, y = N)) +
      geom_line(colour = "#2166ac", linewidth = 1.2) +
      geom_hline(yintercept = input$ref_n, linetype = "dashed",
                 colour = "#e07b39", linewidth = 0.9) +
      geom_hline(yintercept = 1000, linetype = "dotted",
                 colour = "#888", linewidth = 0.8) +
      geom_vline(xintercept = input$sigma, linetype = "dotted",
                 colour = "#d6604d", linewidth = 0.9) +
      annotate("text", x = input$sigma + 0.06, y = 200,
               label = paste0("Current SD = ", input$sigma),
               colour = "#d6604d", hjust = 0, size = 3.5,
               family = "serif") +
      annotate("text", x = 0.6, y = input$ref_n + 60,
               label = paste0("Proposed N = ", input$ref_n),
               colour = "#e07b39", hjust = 0, size = 3.5,
               family = "serif") +
      annotate("text", x = 0.6, y = 1060,
               label = "Feasibility limit ~1000",
               colour = "#888", hjust = 0, size = 3.2,
               family = "serif") +
      scale_y_continuous(labels = comma, limits = c(0, 2200)) +
      labs(x = "Assumed SD (σ)", y = "Total Enrolled N (post-attrition)",
           title = "Required Enrolled N for Different SD Assumptions",
           subtitle = paste0("δ = ", input$delta,
                             ",  power = ", round(input$power_target*100), "%",
                             ",  α = ", input$alpha,
                             ",  ratio = ", input$k_ratio, ":1")) +
      theme_minimal(base_family = "serif", base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(colour = "#2166ac", face = "bold"))
  })

  # ── Plot 3: Minimum Detectable Difference vs N ───────────
  output$plot_mdd <- renderPlot({
    N_seq        <- seq(10, min(req_n() * 2.5, 1500), length.out = 400)
    sd_scenarios <- c(1.0, 2.0, 2.1, 2.6)
    sd_labels    <- c("SD = 1.0", "SD = 2.0", "SD = 2.1", "SD = 2.6")

    df <- do.call(rbind, lapply(seq_along(sd_scenarios), function(i) {
      data.frame(
        N   = N_seq,
        MDD = calc_mdd(N_seq, sd_scenarios[i],
                       input$alpha, input$power_target, input$k_ratio),
        SD  = sd_labels[i]
      )
    }))
    df$SD <- factor(df$SD, levels = sd_labels)

    ggplot(df, aes(x = N, y = MDD, colour = SD)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = input$delta, linetype = "dashed",
                 colour = "black", linewidth = 0.8) +
      geom_vline(xintercept = input$ref_n, linetype = "dotted",
                 colour = "#e07b39", linewidth = 1.0) +
      annotate("text", x = input$ref_n + 12, y = max(df$MDD) * 0.9,
               label = paste0("N = ", input$ref_n),
               colour = "#e07b39", hjust = 0, size = 3.5,
               family = "serif") +
      annotate("text", x = 20, y = input$delta + 0.08,
               label = paste0("Target δ = ", input$delta),
               hjust = 0, size = 3.5, family = "serif") +
      scale_colour_manual(values = PAL) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(x = "Total Enrolled N",
           y = "Minimum Detectable Difference (NRS units)",
           title = "Minimum Detectable Difference vs Sample Size",
           subtitle = paste0("Power = ", round(input$power_target*100), "%",
                             ",  α = ", input$alpha,
                             ",  ratio = ", input$k_ratio, ":1"),
           colour = NULL) +
      theme_minimal(base_family = "serif", base_size = 12) +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(colour = "#2166ac", face = "bold"))
  })

  # ── Scenario table ────────────────────────────────────────
  output$scenario_table <- renderTable({
    sds <- c(1.0, 1.5, 2.0, 2.1, 2.6)
    z_a <- qnorm(1 - input$alpha / 2)
    z_b <- qnorm(input$power_target)
    k   <- input$k_ratio
    inf <- 1 / ((1 - input$loss) * (1 - input$withdraw) *
                (1 - input$noncompliance))

    rows <- lapply(sds, function(s) {
      d_val    <- input$delta / s
      n_ctrl   <- (z_a + z_b)^2 * s^2 * (1 + 1/k) / input$delta^2
      n_trt    <- k * n_ctrl
      base     <- ceiling(n_ctrl + n_trt)
      alloc    <- ceiling(base * ((k + 1)^2 / (4 * k)))
      enrolled <- ceiling(base * inf)
      pwr      <- calc_power(input$ref_n, input$delta, s,
                             input$alpha, k)
      mdd      <- calc_mdd(input$ref_n, s, input$alpha,
                           input$power_target, k)
      data.frame(
        "SD (σ)"             = s,
        "Cohen's d"          = round(d_val, 2),
        "Base N"             = base,
        "Alloc-adj N"        = alloc,
        "Enrolled N"         = enrolled,
        "Power at ref N"     = paste0(round(pwr * 100, 1), "%"),
        "MDD at ref N"       = round(mdd, 2),
        check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, bordered = TRUE,
     align = "c", digits = 0)
}

# ============================================================
shinyApp(ui, server)
