# ─────────────────────────────────────────────────────────────────────────────
#  IES 2022/23 — Household Income Explorer
#  Mid-Century Modern · Browns & Warm Tones
#  Stats SA · Fact_IES2023_Households.csv
# ─────────────────────────────────────────────────────────────────────────────

library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(scales)

# ── Label look-ups ─────────────────────────────────────────────────────────────
sex_lbl     <- c("1"="Male",          "2"="Female")
pop_lbl     <- c("1"="Black African", "2"="Coloured","3"="Indian/Asian","4"="White")
status_lbl  <- c("1"="Wealthy",       "2"="Very comfortable",
                 "3"="Reasonably comfortable",
                 "4"="Just getting along","5"="Poor","6"="Very poor")
dwell_lbl   <- c("1"="Formal house/brick","2"="Traditional hut",
                 "3"="Flat/apartment",    "4"="Cluster house",
                 "5"="Town house",        "6"="Semi-detached",
                 "7"="House/flat in backyard","8"="Informal shack (backyard)",
                 "9"="Informal shack (other)","10"="Room/granny flat",
                 "11"="Caravan/tent","12"="Other")
elec_lbl    <- c("1"="Yes","2"="No")

# MCM palette
PAL <- c("#A0522D","#D4783E","#C4956A","#8B7D6B",
         "#6B4226","#3D2B1F","#EDD9C0","#D4A76A","#B8865C","#7A5C3E")

#estimates the most common income level for a population group by smoothing the data into a 
#continuous curve and identifying the point where that curve reaches its highest peak.
density_mode <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  d <- density(x)
  d$x[which.max(d$y)]
}

# ── ggplot theme ───────────────────────────────────────────────────────────────
theme_mcm <- function() {
  theme_minimal(base_family = "serif") +
    theme(
      plot.background   = element_rect(fill="#F5E6D3", colour=NA),
      panel.background  = element_rect(fill="#EDD9C0", colour=NA),
      panel.grid.major  = element_line(colour="#C4956A50"),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(colour="#3D2B1F", size=10),
      axis.title        = element_text(colour="#3D2B1F", size=11, face="bold"),
      plot.title        = element_text(colour="#3D2B1F", size=13, face="bold"),
      plot.subtitle     = element_text(colour="#6B4226", size=9),
      legend.background = element_rect(fill="#F5E6D3", colour=NA),
      legend.text       = element_text(colour="#3D2B1F", size=9),
      legend.title      = element_text(colour="#3D2B1F", face="bold", size=9),
      strip.background  = element_rect(fill="#6B4226", colour=NA),
      strip.text        = element_text(colour="#FAF3EA", face="bold")
    )
}

# ── CSS ────────────────────────────────────────────────────────────────────────
mcm_css <- "
@import url('https://fonts.googleapis.com/css2?family=Playfair+Display:wght@400;700&family=Lato:wght@300;400;700&display=swap');

* { box-sizing: border-box; }

body {
  background-color: #F5E6D3;
  font-family: 'Lato', sans-serif;
  color: #3D2B1F;
  margin: 0; padding: 0;
}

/* ── Header ── */
.app-header {
  background: linear-gradient(135deg, #3D2B1F 0%, #6B4226 60%, #A0522D 100%);
  color: #FAF3EA;
  padding: 22px 32px 18px;
  border-bottom: 4px solid #D4783E;
  margin-bottom: 0;
}
.app-header h1 {
  font-family: 'Playfair Display', Georgia, serif;
  font-size: 26px;
  margin: 0 0 4px 0;
  letter-spacing: 0.5px;
}
.app-header p {
  font-size: 11px;
  margin: 0;
  color: #C4956A;
  letter-spacing: 1.5px;
  text-transform: uppercase;
}

/* ── Sidebar ── */
.sidebar-wrap {
  background-color: #EDD9C0;
  border-right: 3px solid #C4956A;
  padding: 20px 16px 40px;
  min-height: calc(100vh - 80px);
}
.sidebar-wrap h4 {
  font-family: 'Playfair Display', serif;
  color: #3D2B1F;
  border-bottom: 2px solid #A0522D;
  padding-bottom: 5px;
  margin: 20px 0 10px;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 1.2px;
}
.sidebar-wrap h4:first-child { margin-top: 0; }

/* ── Path bar ── */
.path-bar {
  background-color: #FAF3EA;
  border: 2px dashed #C4956A;
  padding: 12px;
  border-radius: 5px;
  margin-bottom: 18px;
}
.path-bar label { font-size: 11px; font-weight: 700; color: #3D2B1F; }

/* ── Inputs ── */
.selectize-input, .selectize-dropdown {
  background-color: #FAF3EA !important;
  border: 1px solid #C4956A !important;
  color: #3D2B1F !important;
}
.selectize-input.focus {
  border-color: #A0522D !important;
  box-shadow: 0 0 0 2px #D4783E30 !important;
}
.form-control {
  background-color: #FAF3EA;
  border: 1px solid #C4956A;
  color: #3D2B1F;
  font-size: 13px;
}
.form-control:focus {
  border-color: #A0522D;
  box-shadow: 0 0 0 2px #D4783E30;
}
.irs--shiny .irs-bar { background:#A0522D; border-top:1px solid #A0522D; border-bottom:1px solid #A0522D; }
.irs--shiny .irs-handle { background:#D4783E; border:2px solid #A0522D; }
.irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background:#6B4226; }

/* ── Buttons ── */
.btn-mcm {
  background-color: #A0522D;
  color: #FAF3EA;
  border: none;
  font-family: 'Lato', sans-serif;
  font-weight: 700;
  letter-spacing: 0.5px;
  padding: 8px 16px;
  border-radius: 4px;
  width: 100%;
  margin-top: 8px;
  cursor: pointer;
}
.btn-mcm:hover { background-color: #6B4226; color: #FAF3EA; }
.btn-reset {
  background-color: #8B7D6B !important;
}
.btn-reset:hover { background-color: #6B4226 !important; }

/* ── Tabs ── */
.nav-tabs { border-bottom: 3px solid #A0522D; }
.nav-tabs > li > a {
  font-family: 'Lato', sans-serif;
  font-weight: 700;
  color: #6B4226;
  letter-spacing: 0.4px;
  border: 1px solid transparent;
  border-radius: 4px 4px 0 0;
}
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  background-color: #A0522D;
  color: #FAF3EA;
  border-color: #A0522D;
}
.nav-tabs > li > a:hover {
  background-color: #EDD9C0;
  border-color: #C4956A;
  color: #3D2B1F;
}
.tab-content {
  background-color: #FAF3EA;
  border: 1px solid #C4956A;
  border-top: none;
  padding: 24px;
  border-radius: 0 0 4px 4px;
}

/* ── Stat cards ── */
.stat-card {
  background: linear-gradient(135deg, #6B4226, #A0522D);
  border-radius: 6px;
  padding: 15px 18px;
  color: #FAF3EA;
  margin-bottom: 14px;
  box-shadow: 2px 3px 10px #3D2B1F25;
}
.stat-card-dark { background: linear-gradient(135deg, #3D2B1F, #6B4226); }
.stat-card-tan  { background: linear-gradient(135deg, #C4956A, #D4783E); }
.stat-card-tan .stat-label { color: #3D2B1F80; }
.stat-card-tan .stat-value { color: #3D2B1F; }
.stat-label {
  font-size: 9px;
  text-transform: uppercase;
  letter-spacing: 1.8px;
  color: #C4956A;
  margin-bottom: 5px;
}
.stat-value {
  font-family: 'Playfair Display', serif;
  font-size: 22px;
  font-weight: 700;
  line-height: 1.2;
}

/* ── Section headings ── */
.section-head {
  font-family: 'Playfair Display', serif;
  color: #3D2B1F;
  font-size: 16px;
  border-left: 4px solid #D4783E;
  padding-left: 10px;
  margin: 0 0 16px 0;
}

/* ── Badge ── */
.n-badge {
  display: inline-block;
  background-color: #D4783E;
  color: #FAF3EA;
  font-size: 11px;
  font-weight: 700;
  padding: 2px 11px;
  border-radius: 12px;
  margin-left: 8px;
  vertical-align: middle;
}
.weighted-note {
  margin-left: 12px;
  font-size: 11px;
  color: #6B4226;
  vertical-align: middle;
}

/* ── Data table ── */
.dataTables_wrapper { font-family: 'Lato', sans-serif; font-size: 12px; }
table.dataTable thead th {
  background-color: #6B4226 !important;
  color: #FAF3EA !important;
  border-bottom: 2px solid #A0522D !important;
}
table.dataTable tbody tr { background-color: #FAF3EA !important; }
table.dataTable tbody tr:nth-child(even) { background-color: #F5E6D3 !important; }
table.dataTable tbody tr:hover td { background-color: #EDD9C0 !important; }
.dataTables_filter input { border: 1px solid #C4956A; background: #FAF3EA; color: #3D2B1F; }
.dataTables_length select { border: 1px solid #C4956A; background: #FAF3EA; }

/* ── Summary table ── */
table.summary-tbl {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}
table.summary-tbl th {
  background-color: #6B4226;
  color: #FAF3EA;
  padding: 8px 12px;
  text-align: left;
  letter-spacing: 0.5px;
  font-size: 11px;
  text-transform: uppercase;
}
table.summary-tbl td {
  padding: 7px 12px;
  border-bottom: 1px solid #C4956A30;
  color: #3D2B1F;
}
table.summary-tbl tr:nth-child(even) td { background-color: #F5E6D3; }
table.summary-tbl tr:hover td { background-color: #EDD9C0; }

/* ── Main content ── */
.main-content { padding: 18px 22px; }
.record-bar { margin-bottom: 16px; padding: 8px 0; border-bottom: 1px solid #C4956A40; }
"

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  title = "IES 2022/23 Income Explorer",
  tags$head(tags$style(HTML(mcm_css))),

  # ── App header ──
  div(class = "app-header",
    tags$h1("\U0001F3E0  IES 2022/23 — Household Income Explorer"),
    tags$p("Income and Expenditure Survey · Statistics South Africa · Microdata Explorer")
  ),

  fluidRow(

    # ── Sidebar ──────────────────────────────────────────────────────────────
    column(3,
      div(class = "sidebar-wrap",

        # File path loader
        div(class = "path-bar",
          tags$label("📂  CSV File Path"),
          textInput("csv_path",
                    label    = NULL,
                    value    = "C:/Users/arabe/OneDrive/Documents/Research Project/Fact_IES2023_Households.csv",
                    placeholder = "Full path to Fact_IES2023_Households.csv"),
          actionButton("load_btn", "⬆ Load Data", class = "btn-mcm")
        ),

        uiOutput("filter_ui")
      )
    ),

    # ── Main panel ────────────────────────────────────────────────────────────
    column(9,
      div(class = "main-content",

        uiOutput("record_bar"),

        tabsetPanel(id = "tabs",

          # ── Tab 1: Summary stats ──
          tabPanel("📊  Summary Statistics",
            br(),
            p(class = "section-head", "Key Income & Expenditure Metrics"),

            fluidRow(
              column(4, uiOutput("card1")),
              column(4, uiOutput("card2")),
              column(4, uiOutput("card3"))
            ),
            fluidRow(
              column(4, uiOutput("card4")),
              column(4, uiOutput("card5")),
              column(4, uiOutput("card6"))
            ),

            tags$hr(style = "border-color:#C4956A; margin: 8px 0 22px;"),

            fluidRow(
              column(6,
                p(class = "section-head", "Income by Decile"),
                tableOutput("tbl_decile")
              ),
              column(6,
                p(class = "section-head", "Income by Population Group"),
                tableOutput("tbl_pop")
              )
            )
          ),

          # ── Tab 2: Data table ──
          tabPanel("📋  Data Table",
            br(),
            DTOutput("main_dt")
          ),

          # ── Tab 3: Visualisations ──
          tabPanel("📈  Visualisations",
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Income Distribution"),
                plotOutput("plt_hist", height = "270px")
              ),
              column(6,
                p(class = "section-head", "Mode Income by Decile"),
                plotOutput("plt_decile", height = "270px")
              )
            ),
            br(),
            fluidRow(
              column(6,
                p(class = "section-head", "Income by Population Group"),
                plotOutput("plt_pop_box", height = "290px")
              ),
              column(6,
                p(class = "section-head", "Income per Capita by Household Size"),
                plotOutput("plt_hsize", height = "290px")
              )
            ),
            br(),
            fluidRow(
              column(12,
                p(class = "section-head", "Mode Income by Self-Reported Status"),
                plotOutput("plt_status", height = "250px")
              )
            )
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  raw  <- reactiveVal(NULL)

  # ── Helper: load CSV ────────────────────────────────────────────────────────
  load_csv <- function(path) {
    if (!file.exists(path)) {
      showNotification(paste0("File not found:\n", path),
                       type = "error", duration = 7)
      return(NULL)
    }
    withProgress(message = "Loading data…", value = 0.3, {
      df <- tryCatch(
        read.csv(path, stringsAsFactors = FALSE),
        error = function(e) {
          showNotification(paste("Read error:", e$message), type="error")
          NULL
        }
      )
      if (!is.null(df)) {
        names(df) <- toupper(names(df))
        setProgress(1)
        showNotification(
          paste0("✓  Loaded ", format(nrow(df), big.mark=","), " households"),
          type = "message", duration = 4
        )
      }
      df
    })
  }

  # Auto-load on startup if default path exists
  observe({
    if (is.null(isolate(raw()))) {
      p <- isolate(input$csv_path)
      if (nchar(p) > 0 && file.exists(p)) raw(load_csv(p))
    }
  })

  # Manual load button
  observeEvent(input$load_btn, {
    raw(load_csv(trimws(input$csv_path)))
  })

  # ── Dynamic filter UI ───────────────────────────────────────────────────────
  output$filter_ui <- renderUI({
    req(raw())
    df <- raw()

    age_max  <- if ("HEAD_AGE"  %in% names(df)) max(df$HEAD_AGE,  na.rm=TRUE) else 104
    size_max <- if ("HSIZE"     %in% names(df)) min(max(df$HSIZE, na.rm=TRUE), 20) else 20

    tagList(
      tags$h4("Income"),
      sliderInput("f_decile",   "Income Decile",   1, 10,       c(1,10),       step=1, ticks=FALSE),
      selectInput("f_quintile", "Income Quintile",
                  c("All"="all","1 — Lowest"="1","2"="2","3"="3",
                    "4"="4","5 — Highest"="5"), selected="all"),

      tags$h4("Head of Household"),
      selectInput("f_sex", "Sex",
                  c("All"="all","Male"="1","Female"="2"), selected="all"),
      selectInput("f_pop", "Population Group",
                  c("All"="all","Black African"="1","Coloured"="2",
                    "Indian/Asian"="3","White"="4"), selected="all"),
      sliderInput("f_age", "Age", 15, age_max, c(15, age_max),
                  step=1, ticks=FALSE),

      tags$h4("Household"),
      sliderInput("f_hsize", "Size", 1, size_max, c(1, size_max),
                  step=1, ticks=FALSE),
      selectInput("f_dwell", "Dwelling Type",
                  c("All"="all","Formal house/brick"="1","Traditional hut"="2",
                    "Flat/apartment"="3","Cluster house"="4",
                    "Town house"="5","Semi-detached"="6",
                    "House/flat in backyard"="7",
                    "Informal shack (backyard)"="8",
                    "Informal shack (other)"="9",
                    "Room/granny flat"="10","Caravan/tent"="11",
                    "Other"="12"), selected="all"),
      selectInput("f_status", "Present Status",
                  c("All"="all","Wealthy"="1","Very comfortable"="2",
                    "Reasonably comfortable"="3","Just getting along"="4",
                    "Poor"="5","Very poor"="6"), selected="all"),
      selectInput("f_elec", "Electricity Access",
                  c("All"="all","Yes"="1","No"="2"), selected="all"),

      actionButton("reset_btn", "\u21ba Reset All Filters",
                   class = "btn-mcm btn-reset")
    )
  })

  # Reset
  observeEvent(input$reset_btn, {
    df <- raw(); req(df)
    age_max  <- if ("HEAD_AGE" %in% names(df)) max(df$HEAD_AGE,  na.rm=TRUE) else 104
    size_max <- if ("HSIZE"    %in% names(df)) min(max(df$HSIZE, na.rm=TRUE), 20) else 20
    updateSliderInput(session, "f_decile",  value=c(1,10))
    updateSelectInput(session, "f_quintile",selected="all")
    updateSelectInput(session, "f_sex",     selected="all")
    updateSelectInput(session, "f_pop",     selected="all")
    updateSliderInput(session, "f_age",     value=c(15,age_max))
    updateSliderInput(session, "f_hsize",   value=c(1,size_max))
    updateSelectInput(session, "f_dwell",   selected="all")
    updateSelectInput(session, "f_status",  selected="all")
    updateSelectInput(session, "f_elec",    selected="all")
  })

  # ── Filtered data ───────────────────────────────────────────────────────────
  filt <- reactive({
    req(raw())
    df <- raw()

    # Income decile range
    if (!is.null(input$f_decile) && "INCOME_DECILE" %in% names(df))
      df <- df[df$INCOME_DECILE >= input$f_decile[1] &
               df$INCOME_DECILE <= input$f_decile[2], ]

    # Income quintile
    if (!is.null(input$f_quintile) && input$f_quintile != "all" &&
        "INCOME_QUINTILE" %in% names(df))
      df <- df[df$INCOME_QUINTILE == as.integer(input$f_quintile), ]

    # Head sex
    if (!is.null(input$f_sex) && input$f_sex != "all" && "HEAD_SEX" %in% names(df))
      df <- df[df$HEAD_SEX == as.integer(input$f_sex), ]

    # Population group
    if (!is.null(input$f_pop) && input$f_pop != "all" &&
        "HEAD_POPULATION" %in% names(df))
      df <- df[df$HEAD_POPULATION == as.integer(input$f_pop), ]

    # Head age
    if (!is.null(input$f_age) && "HEAD_AGE" %in% names(df))
      df <- df[df$HEAD_AGE >= input$f_age[1] & df$HEAD_AGE <= input$f_age[2], ]

    # Household size
    if (!is.null(input$f_hsize) && "HSIZE" %in% names(df))
      df <- df[df$HSIZE >= input$f_hsize[1] & df$HSIZE <= input$f_hsize[2], ]

    # Dwelling type
    if (!is.null(input$f_dwell) && input$f_dwell != "all" &&
        "IRD_MAIND" %in% names(df))
      df <- df[df$IRD_MAIND == as.integer(input$f_dwell), ]

    # Present status
    if (!is.null(input$f_status) && input$f_status != "all" &&
        "PRESENT_STATUS" %in% names(df))
      df <- df[df$PRESENT_STATUS == as.integer(input$f_status), ]

    # Electricity
    if (!is.null(input$f_elec) && input$f_elec != "all" &&
        "ENG_ACCESS" %in% names(df))
      df <- df[df$ENG_ACCESS == as.integer(input$f_elec), ]

    df
  })

  # ── Record bar ──────────────────────────────────────────────────────────────
  output$record_bar <- renderUI({
    req(raw())
    df  <- filt()
    n   <- nrow(df)
    tot <- nrow(raw())
    wt  <- if ("HHOLD_WGT" %in% names(df))
             format(round(sum(df$HHOLD_WGT, na.rm=TRUE)), big.mark=",")
           else NULL
    div(class = "record-bar",
      tags$span(
        style = "font-family:'Playfair Display',serif;font-size:14px;color:#3D2B1F;",
        "Filtered Households:"
      ),
      tags$span(class = "n-badge", format(n, big.mark=",")),
      tags$span(
        class = "weighted-note",
        paste0("of ", format(tot, big.mark=","), " total"),
        if (!is.null(wt)) paste0("  ·  ≈ ", wt, " weighted households")
      )
    )
  })

  # ── Helpers ─────────────────────────────────────────────────────────────────
  fmt_r <- function(x) paste0("R\u00A0", format(round(x), big.mark=","))

  stat_card <- function(label, value, style="default") {
    cls <- switch(style,
      dark  = "stat-card stat-card-dark",
      tan   = "stat-card stat-card-tan",
      "stat-card"
    )
    div(class = cls,
      div(class = "stat-label", label),
      div(class = "stat-value", value)
    )
  }

  # ── Stat cards ───────────────────────────────────────────────────────────────
  output$card1 <- renderUI({
    req(filt())
    stat_card("Mode Annual Income",
              fmt_r(density_mode(filt()$INCOME, na.rm=TRUE)))
  })
  output$card2 <- renderUI({
    req(filt())
    stat_card("Mode Annual Expenditure",
              fmt_r(density_mode(filt()$EXPENDITURE, na.rm=TRUE)), "dark")
  })
  output$card3 <- renderUI({
    req(filt())
    stat_card("Households (sample)",
              format(nrow(filt()), big.mark=","), "tan")
  })
  output$card4 <- renderUI({
    req(filt())
    stat_card("Mean Annual Income",
              fmt_r(mean(filt()$INCOME, na.rm=TRUE)))
  })
  output$card5 <- renderUI({
    req(filt())
    stat_card("Mean Annual Expenditure",
              fmt_r(mean(filt()$EXPENDITURE, na.rm=TRUE)), "dark")
  })
  output$card6 <- renderUI({
    req(filt())
    stat_card("Mode Income Per Capita",
              fmt_r(density_mode(filt()$INCOME_PCP, na.rm=TRUE)), "tan")
  })

  # ── Summary: Decile table ────────────────────────────────────────────────────
  output$tbl_decile <- renderTable({
    req(filt())
    filt() %>%
      group_by(Decile = INCOME_DECILE) %>%
      summarise(
        `Households` = n(),
        `Mode Income (R)` = round(density_mode(INCOME, na.rm=TRUE)),
        `Mean Income (R)`   = round(mean(INCOME,   na.rm=TRUE)),
        .groups = "drop"
      ) %>%
      mutate(across(c(`Mode Income (R)`,`Mean Income (R)`),
                    ~format(.x, big.mark=",")))
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="cccc",
  width="100%")

  # ── Summary: Population group table ─────────────────────────────────────────
  output$tbl_pop <- renderTable({
    req(filt())
    filt() %>%
      mutate(
        Group = dplyr::recode(as.character(HEAD_POPULATION),
                              !!!pop_lbl, .default="Other")
      ) %>%
      group_by(`Pop Group` = Group) %>%
      summarise(
        `N`                 = n(),
        `Mode Income (R)` = format(round(density_mode(INCOME,       na.rm=TRUE)), big.mark=","),
        `Mode Exp. (R)`   = format(round(density_mode(EXPENDITURE,  na.rm=TRUE)), big.mark=","),
        `Mode Inc/Cap (R)` = format(round(density_mode(INCOME_PCP,  na.rm=TRUE)), big.mark=","),
        .groups = "drop"
      )
  },
  striped=TRUE, hover=TRUE, bordered=TRUE, rownames=FALSE,
  align="lcccc",
  width="100%")

  # ── Data table ───────────────────────────────────────────────────────────────
  output$main_dt <- renderDT({
    req(filt())

    display_cols <- c("UQNO","HEAD_SEX","HEAD_AGE","HEAD_POPULATION","HSIZE",
                      "PRESENT_STATUS","ENG_ACCESS","INCOME","EXPENDITURE",
                      "INCOME_PCP","EXPENDITURE_PCP",
                      "INCOME_DECILE","INCOME_QUINTILE","HHOLD_WGT")

    avail <- display_cols[display_cols %in% names(filt())]
    df    <- filt()[, avail, drop=FALSE]

    # Decode factors
    if ("HEAD_SEX"        %in% names(df))
      df$HEAD_SEX        <- dplyr::recode(as.character(df$HEAD_SEX),        !!!sex_lbl,    .default="?")
    if ("HEAD_POPULATION" %in% names(df))
      df$HEAD_POPULATION <- dplyr::recode(as.character(df$HEAD_POPULATION), !!!pop_lbl,    .default="?")
    if ("PRESENT_STATUS"  %in% names(df))
      df$PRESENT_STATUS  <- dplyr::recode(as.character(df$PRESENT_STATUS),  !!!status_lbl, .default="?")
    if ("ENG_ACCESS"      %in% names(df))
      df$ENG_ACCESS      <- dplyr::recode(as.character(df$ENG_ACCESS),      !!!elec_lbl,   .default="?")

    names(df) <- gsub("_", " ", names(df))

    num_cols <- c("INCOME","EXPENDITURE","INCOME PCP","EXPENDITURE PCP")
    num_cols <- num_cols[num_cols %in% names(df)]

    dt <- datatable(df,
      rownames  = FALSE,
      filter    = "top",
      class     = "cell-border",
      options   = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "lftip",
        columnDefs = list(list(className="dt-center", targets="_all"))
      )
    )

    if (length(num_cols) > 0)
      dt <- formatCurrency(dt, num_cols,
                           currency="R ", interval=3, mark=",", digits=0)

    dt
  })

  # ── Plot: histogram ──────────────────────────────────────────────────────────
  output$plt_hist <- renderPlot({
    req(filt()); df <- filt()
    cap  <- quantile(df$INCOME, 0.95, na.rm=TRUE)
    med  <- density_mode(df$INCOME[df$INCOME <= cap], na.rm=TRUE)
    df2  <- df[df$INCOME <= cap & !is.na(df$INCOME), ]

    ggplot(df2, aes(x=INCOME)) +
      geom_histogram(fill="#A0522D", colour="#3D2B1F", bins=45, alpha=0.9) +
      geom_vline(xintercept=med, colour="#D4783E",
                 linewidth=1.3, linetype="dashed") +
      annotate("text", x=med*1.05, y=Inf, vjust=1.4, hjust=0,
               label=paste0("Mode\nR", format(round(med), big.mark=",")),
               colour="#D4783E", size=3, family="serif", fontface="bold") +
      scale_x_continuous(labels=label_dollar(prefix="R ", big.mark=",", scale=1e-3,
                                              suffix="k")) +
      scale_y_continuous(labels=comma) +
      labs(x="Annual Household Income", y="Households",
           subtitle="Capped at 95th percentile") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: decile bar ─────────────────────────────────────────────────────────
  output$plt_decile <- renderPlot({
    req(filt())
    filt() %>%
      group_by(Decile = factor(INCOME_DECILE)) %>%
      summarise(Mode=density_mode(INCOME, na.rm=TRUE), n=n(), .groups="drop") %>%
      ggplot(aes(x=Decile, y=Mode, fill=Decile)) +
      geom_col(colour="#3D2B1F", width=0.72, alpha=0.9) +
      geom_text(aes(label=paste0("R",format(round(Mode/1000), big.mark=","),"k")),
                vjust=-0.4, size=2.8, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=PAL, guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k"),
                         expand=expansion(mult=c(0,0.12))) +
      labs(x="Income Decile", y="Mode Annual Income") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: population group boxplot ───────────────────────────────────────────
  output$plt_pop_box <- renderPlot({
    req(filt()); df <- filt()
    cap <- quantile(df$INCOME, 0.95, na.rm=TRUE)
    df  %>%
      filter(INCOME <= cap) %>%
      mutate(Pop = dplyr::recode(as.character(HEAD_POPULATION),
                                 !!!pop_lbl, .default="Other")) %>%
      ggplot(aes(x=reorder(Pop, INCOME, FUN=density_mode), y=INCOME, fill=Pop)) +
      geom_boxplot(colour="#3D2B1F", outlier.colour="#D4783E",
                   outlier.size=0.6, alpha=0.85, width=0.6) +
      scale_fill_manual(
        values=c("Black African"="#A0522D","Coloured"="#D4783E",
                 "Indian/Asian"="#C4956A","White"="#8B7D6B","Other"="#6B4226"),
        guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k")) +
      coord_flip() +
      labs(x=NULL, y="Annual Income", subtitle="Capped at 95th percentile") +
      theme_mcm()
  }, bg="#F5E6D3")

  # ── Plot: income/capita by household size ────────────────────────────────────
  output$plt_hsize <- renderPlot({
    req(filt())
    filt() %>%
      filter(HSIZE <= 12, !is.na(INCOME_PCP)) %>%
      group_by(Size = HSIZE) %>%
      summarise(med=density_mode(INCOME_PCP, na.rm=TRUE), n=n(), .groups="drop") %>%
      ggplot(aes(x=Size, y=med)) +
      geom_area(fill="#A0522D30") +
      geom_line(colour="#A0522D", linewidth=1.4) +
      geom_point(aes(size=n), colour="#D4783E", fill="#FAF3EA",
                 shape=21, stroke=1.8) +
      scale_x_continuous(breaks=1:12) +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k")) +
      scale_size_continuous(range=c(3,11), name="N households") +
      labs(x="Household Size", y="Mode Income Per Capita",
           subtitle="Households with size 1–12 | point size = n households") +
      theme_mcm() +
      theme(legend.position="bottom")
  }, bg="#F5E6D3")

  # ── Plot: status bar ─────────────────────────────────────────────────────────
  output$plt_status <- renderPlot({
    req(filt())
    lvls <- c("Wealthy","Very comfortable","Reasonably comfortable",
              "Just getting along","Poor","Very poor")
    clrs <- c("#3D2B1F","#6B4226","#A0522D","#C4956A","#D4783E","#EDD9C0")

    filt() %>%
      mutate(Status = dplyr::recode(as.character(PRESENT_STATUS),
                                    !!!status_lbl, .default="Unknown")) %>%
      filter(Status %in% lvls) %>%
      group_by(Status) %>%
      summarise(med=density_mode(INCOME, na.rm=TRUE), n=n(), .groups="drop") %>%
      mutate(Status = factor(Status, levels=lvls)) %>%
      ggplot(aes(x=Status, y=med, fill=Status)) +
      geom_col(colour="#3D2B1F", width=0.65, alpha=0.92) +
      geom_text(aes(label=paste0("n=", format(n, big.mark=","))),
                vjust=-0.5, size=3, colour="#3D2B1F", family="serif") +
      scale_fill_manual(values=setNames(clrs, lvls), guide="none") +
      scale_y_continuous(labels=label_dollar(prefix="R ", big.mark=",",
                                              scale=1e-3, suffix="k"),
                         expand=expansion(mult=c(0,0.14))) +
      labs(x=NULL, y="Mode Annual Income") +
      theme_mcm() +
      theme(axis.text.x=element_text(angle=22, hjust=1, size=10))
  }, bg="#F5E6D3")
}

shinyApp(ui, server)
