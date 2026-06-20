library(shiny)
library(bslib)
library(readr)
library(DT)

# --- GLOBAL SCOPE: Load data once ---
# This runs only when the app starts
households_data <- read_csv("C:/Users/arabe/OneDrive/Documents/Research_Project/Fact_IES2023_Households.csv")

#Mid-Century Modern Theme
my_theme <- bs_theme(
  version = 5,
  bootswatch = "lux",
  primary = "#C5705D",
  secondary = "#4D6B6E",
  bg = "#F9F7F2",
  fg = "#2F2F2F",
  "input-border-color" = "#D6C6A0"
)

ui <- page_navbar(
  title = "Analytics",
  theme = my_theme,
  nav_panel("Data View",
            layout_sidebar(
              sidebar = sidebar(
                helpText("Data source: Fact_IES2023_Households.csv"),
                helpText("File loaded into memory at startup.")
              ),
              card(
                card_header("Household Dataset"),
                DTOutput("table")
              )
            )
  ),
  tags$head(
    tags$style(HTML("
      /* Force container to stack elements */
      .dataTables_scroll {
        display: flex !important;
        flex-direction: column !important;
      }
      /* Ensure header is on top and scrollable */
      .dataTables_scrollHead {
        order: 1 !important;
        overflow-x: auto !important;
        overflow-y: hidden !important;
      }
      /* Ensure body is below and scrollable */
      .dataTables_scrollBody {
        order: 2 !important;
        overflow-x: auto !important;
      }
    "))
  ),
  tags$script(HTML("
    $(document).on('shiny:connected', function() {
      // Use a small delay to ensure DT is fully initialized
      setTimeout(function() {
        var scrollHead = $('.dataTables_scrollHead');
        var scrollBody = $('.dataTables_scrollBody');

        scrollHead.on('scroll', function() {
          scrollBody.scrollLeft($(this).scrollLeft());
        });

        scrollBody.on('scroll', function() {
          scrollHead.scrollLeft($(this).scrollLeft());
        });
      }, 500); 
    });
  "))
)

server <- function(input, output, session) {
  
  # Simply reference the globally loaded object
  output$table <- renderDT({
    datatable(
      households_data, 
      options = list(
        pageLength = 10, 
        scrollX = TRUE,
        scrollCollapse = TRUE
        # Do NOT add scrollY here if you want the table to expand naturally, 
        # or keep it if you want a fixed height container.
      )
    )
  })
}

shinyApp(ui, server)
