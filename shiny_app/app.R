
library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(reshape2)
library(here)
library(ggborderline)
library(bslib)
library(htmltools)

theme_set(theme_bw())

figure_height <- 600


plot_spaghetti <- function(x, sims, OM_name = NULL, MP_name = NULL, alpha = 0.4, by_origin = FALSE) {

  if (by_origin) {
    require(ggborderline)

    meds <- summarise(x, value = median(value), .by = c(Year, var_name, Scenario, scenario, Origin))

    g <- x %>%
      mutate(gr = paste(Simulation, Origin)) %>%
      ggplot(aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line(alpha = alpha, aes(colour = Origin, group = factor(gr))) +
      #geom_line(data = meds, colour = "black", aes(group = Origin), linewidth = 1.5) +
      ggborderline::geom_borderline(data = meds, aes(colour = Origin), bordercolour = "grey40", linewidth = 1) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, colour = NULL) +
      ggtitle(OM_name, subtitle = MP_name)

  } else if (missing(sims)) { # All simulations

    meds <- summarise(x, value = median(value), .by = c(Year, var_name, Scenario, scenario))

    g <- ggplot(x, aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line(alpha = alpha, colour = "grey40", aes(group = factor(Simulation))) +
      geom_line(data = meds, colour = "blue", linewidth = 1) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL) +
      ggtitle(OM_name, subtitle = MP_name)

  } else {

    val_plot <- dplyr::filter(x, Simulation %in% sims)
    g <- ggplot(val_plot, aes(Year, value, colour = factor(Simulation))) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line() +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, colour = "Simulation") +
      scale_colour_brewer(palette = "Dark2") +
      ggtitle(OM_name, subtitle = MP_name)
  }
  g
}


ui <- fluidPage(
  theme = bslib::bs_theme(),
  titlePanel(
    h2("Sarita visualization app"),
  ),

  mainPanel(
    h4("Filter results:"),
    inputPanel(
      selectInput("OM", label = "Scenario", choices = ""),
      selectInput("MP", label = "Management option", choices = "")
    ),

    conditionalPanel("output.Panel_ready == 1",
      tabsetPanel(

        tabPanel(
          h4("Spaghetti plot"),
          h5("All simulations"),
          plotOutput("spaghetti", height = figure_height),
          h5("All simulations, by origin"),
          plotOutput("spaghetti_origin", height = figure_height),
          h5("Subset of simulations"),
          inputPanel(
            pickerInput("Sim_spaghetti", label = "Simulation Number\n(up to 8 total)", multiple = TRUE, choices = "")
          ),
          plotOutput("spaghetti_subset", height = figure_height),
        ),

        tabPanel(
          h4("Barplots"),
          inputPanel(
            selectInput("Sim_barplot", label = "Simulation Number", choices = "")
          ),
          plotOutput("barplot", height = figure_height),
        )

      )

    )
  )
)

server <- function(input, output, session) {

  Panel_ready <- reactiveVal(0)
  output$Panel_ready <- reactive({Panel_ready()})
  outputOptions(output, "Panel_ready", suspendWhenHidden = FALSE)

  ### Grab array and convert to data frame, download if necessary ----
  if (dir.exists(here::here("tables")) && file.exists(here::here("tables/Sarita_app_results.rds"))) {
    val_array <- readRDS(here::here("tables/Sarita_app_results.rds"))
  } else {
    file_name <- file.path(tempdir(), "Sarita_app_results.rds")
    download.file(
      "https://github.com/Blue-Matter/WCVI_Chinook/raw/refs/heads/main/tables/Sarita_app_results.rds",
      file_name
    )
    val_array <- readRDS(file_name)
  }

  var_plot <- c("Egg (10^6)", "Outmigrating juvenile (10^5)", #"Hatchery Release (10^5)",
                "Preterminal catch", "Mature Return", "Terminal catch", "In-river Return",
                "Brood", "In-river Catch", "Natural Spawners")

  val <- reshape2::melt(val_array) %>%
    mutate(var_name = factor(var_name, var_plot))

  updateSelectInput(session, "OM", choices = unique(val$Scenario))
  updateSelectInput(session, "MP", choices = unique(val$scenario))

  updatePickerInput(session, "Sim_spaghetti", choices = unique(val$Simulation), selected = seq(1, min(val$Simulation, 3)))

  updateSelectInput(session, "Sim_barplot", choices = unique(val$Simulation), selected = unique(val$Simulation)[1])

  Panel_ready(1)

  val_plot_origin <- reactive({
    dplyr::filter(
      val,
      Scenario == input$OM,
      scenario == input$MP,
      ifelse(var_name == "Outmigrating juvenile (10^5)", Year > 1, Year < max(.data$Year))
    )
  })

  val_plot <- reactive({
    val_plot_origin() %>%
      summarise(value = sum(value, na.rm = TRUE), .by = c(Year, var_name, Scenario, scenario, Simulation))
  })

  output$spaghetti <- renderPlot({
    plot_spaghetti(val_plot(), OM_name = input$OM, MP_name = input$MP)
  }, res = 100, height = figure_height)

  output$spaghetti_origin <- renderPlot({
    plot_spaghetti(val_plot_origin(), OM_name = input$OM, MP_name = input$MP, by_origin = TRUE, alpha = 0.1)
  }, res = 100, height = figure_height)

  output$spaghetti_subset <- renderPlot({
    plot_spaghetti(val_plot(), sims = input$Sim_spaghetti, OM_name = input$OM, MP_name = input$MP)
  }, res = 100, height = figure_height)


  val_plot_origin_sim <- reactive({
    dplyr::filter(
      val_plot_origin(),
      Simulation == input$Sim_barplot
    )
  })

  output$barplot <- renderPlot({
    val_plot_origin_sim() %>%
      ggplot(aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_col(colour = "black", linewidth = 0.1, width = 1, alpha = 0.75, aes(fill = Origin)) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, fill = NULL, title = input$OM, subtitle = input$MP, caption = paste("Simulation", input$Sim_barplot))

  }, res = 100, height = figure_height)

}

shinyApp(ui, server)
