

if (FALSE) { # Run once, clean up spreadsheet

  val <- readr::read_csv("tables/Sarita_outcomes_sim_year_app.csv") %>%
    mutate(var_name = variable)
  table(val$var_name)

  val$var_name[grepl("Egg", val$variable)] <- "Egg (10^6)"

  val$var_name[grepl("Smolt_NOS|Smolt_HOS", val$variable)] <- "Outmigrating juvenile (10^5)"
  val$var_name[grepl("Smolt_Rel", val$variable)] <- "Outmigrating juvenile (10^5)"

  val$var_name[grepl("KPT", val$variable)] <- "Preterminal catch"
  val$var_name[grepl("Return", val$variable)] <- "Mature Return"
  val$var_name[grepl("KT", val$variable)] <- "Terminal catch"

  val$var_name[grepl("Escapement", val$variable)] <- "In-river Return"


  val$var_name[grepl("NOB|HOB", val$variable)] <- "Brood"

  val$var_name[grepl("IR_Catch", val$variable)] <- "In-river Catch"

  val$var_name[val$variable == "NOS" | val$variable == "HOS"] <- "Natural Spawners"

  val$Origin <- "Hatchery origin"
  val$Origin[grepl(".NO.", val$variable)] <- "Natural origin"
  val$Origin[val$variable == "Smolt_HOS"] <- "Natural origin"
  val$Origin[val$variable == "NOB"] <- "Natural origin"
  val$Origin[val$variable == "NOS"] <- "Natural origin"

  val$value[grepl("10\\^5", val$var_name)] <- val$value[grepl("10\\^5", val$var_name)]/1e5
  val$value[grepl("10\\^6", val$var_name)] <- val$value[grepl("10\\^6", val$var_name)]/1e6

  # State variables with NO and HO
  var_plot <- c("Egg (10^6)", "Outmigrating juvenile (10^5)", #"Hatchery Release (10^5)",
                "Preterminal catch", "Mature Return", "Terminal catch", "In-river Return",
                "Brood", "In-river Catch", "Natural Spawners")

  val_plot <- filter(val, var_name %in% var_plot) %>%
    summarise(value = sum(value), .by = c(Simulation, Year, Scenario, scenario, var_name, Origin))
  val_array <- reshape2::acast(
    val_plot,
    list("Simulation", "Year", "Scenario", "scenario", "var_name", "Origin"),
    fill = NA_real_,
    value.var = "value"
  )
  names(dimnames(val_array)) <- c("Simulation", "Year", "Scenario", "scenario", "var_name", "Origin")

  # Hatchery-specific state variables with total only
  var_plot_hatchery <- c("PNI", "pHOS_effective", "Outmigrating juvenile (10^5)")

  val_plot_hatchery <- filter(val, var_name %in% var_plot_hatchery, Origin == "Hatchery origin") %>%
    summarise(value = sum(value), .by = c(Simulation, Year, Scenario, scenario, var_name))

  val_plot_hatchery$var_name[val_plot_hatchery$var_name == "Outmigrating juvenile (10^5)"] <- "Releases"
  val_plot_hatchery$var_name[val_plot_hatchery$var_name == "pHOS_effective"] <- "pHOSeff"

  val_array_hatchery <- reshape2::acast(
    val_plot_hatchery,
    list("Simulation", "Year", "Scenario", "scenario", "var_name"),
    fill = NA_real_,
    value.var = "value"
  )
  names(dimnames(val_array_hatchery)) <- c("Simulation", "Year", "Scenario", "scenario", "var_name")

  output <- list(
    State = val_array,
    Hatchery = val_array_hatchery
  )
  saveRDS(output, "tables/Sarita_app_results.rds")

}



var_plot <- c("Egg (10^6)", "Outmigrating juvenile (10^5)", #"Hatchery Release (10^5)",
              "Preterminal catch", "Mature Return", "Terminal catch", "In-river Return",
              "Brood", "In-river Catch", "Natural Spawners")

output <- readRDS("tables/Sarita_app_results.rds")
val <- reshape2::melt(output$State) %>%
  mutate(var_name = factor(var_name, var_plot))

OM_plot <- "A. 10% freshwater survival"
MP_plot <- c("(1) IRER1300 = 0.25, pNOB = 0.5", "(9) IRER1300 = 0.75, pNOB = 1")
sim <- 1:3

# Line figure
source("99-Sarita-results-functions.R")

for (i in 1:length(MP_plot)) {

  #### Lineplot sum NO + HO
  val_plot <- filter(
    val,
    Scenario == OM_plot,
    scenario == MP_plot[i],
    #Year < 30,
    ifelse(var_name == "Outmigrating juvenile (10^5)", Year > 1, Year < 30)
  ) %>%
    summarise(value = sum(value, na.rm = TRUE), .by = c(Year, var_name, Scenario, scenario, Simulation))

  # What causes population to crash in subset of simulations
  if (FALSE && i == 2) {
    sim_low <- val_plot %>% filter(var_name == "Egg (10^6)", Year == 29, value < 0.5) %>% pull(Simulation)
    SMSE <- readRDS("SMSE/Sarita9.rds")

    surv <- SMSE@Misc$SOM@Habitat[[1]]@fry_sdev
    matplot(surv[sim_low, ], type = "l", lty = 1, col = 2, ylim = c(0, 1))
    matlines(surv[-sim_low, ], type = "l", lty = 1, col = 1)
  }

  # All simulations spaghetti plot
  g <- val_plot %>%
    plot_spaghetti(OM_name = OM_plot, MP_name = MP_plot[i])
  ggsave(paste0("figures/SMSE/scenario/lineplot_allsims_", i, ".png"), height = 6, width = 6)

  #### Lineplot separated by NO + HO
  val_plot_origin <- filter(
    val,
    Scenario == OM_plot,
    scenario == MP_plot[i],
    #Year < 30,
    ifelse(var_name == "Outmigrating juvenile (10^5)", Year > 1, Year < 30)
  )

  # All simulations spaghetti plot
  g <- val_plot_origin %>%
    plot_spaghetti(OM_name = OM_plot, MP_name = MP_plot[i], by_origin = TRUE, alpha = 0.1)
  ggsave(paste0("figures/SMSE/scenario/lineplot_allsims_NOHO_", i, ".png"), height = 6, width = 6)

  # Subset of three simulations
  g <- val_plot %>%
    plot_spaghetti(sims = sim, OM_name = OM_plot, MP_name = MP_plot[i])
  ggsave(paste0("figures/SMSE/scenario/lineplot_3sim_", i, ".png"), height = 6, width = 6)

  #### Bar figure - individual simulations
  for (x in sim) {
    val_plot_sim <- filter(
      val,
      Scenario == OM_plot,
      scenario == MP_plot,
      Simulation == x,
      ifelse(var_name == "Outmigrating juvenile (10^5)", Year > 1, Year < 30)
    ) %>%
      summarise(value = sum(value), .by = c(Year, var_name, Origin, Scenario, scenario))

    g <- ggplot(val_plot_sim, aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_col(colour = "black", linewidth = 0.1, width = 1, alpha = 0.75, aes(fill = Origin)) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, fill = NULL, title = OM_plot, subtitle = MP_plot, caption = paste("Simulation", x))
    ggsave(paste0("figures/SMSE/scenario/barplot_scenario_", i, "_sim_", x, ".png"), height = 6, width = 6)
  }
}



