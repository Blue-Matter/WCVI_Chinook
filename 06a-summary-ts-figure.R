



val <- readr::read_csv("tables/Sarita_outcomes_sim_year2.csv") %>%
  mutate(var_name = variable)
table(val$var_name)

val$var_name[grepl("Egg", val$variable)] <- "Egg (10^6)"

val$var_name[grepl("Smolt_NOS|Smolt_HOS", val$variable)] <- "Outmigrating juvenile (10^5)"
val$var_name[grepl("Smolt_Rel", val$variable)] <- "Outmigrating juvenile (10^5)"

val$var_name[grepl("KPT", val$variable)] <- "Preterminal catch"
val$var_name[grepl("Return", val$variable)] <- "Return"
val$var_name[grepl("KT", val$variable)] <- "Terminal catch"

val$var_name[grepl("Escapement", val$variable)] <- "Escapement"


val$var_name[grepl("NOB|HOB", val$variable)] <- "Brood"

val$var_name[grepl("Catch", val$variable)] <- "In-river catch"

val$var_name[val$variable == "NOS" | val$variable == "HOS"] <- "Spawners"

val$Origin <- "Hatchery origin"
val$Origin[grepl(".NO.", val$variable)] <- "Natural origin"
val$Origin[val$variable == "Smolt_HOS"] <- "Natural origin"
val$Origin[val$variable == "NOB"] <- "Natural origin"
val$Origin[val$variable == "NOS"] <- "Natural origin"


val$value[grepl("10\\^5", val$var_name)] <- val$value[grepl("10\\^5", val$var_name)]/1e5
val$value[grepl("10\\^6", val$var_name)] <- val$value[grepl("10\\^6", val$var_name)]/1e6

var_plot <- c("Egg (10^6)", "Outmigrating juvenile (10^5)", #"Hatchery Release (10^5)",
              "Preterminal catch", "Return", "Terminal catch", "Escapement", "Brood", "In-river catch", "Spawners")

OM_plot <- "A. Recent ocean ER"
#MP_plot <- "(1) ER = 0.5, pNOB = 0.5"
MP_plot <- "(9) ER = 1, pNOB = 1"
sim <- 1:3

# Quantile figures (all simulations)
#val_plot <- filter(
#  val,
#  Scenario == OM_plot,
#  scenario == MP_plot,
#  Year < 30
#) %>%
#  filter(var_name %in% var_plot) %>%
#  summarise(median = median(value),
#            lower = quantile(value, 0.025),
#            upper = quantile(value, 0.975),
#            .by = c(Year, var_name, Scenario, scenario)) %>%
#  mutate(var_name = factor(var_name, var_plot))
#
#g <- ggplot(val_plot, aes(Year, median)) +
#  facet_wrap(vars(var_name), scales = "free_y") +
#  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.75) +
#  geom_line() +
#  expand_limits(y = 0) +
#  theme(strip.background = element_blank(), legend.position = "bottom") +
#  labs(y = NULL, colour = "Simulation") +
#  ggtitle(OM_plot, subtitle = MP_plot)
#ggsave("figures/SMSE/example_quantileplot_ts.png", height = 6, width = 6)



# Line figure
val_plot <- filter(
  val,
  Scenario == OM_plot,
  scenario == MP_plot,
  Simulation %in% sim,
  Year < 30
) %>%
  filter(var_name %in% var_plot) %>%
  summarise(value = sum(value), .by = c(Year, var_name, Scenario, scenario, Simulation)) %>%
  mutate(var_name = factor(var_name, var_plot))

g <- ggplot(val_plot, aes(Year, value, colour = factor(Simulation))) +
  facet_wrap(vars(var_name), scales = "free_y") +
  geom_line() +
  expand_limits(y = 0) +
  theme(strip.background = element_blank(), legend.position = "bottom") +
  labs(y = NULL, colour = "Simulation") +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle(OM_plot, subtitle = MP_plot)
ggsave("figures/SMSE/example_lineplot_ts.png", height = 6, width = 6)

# Bar figure
for (x in sim) {
  val_plot <- filter(
    val,
    Scenario == OM_plot,
    scenario == MP_plot,
    Simulation == x,
    Year < 30
  ) %>%
    filter(var_name %in% var_plot) %>%
    summarise(value = sum(value), .by = c(Year, var_name, Origin, Scenario, scenario)) %>%
    mutate(var_name = factor(var_name, var_plot))

  g <- ggplot(val_plot, aes(Year, value)) +
    facet_wrap(vars(var_name), scales = "free_y") +
    geom_col(colour = "black", linewidth = 0.1, width = 1, alpha = 0.75, aes(fill = Origin)) +
    expand_limits(y = 0) +
    theme(strip.background = element_blank(), legend.position = "bottom") +
    labs(y = NULL, fill = NULL, title = OM_plot, subtitle = MP_plot, caption = paste("Simulation", x))
  ggsave(paste0("figures/SMSE/example_barplot_ts_sim_", x, ".png"), height = 6, width = 6)
}


