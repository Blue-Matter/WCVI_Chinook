

library(salmonMSE)
library(tidyverse)
library(reshape2)

source("99-Sarita-results-functions.R")

# Identify scenarios and management options ----
g_init <- expand.grid(
  ER = c(0.5, 0.75, 1),
  pNOB_target = c(0.5, 0.75, 1),
  ocean_ER_scalar = 1,
  surv = c("bootstrap", "high", "low"),
  fec = "constant"
)

gr <- rbind(
  g_init %>%
    mutate(Scenario = ifelse(surv == "bootstrap", "A. Recent ocean ER",
                             ifelse(surv == "high", "B. High freshwater survival", "C. Low freshwater survival"))),
  g_init %>% filter(surv == "bootstrap") %>% mutate(ocean_ER_scalar = 0.75) %>%
    mutate(Scenario = "D. Lower ocean ER"),
  g_init %>% filter(surv == "bootstrap") %>% mutate(ocean_ER_scalar = 0.75, fec = "decline") %>%
    mutate(Scenario = "E. Decline mat & fec")
) %>%
  mutate(
    Option = paste0("ER = ", ER, ", pNOB = ", pNOB_target),
    n = 1:n()
  )

nOM <- nrow(gr)

# Grab all model output ----
scenario_unique <- unique(gr$Scenario)

SMSE_list <- lapply(gr$n, function(i) {
  SMSE <- readRDS(file.path("SMSE", paste0("Sarita", i, ".rds")))
  return(SMSE)
})

# Performance metrics (all simulations at end of projection) ----
y <- 29

PNI <- sapply(SMSE_list, function(x) x@PNI[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, PNI = value, n = Var2)

NS <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  rowSums(NS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Natural Spawners` = value, n = Var2)

MA <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  MA <- apply(NS_a, 1, function(w) weighted.mean(x = 1:5, w = w))
  return(MA)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Mean age` = value, n = Var2)

p_wild <- sapply(SMSE_list, function(x) x@p_wild[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pWILD = value, n = Var2)

pHOSeff <- sapply(SMSE_list, function(x) x@pHOS_effective[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pHOSeff = value, n = Var2)

pNOBeff <- sapply(SMSE_list, function(x) x@pNOB[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pNOBeff = value, n = Var2)

Brood <- sapply(SMSE_list, function(x) x@NOB[, , y] + x@HOB[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Brood = value, n = Var2)

IRR <- sapply(SMSE_list, function(x) {
  apply(x@Escapement_NOS[, , , y] + x@Escapement_HOS[, , , y], 1, sum)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, IR_Return = value, n = Var2)

IRCatch <- sapply(SMSE_list, function(x) {
  x@Misc$inriver_catch[, y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, IR_Catch = value, n = Var2)

Egg <- sapply(SMSE_list, function(x) x@Egg_NOS[, 1, y] + x@Egg_HOS[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Egg = value, n = Var2) %>%
  mutate(Egg = Egg/1e6)

Rel <- sapply(SMSE_list, function(x) x@Smolt_Rel[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Releases = value, n = Var2) %>%
  mutate(Releases = Releases/1e5)

Option_name <- data.frame(
  Option = gr$Option[1:9]
) %>%
  mutate(scenario = paste0("(", 1:9, ") ", Option))

# All in one data frame
val_sim_all <- list(PNI, NS, pNOBeff, pHOSeff, p_wild, MA,
                    Brood, IRCatch, IRR, Egg, Rel) %>%
  Reduce(left_join, .) %>%
  left_join(gr %>% select(Scenario, Option, n), by = "n") %>%
  reshape2::melt(id.vars = c("Simulation", "Option", "Scenario", "n")) %>%
  left_join(Option_name)
readr::write_csv(val_sim_all, file = "tables/Sarita_outcomes_sim.csv") # Save for Slick object

val_sim <- val_sim_all %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c("Option", "Scenario", "n", "variable")) %>%
  left_join(Option_name, by = "Option") %>%
  mutate(scenario = factor(scenario, rev(Option_name$scenario)))

# Probability (end of projection) ----
P_Sgen_NOS <- sapply(SMSE_list, function(x, Sgen = 250) {
  val <- rowSums(x@NOS[, 1, , y]) >= Sgen
  mean(val)
})

P_Smsy85_NOS <- sapply(SMSE_list, function(x, SMSY = 560) {
  val <- rowSums(x@NOS[, 1, , y]) >= 0.85 * SMSY
  mean(val)
})

P_Sgen <- sapply(SMSE_list, function(x, Sgen = 250) {
  val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >= Sgen
  mean(val)
})

P_Smsy85 <- sapply(SMSE_list, function(x, SMSY = 560) {
  val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >= 0.85 * SMSY
  mean(val)
})

P_1300 <- sapply(SMSE_list, function(x) {
  val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >= 1300
  mean(val)
})

P_PNI50 <- sapply(SMSE_list, function(x) {
  NOB <- x@NOB[, 1, y]
  PNI <- x@PNI[, 1, y]
  PNI[is.na(PNI) & NOB == 0] <- 0
  mean(PNI >= 0.5)
})

val_prob <- data.frame(n = 1:nrow(gr)) %>%
  mutate(
    P_PNI50 = P_PNI50,
    P_250_NOS = P_Sgen_NOS,
    P_476_NOS = P_Smsy85_NOS,
    P_250_NS = P_Sgen,
    P_476_NS = P_Smsy85,
    P_1300_NS = P_1300
  ) %>%
  reshape2::melt(id.vars = "n") %>%
  left_join(gr) %>%
  left_join(Option_name)
readr::write_csv(val_prob, file = "tables/Sarita_outcomes_prob.csv") # Save for Slick object

# Big data frame of state variables for each simulation and year (for Slick)
state_var <- c("PNI", "Total Spawners", "pNOB", "pHOS_effective", "p_wild", "Mean age", "Brood", "IR_Catch", "IR_Return", "Egg", "Smolt_Rel")
state_names <- state_var
state_names[2] <- "Natural Spawners"
state_names[3] <- "pNOBeff"
state_names[4] <- "pHOSeff"

df <- lapply(1:length(state_var), function(j) {
  lapply(1:length(SMSE_list), function(i) {
    .ts_fn(SMSE_list[[i]], var = state_var[j], all_sims = TRUE) %>%
      reshape2::melt() %>%
      rename(Simulation = Var1, Year = Var2) %>%
      mutate(variable = state_names[j], n = gr$n[i])
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  left_join(select(gr, Scenario, Option, n), by = "n") %>%
  left_join(Option_name)
readr::write_csv(df, file = "tables/Sarita_outcomes_sim_year.csv") # Save for Slick object
rm(df)

# Big data frame of state variables for each simulation and year (for salmonMSE app figure)
state_var2 <- c("Egg_NOS", "Egg_HOS", "Smolt_Rel", "Smolt_NOS", "Smolt_HOS", "KPT_NOS", "KPT_HOS", "Return_NOS", "Return_HOS",
                "KT_NOS", "KT_HOS",
                "Escapement_NOS", "Escapement_HOS", "NOB", "HOB", "IR_Catch", "NOS", "HOS")

df <- lapply(1:length(state_var2), function(j) {
  lapply(1:length(SMSE_list), function(i) {
    .ts_fn(SMSE_list[[i]], var = state_var2[j], all_sims = TRUE) %>%
      reshape2::melt() %>%
      rename(Simulation = Var1, Year = Var2) %>%
      mutate(variable = state_var2[j], n = gr$n[i])
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  left_join(select(gr, Scenario, Option, n), by = "n") %>%
  left_join(Option_name)
readr::write_csv(df, file = "tables/Sarita_outcomes_sim_year2.csv") # Save for Slick object
rm(df)

#### Run loop over each of 4 scenarios for performance metric tables ----
pm_primary <- c("PNI", "Natural Spawners", "P_PNI50", "P_250_NOS", "P_476_NOS", "P_250_NS", "P_476_NS", "P_1300_NS")
pm_ancillary <- c("IR_Return", "IR_Catch", "Brood", "Egg", "Releases",
                  "pNOBeff", "pHOSeff", "pWILD") #"Mean age")

for (i in 1:length(scenario_unique)) {

  ind <- gr$Scenario == scenario_unique[i]

  dir <- file.path("figures", "SMSE", paste0("Set", i, "_"))

  #### Time series ribbon plots (annual medians + 95% intervals)
  #g1 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "PNI") +
  #  coord_cartesian(ylim = c(0, 1))
#
  #g2 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Total Spawners") +
  #  coord_cartesian(ylim = c(0, 2500)) +
  #  guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
#
  #g3 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "pHOS_effective") +
  #  coord_cartesian(ylim = c(0, 1)) +
  #  labs(y = expression(pHOS[eff]))
#
  #g4 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "p_wild") +
  #  coord_cartesian(ylim = c(0, 1)) +
  #  labs(y = "pWILD")
#
  #g5 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Brood") +
  #  coord_cartesian(ylim = c(0, 700))
#
  #g6 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "pNOB") +
  #  coord_cartesian(ylim = c(0, 1))
#
  #g7 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Egg") +
  #  coord_cartesian(ylim = c(0, 2e6)) +
  #  labs(y = "Egg production (natural environment)")
#
  #g8 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Smolt_Rel") +
  #  coord_cartesian(ylim = c(0, 5e5)) +
  #  labs(y = "Hatchery releases")
#
  #g <- ggpubr::ggarrange(g1, g2, g3, g4, g5, g6, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
  #ggsave(paste0(dir, "ts.png"), g, height = 7, width = 6)

  #### Time series barplots (annual medians for by scenario for each management option)
  ind2 <- which(ind)

  png(paste0(dir, "spawners_prop.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_spawners(SMSE_list[[jj]])
    box()
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "spawners.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_spawners(SMSE_list[[jj]], prop = FALSE, ylim = c(0, 4000))
    box()
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "p_brood.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_escapement(SMSE_list[[jj]], ylim = c(0, 1))
    box()
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "pHOS_fitness.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_fitness(SMSE_list[[jj]], ylim = c(0, 1))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "RS_HOS.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "HOS", type = "abs",
            name = c("Fed Fry", "Traditionals"), ylim = c(0, 4000))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "RS_Rel.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Fed Fry", "Traditionals"),
            ylab = "Hatchery releases")
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "RS_Esc.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Fed Fry", "Traditionals"),
            ylab = "HO Escapement",
            ylim = c(0, 5000))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_NOS.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "NOS", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NOS",
             ylim = c(0, 700))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_Esc.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO Escapement",
             ylim = c(0, 900))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_Smolt.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO outmigrating juveniles",
             ylim = c(0, 4e5))
    title(Option_name$scenario[j])
  }
  dev.off()

  # Histograms across simulations of performance metrics at end of projection
  PNI_hist <- val_sim_all %>%
    filter(Scenario == scenario_unique[i]) %>%
    left_join(Option_name) %>%
    filter(variable == "PNI")
  g <- PNI_hist %>%
    plot_histogram() +
    ggtitle(scenario_unique[i]) +
    coord_cartesian(xlim = c(0, 1))
  ggsave(paste0(dir, "histogram_PNI.png"), g, width = 6, height = 4)

  ref <- data.frame(
    Lower = 250,
    Upper = 476,
    Srep = 1300
  ) %>%
    reshape2::melt()

  g <- val_sim_all %>%
    filter(Scenario == scenario_unique[i]) %>%
    left_join(Option_name) %>%
    filter(variable == "Natural Spawners") %>%
    plot_histogram("Natural Spawners", binwidth = 100) +
    ggtitle(scenario_unique[i]) +
    geom_vline(data = ref, aes(xintercept = value, linetype = variable)) +
    theme(legend.position = "bottom") +
    labs(linetype = "Benchmark") +
    scale_linetype_manual(values = c(2, 3, 4))
  ggsave(paste0(dir, "histogram_NS.png"), g, width = 6, height = 4.5)

  # Make figure of performance metrics (all simulations at end of projection) ----
  g <- val_sim %>%
    filter(Scenario == scenario_unique[i]) %>%
    left_join(gr) %>%
    #filter(variable %in% pm_primary) %>%
    plot_dotplot() +
    scale_shape_manual(values = c(1, 4, 16)) +
    #theme(strip.placement = "outside") +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    labs(x = NULL, y = NULL, shape = "ER", colour = "pNOB target") +
    ggtitle(scenario_unique[i]) +
    scale_x_discrete(labels = bold_scenario) +
    geom_vline(xintercept = c(3, 6) + 0.5)
  ggsave(paste0(dir, "performance_metrics.png"), g, width = 6, height = 7)

  # Performance metrics (medians at end of the projection) ----
  d <- rbind(
    val_sim %>% select(Scenario, variable, median, scenario),
    val_prob %>% select(Scenario, variable, value, scenario) %>% rename(median = value)
  ) %>%
    filter(Scenario == scenario_unique[i],
           variable %in% c(pm_primary, pm_ancillary)) %>%
    mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

  g <- plot_table(d) +
    geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
    ggtitle(scenario_unique[i]) +
    scale_y_discrete(labels = bold_scenario) +
    geom_hline(yintercept = c(3, 6) + 0.5, linewidth = 1)
  ggsave(paste0(dir, "performance_table_full.png"), g, width = 7.5, height = 3)
}

# Full performance table all scenarios
d <- rbind(
  val_sim %>% select(Scenario, variable, median, scenario),
  val_prob %>% select(Scenario, variable, value, scenario) %>% rename(median = value)
) %>%
  filter(variable %in% c(pm_primary, pm_ancillary)) %>%
  mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))


g <- plot_table(d) +
  geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
  facet_wrap(vars(Scenario), ncol = 1) +
  scale_y_discrete(labels = bold_scenario) +
  geom_hline(yintercept = c(3, 6) + 0.5, linewidth = 1)
ggsave("figures/SMSE/performance_table_full.png", g, width = 7.5, height = 9)

#### Decision table for all scenarios ----
y <- 29

dir_dt <- file.path("figures", "SMSE")

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "PNI") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid(
    ncol = 2,
    title = "Median PNI",
    fill_scheme =
      scale_fill_gradientn(
        values = c(0, 0.7, 1),
        colours = c("deeppink", "lightgreen", "green4")
      )
  )
ggsave(file.path(dir_dt, "decisiontable_PNI.png"), g, width = 4, height = 5)

val <- seq(0, 1, 0.01)
cols <- ifelse(val >= 0.5, "lightgreen", "deeppink")
g <- val_prob %>%
  left_join(gr) %>%
  filter(variable == "P_PNI50") %>%
  select(Scenario, ER, pNOB_target, value) %>%
  decision_table_grid(
    ncol = 2,
    title = "Probability PNI > 0.5",
    fill_scheme =
      scale_fill_gradientn(
        values = val,
        colours = cols
      )
  )
ggsave(file.path(dir_dt, "decisiontable_PNI50.png"), g, width = 4, height = 5)

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "Natural Spawners") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  mutate(value = round(median)) %>%
  decision_table_grid(ncol = 2, title = "Natural Spawners")
ggsave(file.path(dir_dt, "decisiontable_sp.png"), g, width = 4, height = 5)

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "Releases") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid(ncol = 2, "Hatchery releases (100,000s)")
ggsave(file.path(dir_dt, "decisiontable_rel.png"), g, width = 4, height = 5)

g <- val_prob %>%
  left_join(gr) %>%
  filter(variable == "P_1300_NS") %>%
  select(Scenario, ER, pNOB_target, value) %>%
  decision_table_grid(
    ncol = 2,
    title = "Probabilty > 1300 natural spawners",
    fill_scheme =
      scale_fill_gradientn(
        values = c(0, 0.7, 1),
        colours = c("deeppink", "lightgreen", "green4")
      )
  )
ggsave(file.path(dir_dt, "decisiontable_P_1300.png"), g, width = 4, height = 5)

# Tradeoff figure
g <- val_sim %>%
  left_join(gr) %>%
  tradeoff_grid(xname = "Natural Spawners", yname = "PNI", xlim = c(0, 4000), ylim = c(0, 1), ncol = 2) +
  geom_vline(xintercept = 1300, linetype = 3) +
  geom_hline(yintercept = 0.5, linetype = 3)
ggsave(file.path(dir_dt, "tradeoff_PNI_sp.png"), g, width = 5, height = 5)


g <- val_sim %>%
  left_join(gr) %>%
  tradeoff_grid(xname = "Releases", yname = "PNI", xlim = c(0, 5), ylim = c(0, 1),
                ncol = 2, xlab = "Hatchery releases (100,000s)") +
  geom_hline(yintercept = 0.5, linetype = 3)
ggsave(file.path(dir_dt, "tradeoff_PNI_rel.png"), g, width = 5, height = 5)

