

library(salmonMSE)

# Identify scenarios and management options ----
gr <- expand.grid(
  ER = c(0.5, 0.75, 1),
  pNOB_target = c(0.5, 0.75, 1),
  ocean_ER_scalar = 1,
  surv = c("avg", "high"),
  fec = "constant"
) %>%
  mutate(Scenario = ifelse(surv == "avg", "A. Base (high ocean ER)", "B. High freshwater survival"))

gr <- rbind(
  gr,
  gr %>% filter(surv == "avg") %>%
    mutate(ocean_ER_scalar = 0.75) %>%
    mutate(Scenario = "C. Lower ocean ER"),
  gr %>% filter(surv == "avg") %>%
    mutate(ocean_ER_scalar = 0.75, fec = "decline") %>%
    mutate(Scenario = "D. Lower ocean ER, decline mat & fec")
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

TS <- sapply(SMSE_list, function(x) {
  TS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  rowSums(TS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Total Spawners` = value, n = Var2)

MA <- sapply(SMSE_list, function(x) {
  TS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  MA <- apply(TS_a, 1, function(w) weighted.mean(x = 1:5, w = w))
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

Esc <- sapply(SMSE_list, function(x) {
  apply(x@Escapement_NOS[, , , y] + x@Escapement_HOS[, , , y], 1, sum)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Escapement = value, n = Var2)

K <- sapply(SMSE_list, function(x) {
  x@Misc$ESSR_catch[, y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Catch = value, n = Var2)

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
val_sim_all <- list(PNI, TS, pNOBeff, pHOSeff, p_wild, MA,
                    Brood, K, Esc, Egg, Rel) %>%
  Reduce(left_join, .) %>%
  left_join(gr %>% select(Scenario, Option, n), by = "n") %>%
  reshape2::melt(id.vars = c("Simulation", "Option", "Scenario", "n"))

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

val_prob <- data.frame(n = 1:nrow(gr)) %>%
  mutate(
    P_250_NOS = P_Sgen_NOS,
    P_476_NOS = P_Smsy85_NOS,
    P_250 = P_Sgen,
    P_476 = P_Smsy85,
    P_1300 = P_1300
  ) %>%
  reshape2::melt(id.vars = "n") %>%
  left_join(gr) %>%
  left_join(Option_name)


source("99-Sarita-results-functions.R")
pm_primary <- c("PNI", "Total Spawners", "P_250_NOS", "P_476_NOS", "P_250", "P_476", "P_1300")
pm_ancillary <- c("pNOBeff", "pHOSeff", "pWILD", "Mean age",
                  "Brood", "Catch", "Escapement", "Egg", "Releases")


for (i in 1:length(scenario_unique)) {

  ind <- gr$Scenario == scenario_unique[i]

  dir <- file.path("figures", "SMSE", paste0("Set", i, "_"))

  #### Time series ribbon plots (annual medians + 95% intervals)
  g1 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "PNI") +
    coord_cartesian(ylim = c(0, 1))

  g2 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Total Spawners") +
    coord_cartesian(ylim = c(0, 2500)) +
    guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

  g3 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "pHOS_effective") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = expression(pHOS[eff]))

  g4 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "p_wild") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = "pWILD")

  g5 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Brood") +
    coord_cartesian(ylim = c(0, 700))

  g6 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "pNOB") +
    coord_cartesian(ylim = c(0, 1))

  g7 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Egg") +
    coord_cartesian(ylim = c(0, 2e6)) +
    labs(y = "Egg production (natural environment)")

  g8 <- ts_fn(SMSE_list[ind], Option_name$scenario, var = "Smolt_Rel") +
    coord_cartesian(ylim = c(0, 5e5)) +
    labs(y = "Hatchery releases")

  g <- ggpubr::ggarrange(g1, g2, g3, g4, g5, g6, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
  ggsave(paste0(dir, "ts.png"), g, height = 7, width = 6)

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
    plot_spawners(SMSE_list[[jj]], prop = FALSE, ylim = c(0, 2000))
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
            name = c("Fed Fry", "Traditionals"), ylim = c(0, 2000))
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
            ylim = c(0, 3000))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_NOS.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "NOS", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NOS",
             ylim = c(0, 300))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_Esc.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO Escapement",
             ylim = c(0, 300))
    title(Option_name$scenario[j])
  }
  dev.off()

  png(paste0(dir, "LHG_Smolt.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO outmigrating juveniles",
             ylim = c(0, 3e5))
    title(Option_name$scenario[j])
  }
  dev.off()

  # Histograms across simulations of performance metrics at end of projection
  g <- val_sim_all %>%
    filter(Scenario == scenario_unique[i]) %>%
    left_join(Option_name) %>%
    filter(variable == "PNI") %>%
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
    filter(variable == "Total Spawners") %>%
    plot_histogram("Total Spawners", binwidth = 100) +
    ggtitle(scenario_unique[i]) +
    geom_vline(data = ref, aes(xintercept = value, linetype = variable)) +
    theme(legend.position = "bottom") +
    labs(linetype = "Benchmark") +
    scale_linetype_manual(values = c(2, 3, 4))
  ggsave(paste0(dir, "histogram_TS.png"), g, width = 6, height = 4.5)

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
    ggtitle(scenario_unique[i])
  ggsave(paste0(dir, "performance_metrics.png"), g, width = 6, height = 7)

  # Primary performance metrics (medians at end of the projection) ----
  d <- rbind(
    val_sim %>% select(Scenario, variable, median, scenario),
    val_prob %>% select(Scenario, variable, value, scenario) %>% rename(median = value)
  ) %>%
    filter(Scenario == scenario_unique[i]) %>%
    mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

  #d <- val_sim %>%
  #  filter(Scenario == scenario_unique[i]) %>%
  #  select(scenario, variable, median)

  g <- plot_table(d) +
    geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
    ggtitle(scenario_unique[i])
  ggsave(paste0(dir, "performance_table_full.png"), g, width = 7.5, height = 3)

  # Probability (end of projection) ----
  #d_prob <- val_prob %>%
  #  rename(median = value) %>%
  #  filter(Scenario == scenario_unique[i]) %>%
  #  select(scenario, variable, median) %>%
  #  mutate(variable = factor(variable),
  #         scenario = factor(scenario, rev(Option_name$scenario)))

  #g <- plot_table(d_prob) +
  #  ggtitle(scenario_unique[i])
  #ggsave(paste0(dir, "prob_table.png"), g, width = 5, height = 3.5)

}

#### Decision table for all scenarios ----

y <- 29

dir_dt <- file.path("figures", "SMSE")

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "PNI") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid()
ggsave(file.path(dir_dt, "decisiontable_PNI.png"), g, width = 5, height = 5)

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "Total Spawners") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid("Total Spawners")
ggsave(file.path(dir_dt, "decisiontable_sp.png"), g, width = 5, height = 5)

g <- val_sim %>%
  left_join(gr) %>%
  filter(variable == "Releases") %>%
  select(Scenario, ER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid("Hatchery releases (100,000s)")
ggsave(file.path(dir_dt, "decisiontable_rel.png"), g, width = 5, height = 5)

g <- val_prob %>%
  left_join(gr) %>%
  filter(variable == "P_1300") %>%
  select(Scenario, ER, pNOB_target, value) %>%
  decision_table_grid("P_1300")
ggsave(file.path(dir_dt, "decisiontable_P_1300.png"), g, width = 5, height = 5)


g <- val_sim %>%
  left_join(gr) %>%
  #mutate(lwr = median, upr = median) %>%
  tradeoff_grid(xlim = c(0, 4000), ylim = c(0, 1))
ggsave(file.path(dir_dt, "tradeoff_PNI_sp.png"), g, width = 5, height = 5)

