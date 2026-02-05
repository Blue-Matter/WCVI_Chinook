

library(salmonMSE)
library(tidyverse)
library(reshape2)

source("99-Sarita-results-functions.R")

# Identify scenarios and management options ----
gr <- readr::read_csv(file.path("tables", "Sarita_scenarios.csv"))
nOM <- nrow(gr)

# Grab all model output ----
scenario_unique <- unique(gr$Scenario) # Represented by individual table

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
  x@Misc$inriver_catch$NOS[, y] + x@Misc$inriver_catch$HOS[, y]
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

#Option_name <- data.frame(
#  Option = gr$Option[1:9]
#) %>%
#  mutate(scenario = paste0("(", 1:9, ") ", Option))

# State variables for all simulations in year y in one data frame
val_sim_all <- list(PNI, NS, pNOBeff, pHOSeff, p_wild, MA,
                    Brood, IRCatch, IRR, Egg, Rel) %>%
  Reduce(left_join, .) %>%
  left_join(gr %>% select(Scenario_name, Option_name, n), by = "n") %>%
  rename(Option = Option_name, Scenario = Scenario_name) %>%
  reshape2::melt(id.vars = c("Simulation", "Option", "Scenario", "n"))
readr::write_csv(val_sim_all, file = "tables/Sarita_outcomes_sim.csv") # Save for Slick object

# Median and range for state variables across simulations
val_sim <- val_sim_all %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c("Option", "Scenario", "n", "variable")) #%>%
  #left_join(Option_name, by = "Option") %>%
  #mutate(scenario = factor(scenario, rev(Option_name$scenario)))

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
  left_join(select(gr, Option_name, Scenario_name, n)) %>%
  rename(Option = Option_name, Scenario = Scenario_name)
readr::write_csv(val_prob, file = "tables/Sarita_outcomes_prob.csv") # Save for Slick object

# Big data frame of state variables for each simulation and year (for shinysalmon app)
# CSV file is 200 MB, but we'll have to convert the data frame to arrays and save as an R object (.rds) later to reduce disk space
state_var2 <- c("Egg_NOS", "Egg_HOS", "Smolt_Rel", "Smolt_NOS", "Smolt_HOS", "KPT_NOS", "KPT_HOS", "Return_NOS", "Return_HOS",
                "KT_NOS", "KT_HOS",
                "Escapement_NOS", "Escapement_HOS", "NOB", "HOB", "IR_Catch", "NOS", "HOS", "pNOB", "PNI", "pHOS_effective")

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
  left_join(select(gr, Scenario_name, Option_name, n), by = "n") %>%
  rename(Option = Option_name, Scenario = Scenario_name)
readr::write_csv(df, file = "tables/Sarita_outcomes_sim_year_app.csv") # Save for Slick object
rm(df)

#### Run loop over each scenario to create performance metric tables ----
pm_primary <- c("PNI", "Natural Spawners", "P_PNI50", "P_250_NOS", "P_476_NOS", "P_250_NS", "P_476_NS", "P_1300_NS")
pm_ancillary <- c("IR_Return", "IR_Catch", "Brood", "Egg", "Releases",
                  "pNOBeff", "pHOSeff", "pWILD") #"Mean age")

# Time series medians
for (i in 1:length(scenario_unique)) {

  ind <- gr$Scenario == scenario_unique[i]

  dir <- file.path("figures", "SMSE", paste0("Set", i, "_"))

  #### Time series barplots (annual medians for by scenario for each management option)
  ind2 <- which(ind)

  png(paste0(dir, "spawners_prop.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_spawners(SMSE_list[[jj]])
    box()
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "spawners.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_spawners(SMSE_list[[jj]], prop = FALSE, ylim = c(0, 4000))
    box()
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "p_brood.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_escapement(SMSE_list[[jj]], ylim = c(0, 1))
    box()
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "pHOS_fitness.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_fitness(SMSE_list[[jj]], ylim = c(0, 1))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "RS_HOS.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "HOS", type = "abs",
            name = c("Fed Fry", "Traditionals"), ylim = c(0, 4000))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "RS_Rel.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Fed Fry", "Traditionals"),
            ylab = "Hatchery releases")
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "RS_Esc.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_RS(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Fed Fry", "Traditionals"),
            ylab = "HO Escapement",
            ylim = c(0, 5000))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "LHG_NOS.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "NOS", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NOS",
             ylim = c(0, 700))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "LHG_Esc.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO Escapement",
             ylim = c(0, 900))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

  png(paste0(dir, "LHG_Smolt.png"), height = 6, width = 6, units = "in", res = 400)
  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
  for (j in 1:nrow(Option_name)) {
    jj <- ind2[j]
    plot_LHG(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Early Smalls", "Late Larges"),
             ylab = "NO outmigrating juveniles",
             ylim = c(0, 4e5))
    title(Option_name$scenario[j], cex.main = 0.75)
  }
  dev.off()

}


# Tables of performance metrics by freshwater survival
for (i in LETTERS[1:3]) {

  # Make figure of performance metrics (all simulations at end of projection) ----
  #g <- val_sim %>%
  #  filter(Scenario == scenario_unique[i]) %>%
  #  left_join(gr) %>%
  #  #filter(variable %in% pm_primary) %>%
  #  plot_dotplot() +
  #  scale_shape_manual(values = c(1, 4, 16)) +
  #  #theme(strip.placement = "outside") +
  #  theme(legend.position = "bottom") +
  #  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  #  labs(x = NULL, y = NULL, shape = "IRER 1300", colour = "pNOB target") +
  #  ggtitle(scenario_unique[i]) +
  #  scale_x_discrete(labels = bold_scenario) +
  #  geom_vline(xintercept = c(3, 6) + 0.5)
  #ggsave(paste0(dir, "performance_metrics.png"), g, width = 6, height = 7)


  # Full performance metrics (medians at end of the projection) ----
  d <- rbind(
    val_sim %>% select(Option, Scenario, variable, median),
    val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
  ) %>%
    filter(grepl(i, Scenario)) %>%
    filter(variable %in% c(pm_primary, pm_ancillary)) %>%
    mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

  g <- plot_table(d, ncol = 1) +
    geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
    scale_y_discrete(labels = font_fn, limits = rev) +
    geom_hline(yintercept = 3.5, linewidth = 1)
  ggsave(file.path("figures", "SMSE", paste0("performance_table_", i, ".png")), g, width = 7.5, height = 7)

  # Short performance metrics
  pm_short <- c("PNI", "Natural Spawners", "P_PNI50", "P_1300_NS", "IR_Catch")

  d <- rbind(
    val_sim %>% select(Option, Scenario, variable, median),
    val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
  ) %>%
    filter(grepl(i, Scenario)) %>%
    filter(variable %in% pm_short) %>%
    mutate(variable = factor(variable, pm_short))

  g <- plot_table(d, ncol = 1) +
    #geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
    scale_y_discrete(labels = font_fn, limits = rev) +
    geom_hline(yintercept = 3.5, linewidth = 1)
  ggsave(file.path("figures", "SMSE", paste0("performance_table_", i, ".png")), g, width = 7.5, height = 7)

}

# Full performance table all scenarios
glist <- lapply(LETTERS[1:3], function(i) {
  d <- rbind(
    val_sim %>% select(Option, Scenario, variable, median),
    val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
  ) %>%
    filter(grepl(i, Scenario)) %>%
    filter(variable %in% c(pm_primary, pm_ancillary)) %>%
    mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

  g <- plot_table(d, ncol = 1) +
    geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
    geom_hline(yintercept = 3.5, linewidth = 1)

  if (i == "A") {
    g <- g + scale_y_discrete(labels = font_fn, limits = rev)
  } else {
    g <- g + scale_y_discrete(labels = NULL, limits = rev)
  }

  g
})

g <- ggpubr::ggarrange(plotlist = glist, ncol = 3, widths = c(3, 2, 2))
ggsave("figures/SMSE/performance_table_full.png", g, width = 17, height = 8)

# Short performance metrics
pm_short <- c("PNI", "Natural Spawners", "P_PNI50", "P_1300_NS", "IR_Catch")

glist <- lapply(LETTERS[1:3], function(i) {
  d <- rbind(
    val_sim %>% select(Option, Scenario, variable, median),
    val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
  ) %>%
    filter(grepl(i, Scenario)) %>%
    filter(variable %in% pm_short) %>%
    mutate(variable = factor(variable, pm_short))

  g <- plot_table(d, ncol = 1) +
    geom_vline(xintercept = c(2.5, 4.5), linewidth = 1, linetype = 2) +
    geom_hline(yintercept = 3.5, linewidth = 1)

  if (i == "A") {
    g <- g + scale_y_discrete(labels = font_fn, limits = rev)
  } else {
    g <- g + scale_y_discrete(labels = NULL, limits = rev)
  }

  g
})

g <- ggpubr::ggarrange(plotlist = glist, ncol = 3, widths = c(3, 2, 2))
ggsave("figures/SMSE/performance_table_short.png", g, width = 8, height = 7)

#### Decision table for all scenarios ----
y <- 29

g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  filter(variable == "PNI") %>%
  select(Scenario, IRER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid(
    ncol = 3,
    title = "Median PNI",
    fill_scheme =
      scale_fill_gradientn(
        values = c(0, 0.7, 1),
        colours = c("deeppink", "lightgreen", "green4")
      )
  )
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
ggsave(file.path("figures", "SMSE", "decisiontable_PNI.png"), g, width = 7, height = 5)

val <- seq(0, 1, 0.01)
cols <- ifelse(val >= 0.5, "lightgreen", "deeppink")
g <- val_prob %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  filter(variable == "P_PNI50") %>%
  select(Scenario, IRER, pNOB_target, value) %>%
  decision_table_grid(
    ncol = 3,
    title = "Probability PNI > 0.5",
    fill_scheme =
      scale_fill_gradientn(
        values = val,
        colours = cols
      )
  )
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
ggsave(file.path("figures", "SMSE", "decisiontable_PNI50.png"), g, width = 7, height = 5)

g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  filter(variable == "Natural Spawners") %>%
  select(Scenario, IRER, pNOB_target, median) %>%
  mutate(value = round(median)) %>%
  decision_table_grid(ncol = 3, title = "Natural Spawners")
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
ggsave(file.path("figures", "SMSE", "decisiontable_sp.png"), g, width = 7, height = 5)

g <- val_prob %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  filter(variable == "P_1300_NS") %>%
  select(Scenario, IRER, pNOB_target, value) %>%
  decision_table_grid(
    ncol = 3,
    title = "Probabilty > 1300 natural spawners",
    fill_scheme =
      scale_fill_gradientn(
        values = c(0, 0.7, 1),
        colours = c("deeppink", "lightgreen", "green4")
      )
  )
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
ggsave(file.path("figures", "SMSE", "decisiontable_P_1300.png"), g, width = 7, height = 5)

g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  filter(variable == "Releases") %>%
  select(Scenario, IRER, pNOB_target, median) %>%
  rename(value = median) %>%
  decision_table_grid(ncol = 3, "Hatchery releases (100,000s)")
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
ggsave(file.path("figures", "SMSE", "decisiontable_rel.png"), g, width = 7, height = 5)

# Tradeoff figure
g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  tradeoff_grid(xname = "Natural Spawners", yname = "PNI", xlim = c(0, 12000), ylim = c(0, 1), ncol = 3) +
  geom_vline(xintercept = 1300, linetype = 3) +
  geom_hline(yintercept = 0.5, linetype = 3)
ggsave(file.path("figures", "SMSE", "tradeoff_PNI_sp.png"), g, width = 7, height = 5)

g <- val_prob %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  mutate(lwr = value, median = value, upr = value) %>%
  tradeoff_grid(xname = "P_1300_NS", yname = "P_PNI50", xlim = c(0, 1), ylim = c(0, 1), ncol = 3) +
  geom_vline(xintercept = 0.5, linetype = 3) +
  geom_hline(yintercept = 0.5, linetype = 3)
ggsave(file.path("figures", "SMSE", "tradeoff_prob.png"), g, width = 7, height = 5)

g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  tradeoff_grid(xname = "IR_Catch", yname = "PNI", xlim = c(0, 4000), ylim = c(0, 1), ncol = 3) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  labs(x = "In-river catch")
ggsave(file.path("figures", "SMSE", "tradeoff_PNI_IRC.png"), g, width = 7, height = 5)


g <- val_sim %>%
  left_join(select(gr, IRER, pNOB_target, n)) %>%
  tradeoff_grid(xname = "Releases", yname = "PNI", xlim = c(0, 5), ylim = c(0, 1),
                ncol = 3, xlab = "Hatchery releases (100,000s)") +
  geom_hline(yintercept = 0.5, linetype = 3)
ggsave(file.path(dir_dt, "tradeoff_PNI_rel.png"), g, width = 7, height = 5)


#### Plot Simulated CYER
CYER_PT <- calc_CYER(SMSE_list[[1]], PT = TRUE)
CYER_T <- calc_CYER(SMSE_list[[1]], PT = FALSE)

png("simulated_CYER.png", height = 6, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
matplot(t(CYER_PT), type = "l", lty = 1, xlab = "Projection Year", ylab = "Preterminal CYER")
matplot(t(CYER_T), type = "l", lty = 1, xlab = "Projection Year", ylab = "Terminal CYER")
dev.off()


#### Plot maturity and exploitation rate at age ----
png("figures/SMSE/maturity_ER.png", height = 6, width = 3, units = "in", res = 400)
par(mfrow = c(3, 1), mar = c(5, 4, 1, 1))

SOM <- SMSE_list[[1]]@Misc$SOM

salmonMSE:::plot_Mjuv_RS(SOM@Hatchery[[1]]@p_mature_HOS[, , 1, ],
                         RS_names = c("Fed Fry", "Traditionals"), ylab = "Proportion mature")

SOM@Harvest[[1]]@vulPT <- SMSE_list[[1]]@ExPT_NOS[, 1, , 30]

salmonMSE:::plot_SOM(SOM@Harvest[[1]], "vulPT",
                     type = "age", nsim = SOM@nsim, maxage = SOM@Bio[[1]]@maxage,
                     proyears = SOM@proyears,
                     ylab = "Preterminal exploitation rate")

SOM@Harvest[[1]]@vulT <- SMSE_list[[1]]@ExT_NOS[, 1, , 30]
salmonMSE:::plot_SOM(SOM@Harvest[[1]], "vulT",
                     type = "age", nsim = SOM@nsim, maxage = SOM@Bio[[1]]@maxage,
                     proyears = SOM@proyears,
                     ylab = "Terminal exploitation rate")

dev.off()
