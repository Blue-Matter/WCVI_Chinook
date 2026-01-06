
# Make OM
library(tidyverse)
library(salmonMSE)

maxage <- 5
nsim <- 100
nyears <- 2
proyears <- 30
n_g <- 2

# Load exploitation rate model - Sarita (with Robertson CWT traditionals)
ERM_Sarita <- readRDS("CM/Sarita_RBT_CM_01.05.26_pesc0.75.rds")
report_RBT <- salmonMSE:::get_report(ERM_Sarita) # Removes warmup and thinned samples

# Take maturity average from the 6 most recent brood years (2013-2018)
matt_dev <- readRDS("CM/Sarita_maturity.rds")
matt_avg <- sapply(matt_dev, function(x) {
  apply(x$matt[seq(9, 14), , ], 2:3, mean)
}, simplify = "array")

#apply(matt_avg, 1:2, mean)

set.seed(24)
sim_samp <- sample(seq(1, length(report_RBT)), nsim)

### Natural mortality - NOS ----
Mjuv_NOS <- array(0, c(nsim, maxage, nyears + proyears, n_g))

# Survival from CTC 23-06 p. 9 for ages 2+
M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1))

# Age 1 value by life cycle group
# Life cycle groups only differ in age 1 survival (first year of life in marine environment) from life cycle spreadsheet
surv1 <- 1 - c(0.996, 0.95)
Mjuv_NOS[, 1, , 1] <- -log(surv1[1])
Mjuv_NOS[, 1, , 2] <- -log(surv1[2])

# Age 2-5
Mjuv_NOS[, 2, , ] <- M_CTC[2]
Mjuv_NOS[, 3, , ] <- M_CTC[3]
Mjuv_NOS[, 4, , ] <- M_CTC[4]
Mjuv_NOS[, 5, , ] <- M_CTC[5]

### Maturity ----
# Assume identical maturity from hypothesized fed fry maturity
p_mature <- array(0, c(nsim, maxage, nyears + proyears))
p_mature[] <- matt_avg[, 1, sim_samp] %>% t() %>%
  array(c(nsim, maxage, nyears + proyears))

fec_Sarita <- c(0, 1500, 3000, 3600, 4600) # See Res Doc, but set age 2 fecundity arbitrarily equal to half of age 3
p_female <- c(0, 0.01, 0.1, 0.55, 0.8)
fec <- fec_Sarita * p_female

Bio <- new(
  "Bio",
  maxage = maxage,
  n_g = n_g,
  p_LHG = c(0.95, 0.05),
  p_mature = p_mature,
  Mjuv_NOS = Mjuv_NOS,
  p_female = 1, # p_female specified in fecundity
  fec = fec,
  s_enroute = 1
)

### Harvest, fishery vulnerability ----
vulPT <- sapply(report_RBT[sim_samp], getElement, "vulPT")
vulT <- sapply(report_RBT[sim_samp], getElement, "vulT")

#sapply(report_RBT, getElement, "vulPT") %>% apply(1, mean)
#sapply(report_RBT, getElement, "vulT") %>% apply(1, mean)

#matplot(vulPT, typ = 'l')
#matplot(vulT, typ = 'l')
#g <- salmonMSE:::CM_ER(report_RBT, brood = FALSE, type = "PT", year1 = 1979, at_age = FALSE)
#g <- salmonMSE:::CM_ER(report_RBT, brood = FALSE, type = "T", year1 = 1979, at_age = FALSE)
#g <- salmonMSE:::CM_ER(report_RBT, brood = FALSE, type = "all", year1 = 1979, at_age = FALSE)

Harvest <- new(
  "Harvest",
  type_PT = "u",
  type_T = "u",
  u_preterminal = 0.45,
  u_terminal = 0.24,
  MSF_PT = FALSE,
  MSF_T = FALSE,
  release_mort = c(0.1, 0.1),
  vulPT = t(vulPT),
  vulT = t(vulT)
)


### Hatchery ----
n_r <- 2

# Natural mortality HOS
Mjuv_HOS <- array(0, c(nsim, maxage, nyears + proyears, n_r))

# Age 1 survival
#g <- salmonMSE:::CM_M(report_RBT)
#g <- salmonMSE:::CM_M(report_Q)

surv1_HOS <- 1 - c(0.969, 0.969) # 0.969 from life cycle table, assuming traditional survival is 2x of small
Mjuv_HOS[, 1, , 1] <- -log(surv1_HOS[1])
Mjuv_HOS[, 1, , 2] <- -log(surv1_HOS[2])

# Age 2 survival
Mjuv_HOS[, 2, , ] <- M_CTC[2]
Mjuv_HOS[, 3, , ] <- M_CTC[3]
Mjuv_HOS[, 4, , ] <- M_CTC[4]
Mjuv_HOS[, 5, , ] <- M_CTC[5]

p_mature_RS <- array(0, c(nsim, maxage, nyears + proyears, n_r)) # Traditionals mature earlier
p_mature_RS[] <- array(matt_avg[, , sim_samp], c(maxage, n_r, nsim, nyears + proyears)) %>%
  aperm(c(3, 1, 4, 2))

# Sarita releases
sarita_rel <- readxl::read_excel(
  "data/Sarita/2025-05-08_Sarita_ReleaseRep_2000-2024.xlsx",
  sheet = "Actual Release"
)

# 2023 releases
sarita_rel %>%
  summarise(marked = sum(), n = sum(TotalRelease), .by = c(RELEASE_STAGE_NAME, RELEASE_YEAR)) %>%
  filter(RELEASE_YEAR == 2023)

stray <- c(0, 0.11, 0.16, 0.69, 0.04) * 100 # Proportions based on Sarita CWT escapement of traditionals in 2018
h2 <- EnvStats::rnormTrunc(nsim, 0.25, 0.15, min = 0, max = 0.5)

#### Note brood rule and in-river harvest rule are defined in script 4
Hatchery <- new(
  "Hatchery",
  n_r = n_r,
  n_yearling = c(335000, 165000), # Sarita smalls and traditionals
  n_subyearling = c(0, 0),
  s_prespawn = 0.85,  # Hatchery data, Sarita AHA inputs (Sarita CN AHA inputs.xlsx)
  s_egg_smolt = 0.9,  # Assumed 10 percent mortality shortly after release
  s_egg_subyearling = 1,
  Mjuv_HOS = Mjuv_HOS,
  p_mature_HOS = p_mature_RS,
  #stray_external = matrix(c(rep(0, 5), stray), maxage, 2), # Turn off strays for now
  gamma = 0.8,  # HSRG standard, Sarita AHA inputs
  m = 1,
  #f_brood = f_brood,  # Function defined in script 4
  pmax_esc = 0.33,
  pmax_NOB = 0.5,      # SEP guideline, suggested by Lian
  ptarget_NOB = 0.50,  # Hatchery data, Sarita AHA inputs, will evaluate a grid of 50, 75, 100 percent (see f_brood function in script 4)
  phatchery = NA,
  premove_HOS = 0,   # Stand-in for in-river fishery with HOS exploitation rates of 0.5, 0.7, or 0.99 (see f_catch function in script 4)
  fec_brood = fec, #rep(3625, maxage) is used from Hatchery data, Sarita AHA input
  fitness_type = c("Ford", "none"),
  theta = c(100, 80),
  rel_loss = rep(0, 3),
  zbar_start = c(90, 80),
  fitness_variance = 100,
  phenotype_variance = 10,
  heritability = h2,
  fitness_floor = 0.01
)

### Historical object ----
# Specify the initial spawners from 2023 escapement of 3000 spawners (Res Doc)
# Conditioning model suggests hatchery dominant system, arbitrarily setting pHOS = 2/3
#Historical <- new(
#  "Historical",
#  HistSpawner_NOS = 1/3 * 3000,
#  HistSpawner_HOS = 2/3 * 3000
#)

# Calculate HistNjuv_NOS (which return year?)
nyears_cm <- 45 - 5 # Sample abundance from final year with complete brood life cycle
#Njuv_CM <- sapply(report_RBT[sim_samp], function(x) x$N[nyears_cm + 1, , ],
#                  simplify = "array")

# Assume 50-50 ratio of spawners by life history group and release strategy
NOS <- sapply(report_RBT[sim_samp], function(x) x$syear[seq(nyears_cm - nyears + 1, nyears_cm), , 1],
              simplify = "array") %>%
  aperm(c(3, 2, 1))
HOS <- sapply(report_RBT[sim_samp], function(x) x$syear[seq(nyears_cm - nyears + 1, nyears_cm), , 2],
              simplify = "array") %>%
  aperm(c(3, 2, 1))
#colSums(NOS)
#colSums(HOS)
#colSums(HOS)/(colSums(NOS) + colSums(HOS))
#colSums(NOS) + colSums(HOS)

NOS_g <- sapply(1:Bio@n_g, function(g) {
  0.5 * NOS
}, simplify = "array")

HOS_r <- sapply(1:Hatchery@n_r, function(r) {
  0.5 * HOS
}, simplify = "array")

# Get F
FPT <- sapply(report_RBT[sim_samp], function(x) x$FPT[seq(nyears_cm - nyears + 1, nyears_cm)])
FT <- sapply(report_RBT[sim_samp], function(x) x$FT[seq(nyears_cm - nyears + 1, nyears_cm)])

Njuv_NOS <- sapply(1:Bio@n_g, function(g) {
  N <- sapply(report_RBT[sim_samp], getElement, "N", simplify = "array")

  Njuv <- array(0, c(nsim, maxage, nyears + 1))
  Njuv[, 1, ] <- Bio@p_LHG[g] * t(N[seq(nyears_cm - nyears + 1, nyears_cm + 1), 1, 1, ]) # Age - 1
  Njuv[, -1, 1] <- 0.5 * t(N[nyears_cm - nyears + 1, -1, 1, ]) # Year 1

  for (y in seq(2, nyears + 1)) {
    surv <- exp(-t(vulPT) * FPT[y-1, ] - Bio@Mjuv_NOS[, , y-1, g])
    Njuv[, seq(2, maxage), y] <- Njuv[, seq(1, maxage - 1), y-1] * surv[, seq(1, maxage - 1)] *
      (1 - Bio@p_mature[, seq(1, maxage - 1), y-1])
  }
  return(Njuv)
}, simplify = 'array')

Njuv_HOS <- sapply(1:Hatchery@n_r, function(r) {
  N <- sapply(report_RBT[sim_samp], getElement, "N", simplify = "array")

  Njuv <- array(0, c(nsim, maxage, nyears + 1))
  Njuv[, 1, ] <- 0.5 * t(N[seq(nyears_cm - nyears + 1, nyears_cm + 1), 1, 2, ]) # Age - 1
  Njuv[, -1, 1] <- 0.5 * t(N[nyears_cm - nyears + 1, -1, 2, ]) # Year 1

  for (y in seq(2, nyears + 1)) {
    surv <- exp(-t(vulPT) * FPT[y-1, ] - Hatchery@Mjuv_HOS[, , y-1, r])
    Njuv[, seq(2, maxage), y] <- Njuv[, seq(1, maxage - 1), y-1] * surv[, seq(1, maxage - 1)] *
      (1 - Hatchery@p_mature_HOS[, seq(1, maxage - 1), y-1, r])
  }
  return(Njuv)
}, simplify = 'array')

Historical <- new(
  "Historical",
  HistSpawner_NOS = NOS_g,
  HistSpawner_HOS = HOS_r,
  HistNjuv_NOS = Njuv_NOS,
  HistNjuv_HOS = Njuv_HOS,
  HistFPT = array(t(FPT), c(nsim, nyears, 2)),
  HistFT = array(t(FT), c(nsim, nyears, 2))
)




### Habitat ----
fry_surv <- read.csv("data/Sarita/fry_surv.csv")
fry_surv_year <- read.csv("data/Sarita/fry_surv_year.csv")

#### Get variation in historical freshwater survival
fry_surv_year <- fry_surv_year %>%
  mutate(fpe = y/3900/0.4)

fry_surv_year %>%
  summarise(m = mean(fpe), sd = sd(fpe))

# Log-linear regression to reproduce Figure 3.7 of Brown et al. Res Doc
if (FALSE) {
  fit <- lm(log(y) ~ x, fry_surv_year)
  xpred <- seq(0, 100)
  ypred <- predict(fit, newdata = data.frame(x = xpred), se.fit = TRUE)

  #g <- data.frame(
  #  xpred = xpred,
  #  ypred = ypred
  #) %>%
  #  mutate(lower = ypred.fit - 1.96 * ypred.se.fit,
  #         upper = ypred.fit + 1.96 * ypred.se.fit) %>%
  #  ggplot(aes(xpred, exp(ypred.fit))) +
  #  geom_ribbon(aes(ymin = exp(lower), ymax = exp(upper)), colour = NA, fill = "grey80") +
  #  geom_line() +
  #  geom_point(data = fry_surv_year, aes(x, y)) +
  #  geom_label(data = fry_surv_year, aes(x, y, label = year), vjust = "bottom", hjust = "left") +
  #  labs(x = expression("Percent Days in October with Flow > 5"~~(m^3/sec)),
  #       y = "Fry Per Returning Adult") +
  #  coord_cartesian(xlim = c(0, 105))

  g <- ggplot(fry_surv, aes(x, median)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), colour = NA, fill = "grey80") +
    geom_line() +
    geom_point(data = fry_surv_year, aes(x, y)) +
    geom_label(data = fry_surv_year, aes(x, y, label = year),
               size = 0.75 * GeomLabel$default_aes$size,
               vjust = "bottom", hjust = "left") +
    labs(x = expression("Percent Days in October with Flow > 5"~~(m^3/sec)),
         y = "Fry Per Returning Adult") +
    coord_cartesian(xlim = c(10, 105))
  ggsave("figures/Sarita-fry-survival.png", g, height = 4, width = 6)

  g <- fry_surv_year %>%
    select(year, y, fpe) %>%
    rename(`Fry/Return` = y, Survival = fpe) %>%
    reshape2::melt(id.vars = "year") %>%
    ggplot(aes(year, value)) +
    geom_point() +
    geom_line() +
    expand_limits(y = 0) +
    facet_wrap(vars(variable), ncol = 1, scales = "free_y", strip.position = "left") +
    labs(x = "Year", y = NULL) +
    theme(strip.placement = "outside", strip.background = element_blank())
  ggsave("figures/Sarita-fry-survival-ts.png", g, height = 4, width = 4)
}

#### Code used to stochastically simulate freshwater survival from historical October flow
#### Not used as of January 5, 2026
#### Average conditions 2017-2023 (x = environmental variable, y = fry/spawner)
#mean(fry_surv_year$x)
#get_eggfry_surv <- function(env_series, seed = 342, avg_fec = 3900, p_female = 0.4, nsim = 100) {
#  set.seed(seed)
#
#  fps_sim <- lapply(1:nsim, function(x) {
#    data.frame(
#      year = seq(1, proyears),
#      x = env_series
#    ) %>%
#      mutate(Simulation = .env$x)
#  }) %>%
#    bind_rows() %>%
#    left_join(fry_surv) %>%
#    mutate(fps = rlnorm(nrow(.), log(median), sd_lower)) %>%
#    mutate(fpe = fps/avg_fec/p_female)
#
#  return(fps_sim)
#
#}
#env_series <- rep(60, proyears)
#sim_surv <- get_eggfry_surv(env_series, nsim = nsim)
#fps <- reshape2::acast(sim_surv, list("Simulation", "year"), value.var = "fps")
#fpe <- reshape2::acast(sim_surv, list("Simulation", "year"), value.var = "fpe")
#matplot(t(fps), typ = 'l')
#matplot(t(fpe), typ = 'l')

# CV = 0.79
fry_surv_sd <- fry_surv_year %>%
  summarise(m = mean(qlogis(fpe)), sd = sd(qlogis(fpe))) %>%
  pull(sd) %>%
  round(2)

set.seed(234)
fry_surv_sim <- lapply(c(0.1, 0.2, 0.3), function(m) {
  samp <- rnorm(nsim * proyears, qlogis(m), fry_surv_sd) %>% matrix(nsim, proyears)
  samp_trans <- plogis(samp)
  samp_trans/mean(samp_trans) * m
})

if (FALSE) {

  fry_surv_sim %>% lapply(mean)
  fry_surv_sim %>% lapply(sd)

  png("figures/Sarita-fry-survival-sim.png", height = 6, width = 5, units = "in", res = 400)
  par(mfrow = c(3, 1), mar = c(2, 3, 3, 1), oma = c(3, 1.5, 0, 0))
  matplot(fry_surv_sim[[1]][1:3, ] %>% t(),
          main = "10% mean survival",
          type = "l", lty = 1, panel.first = grid(), ylim = c(0, 1), ylab = "")
  matplot(fry_surv_sim[[2]][1:3, ] %>% t(),
          main = "20% mean survival",
          type = "l", lty = 1, panel.first = grid(), ylim = c(0, 1), ylab = "")
  matplot(fry_surv_sim[[3]][1:3, ] %>% t(),
          main = "30% mean survival",
          type = "l", lty = 1, panel.first = grid(), ylim = c(0, 1), ylab = "")
  mtext("Projection year", side = 1, outer = TRUE, line = 1)
  mtext("Freshwater survival", side = 2, outer = TRUE, line = 0)
  dev.off()

}



Habitat <- new(
  "Habitat",
  use_habitat = TRUE,
  fry_rel = "BH",
  fry_prod = 1,
  fry_capacity = Inf,
  fry_sdev = fry_surv_sim[[1]]
)


SOM_low <- new("SOM",
               Name = "Sarita base, 2 LHG, 2 RS",
               nsim = nsim,
               nyears = nyears,
               proyears = proyears,
               seed = 1,
               Bio = Bio,
               Habitat = Habitat,
               Hatchery = Hatchery,
               Harvest = Harvest,
               Historical = Historical)
saveRDS(SOM_low, "SOM/SOM_lowsurv.rds")

SOM_medium <- SOM_low
SOM_medium@Habitat@fry_sdev <- fry_surv_sim[[2]]
saveRDS(SOM_medium, "SOM/SOM_mediumsurv.rds")

SOM_high <- SOM_low
SOM_high@Habitat@fry_sdev <- fry_surv_sim[[3]]
saveRDS(SOM_high, "SOM/SOM_highsurv.rds")




# Scenario with declining maturity and fecundity
# Get posterior medians
# Do regression vs time (since 1990-2021)
# Apply slope through projection

#g <- salmonMSE:::CM_maturity(report_RBT, salmonMSE:::get_CMdata(ERM_Sarita@.MISC$CMfit), year1 = 1979, brood = FALSE)

matt_med <- sapply(report_RBT, getElement, "matt", simplify = "array") %>%
  qlogis() %>%
  apply(c(1:2), median)
#c(1990, 2018) - 1979 + 1
#matplot(matt_med[seq(12, 40), ], typ = 'l')

matt_slope <- matt_med[seq(12, 40), ] %>%
  reshape2::melt() %>%
  rename(Year = Var1, Age = Var2) %>%
  filter(is.finite(value)) %>%
  summarise(slope = lm(value ~ Year) %>% coef() %>% getElement(2), .by = Age)

SOM_decfec <- SOM_medium
for (i in matt_slope$Age) {
  matt_i <- matrix(qlogis(SOM_decfec@Bio@p_mature[, i, SOM_decfec@nyears]), SOM_decfec@nsim, SOM_decfec@proyears)
  for (y in 2:SOM_decfec@proyears) {
    matt_i[, y] <- matt_i[, 1] + matt_slope$slope[matt_slope$Age == i] * y
  }
  SOM_decfec@Bio@p_mature[, i, SOM_decfec@nyears + seq(1, SOM_decfec@proyears)] <- plogis(matt_i)

  for (r in SOM_decfec@Hatchery@n_r) {
    matt_i <- matrix(qlogis(SOM_decfec@Hatchery@p_mature_HOS[, i, SOM_decfec@nyears, r]), SOM_decfec@nsim, SOM_decfec@proyears)
    for (y in 2:SOM_decfec@proyears) {
      matt_i[, y] <- matt_i[, 1] + matt_slope$slope[matt_slope$Age == i] * y
    }
    SOM_decfec@Hatchery@p_mature_HOS[, i, SOM_decfec@nyears + seq(1, SOM_decfec@proyears), r] <- plogis(matt_i)
  }

}

fec_val <- read.csv("tables/fec_decline.csv")
fec_decline <- array(fec_Sarita * p_female, c(maxage, nsim, nyears + proyears)) %>% aperm(c(2, 1, 3))
for (a in 3:5) {
  fec_decline[, a, nyears + seq(1, proyears)] <- matrix(
    as.numeric(fec_val[a-2, -1]) * fec_Sarita[a]/fec_val[a-2, 2] * p_female[a],
    nsim,
    proyears,
    byrow = TRUE
  )
}
SOM_decfec@Bio@fec <- fec_decline
SOM_decfec@Hatchery@fec_brood <- fec_decline
saveRDS(SOM_decfec, "SOM/SOM_declinematfec.rds")


if (FALSE) {

  #### Base maturity and vulnerability ----
  png("figures/SMSE/maturity_vul.png", height = 6, width = 3, units = "in", res = 400)
  par(mfrow = c(3, 1), mar = c(5, 4, 1, 1))

  salmonMSE:::plot_Mjuv_RS(SOM_low@Hatchery@p_mature_HOS[, , 1, ],
                           RS_names = c("Fed Fry", "Traditionals"), ylab = "Proportion mature")

  salmonMSE:::plot_SOM(SOM_low@Harvest, "vulPT",
                       type = "age", nsim = SOM_low@nsim, maxage = SOM_low@Bio@maxage,
                       nyears = SOM_low@nyears, proyears = SOM_low@proyears,
                       ylab = "Juvenile fishery vulnerability")

  salmonMSE:::plot_SOM(SOM@Harvest, "vulT",
                       type = "age", nsim = SOM_low@nsim, maxage = SOM_low@Bio@maxage,
                       nyears = SOM_low@nyears, proyears = SOM_low@proyears,
                       ylab = "Terminal fishery vulnerability")

  dev.off()

  #### Time-varying maturity and fecundity ----
  matt_new <- lapply(1:length(sim_samp), function(x) {
    out <- report_RBT[[sim_samp[x]]]["matt"]
    matt_proj <- SOM_low@Hatchery@p_mature_HOS[x, , SOM_low@nyears + seq(1, SOM_low@proyears), 2] %>%
      t() %>%
      array(c(SOM_low@proyears, maxage, 1))
    out$matt <- abind::abind(out$matt, matt_proj, along = 1)
    return(out)
  })
  g1 <- salmonMSE:::CM_maturity(matt_new, salmonMSE:::get_CMdata(ERM_Sarita@.MISC$CMfit), year1 = 1979, brood = TRUE) +
    ggtitle("Constant maturity (2013-2018 average)") +
    labs(y = "Traditionals maturity")


  matt_new <- lapply(1:length(sim_samp), function(x) {
    out <- report_RBT[[sim_samp[x]]]["matt"]
    matt_proj <- SOM_decfec@Hatchery@p_mature_HOS[x, , SOM_decfec@nyears + seq(1, SOM_decfec@proyears), 2] %>%
      t() %>%
      array(c(SOM_decfec@proyears, maxage, 1))
    out$matt <- abind::abind(out$matt, matt_proj, along = 1)
    return(out)
  })
  g2 <- salmonMSE:::CM_maturity(matt_new, salmonMSE:::get_CMdata(ERM_Sarita@.MISC$CMfit), year1 = 1979, brood = TRUE) +
    ggtitle("Declining maturity (1990-2018 trend)") +
    labs(y = "Traditionals maturity")

  g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "right")
  ggsave("figures/SMSE/Sarita_proj_maturity.png", g, height = 6, width = 6)

  g <- (SOM_decfec@Bio@fec[1, , ]/p_female) %>%
    reshape2::melt() %>%
    rename(Age = Var1, Year = Var2) %>%
    filter(Age > 1) %>%
    ggplot(aes(Year, value, colour = factor(Age))) +
    geom_line() +
    geom_point() +
    labs(y = "Fecundity", colour = "Age")
  ggsave("figures/SMSE/Sarita_decline_fecundity.png", g, height = 3, width = 5)
}
