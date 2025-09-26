
library(tidyverse)
library(salmonMSE)


#### Data ----
# Sarita escapement (1979 - 2023)
esc <- readr::read_csv("data/R-OUT_infilled_indicators_escapement_timeseries.csv") %>%
  filter(river == "sarita_river") %>%
  arrange(year)

# LGL estimates
esc_lgl <- readr::read_csv("data/Sarita/SaritaCN_annualVAR.csv")

g <- esc_lgl %>%
  mutate(type = "LGL") %>%
  ggplot(aes(Year, MedianEscapement, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = perc15, ymax = perc85), alpha = 0.5) +
  geom_point() +
  geom_line() +
  geom_line(data = mutate(esc, type = "Res Doc"), aes(x = year, y = escapement)) +
  expand_limits(y = 0) +
  labs(y = "Sarita escapement", colour = NULL, fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/Sarita_escapement_comparison.png", g, height = 4, width = 6)

# Sarita hatchery releases and broodtake
#rel <- readxl::read_excel(
#  "data/Sarita/2025-02-04 Sarita Chinook Releases and Removals.xlsx",
#  sheet = "Sarita CN Releases",
#  range = "A10:E50"
#)
rel <- readr::read_csv("data/Sarita/sarita_rel.csv")
rel[is.na(rel)] <- 0
rel <- rel %>%
  rename(Year = `Row Labels`) %>%
  mutate(Total = `Fed Fry` + `Seapen 0+` + `Smolt 0+`, `Smolt 1+`)

#brood <- readxl::read_excel(
#  "data/Sarita/2025-02-04 Sarita Chinook Releases and Removals.xlsx",
#  sheet = "Broodstock Used",
#  range = "A5:D39"
#) %>%
#  rename(Year = `Row Labels`)

brood <- readr::read_csv("data/Sarita/sarita_brood.csv") %>%
  rename(Year = `Row Labels`)
brood$Broodtake <- rowSums(brood[, -1], na.rm = TRUE)

png("figures/Sarita_brood.png", units = "in", height = 4, width = 6, res = 400)
par(mar = c(5, 4, 1, 1))
plot(`Sum of 01Used For Broodstock` ~ Year, brood, typ = "o", ylab = "Brood", ylim = c(0, 600))
lines(Broodtake ~ Year, brood, typ = "o", ylab = "Broodstock", col = "red", pch = 16)
graphics::grid()
legend("topleft", c("Broodtake", "Brood + Holding Mortality + Other"), col = 1:2, pch = c(1, 16), lty = 1, bty = "n")
dev.off()

#g <- sarita_rel %>% reshape2::melt(id.vars = "Year") %>%
#  ggplot(aes(Year, value, colour = variable)) +
#  geom_line() +
#  geom_point()

# CWT data from RBT
# broodyear 1973 - 2021
# runyear 1975-2023
cwt_dat <- readr::read_csv("data/RBT_data_wfisheries.csv")
#problems(dat)
cwt_dat[c(1013, 2273), ]

#### Process data ----

# Full matrix of ages (1-5) and years (1979 - 2023)
full_matrix <- expand.grid(
  BroodYear = 1979:2023,
  Age = 1:5
) %>%
  as.data.frame()
full_year <- data.frame(BroodYear = 1979:2023)

# Escapement  + broodtake
esc_sarita <- filter(esc, year %in% full_year$BroodYear) %>%
  left_join(brood %>% select(Year, Broodtake), by = c("year" = "Year")) %>%
  mutate(Broodtake = ifelse(is.na(Broodtake), 0, Broodtake)) %>%
  rename(spawners = escapement) %>%
  mutate(escapement = spawners + Broodtake,
         p_spawn = 1 - Broodtake/escapement)

# Sarita releases
rel_sarita <- left_join(full_year, rel, by = c("BroodYear" = "Year"))
rel_sarita$Total[is.na(rel_sarita$Total)] <- 0


# CWT releases by tag code
cwt_rel <- cwt_dat %>%
  summarise(CWTMark1Count = unique(CWTMark1Count), .by = c(TagCode, BroodYear)) %>%
  arrange(BroodYear, TagCode)

# Annual CWT releases
cwt_rel_annual <- cwt_rel %>%
  summarise(rel = sum(CWTMark1Count), .by = BroodYear) %>%
  right_join(full_year) %>%
  mutate(rel = ifelse(is.na(rel), 0, rel))

# CWT escapement by brood year, age
cwt_esc <- cwt_dat %>%
  filter(fishery_type == "escapement") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(BroodYear, Age)) %>%
  right_join(full_matrix, by = c("BroodYear", "Age")) %>%
  reshape2::acast(list("BroodYear", "Age"), value.var = "n", fill = 0)

# Plot CWT escapement
g <- cwt_dat %>%
  filter(fishery_type == "escapement") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(BroodYear, Age)) %>%
  right_join(full_matrix, by = c("BroodYear", "Age")) %>%
  arrange(Age, BroodYear) %>%
  mutate(p = n/sum(n, na.rm = TRUE), .by = BroodYear) %>%
  filter(!is.na(p)) %>%
  ggplot(aes(BroodYear, p, fill = factor(Age, levels = 5:2))) +
  geom_col(width = 1, colour = "grey40") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Brood Year", y = "Proportion", fill = "Age", title = "CWT escapement") +
  coord_cartesian(expand = FALSE)
ggsave("figures/CWT_esc_prop.png", g, height = 4, width = 6)

g2 <- g +
  coord_cartesian(expand = FALSE, xlim = c(2013.5, 2020.5))
ggsave("figures/CWT_esc_prop2.png", g2, height = 4, width = 6)

# Preterminal CWT
cwt_pt <- cwt_dat %>%
  filter(fishery_type == "pre-terminal", Age < 7) %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(BroodYear, Age)) %>%
  right_join(full_matrix, by = c("BroodYear", "Age")) %>%
  reshape2::acast(list("BroodYear", "Age"), value.var = "n", fill = 0)

# Terminal CWT
cwt_t <- cwt_dat %>%
  filter(fishery_type == "terminal") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(BroodYear, Age)) %>%
  right_join(full_matrix, by = c("BroodYear", "Age")) %>%
  reshape2::acast(list("BroodYear", "Age"), value.var = "n", fill = 0)



# Data object for model
Ldyr <- nrow(cwt_esc)
Nages <- 5

mat <- c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1)
vulT <- vulPT

# CTC 23-06 p. 8
surv <- c(0.1, 0.7, 0.8, 0.9, 0.9)
M_CTC <- -log(surv)

# See Brown et al. Res Doc p. 25 (section 3.3.3), use average of 3900 when needed
fec_Sarita <- c(0, 1500, 3000, 3600, 4600)

# Expansion factors for CWT
cwt_dat %>%
  mutate(EstimatedNumber = as.numeric(EstimatedNumber),
         Exp = AdjustedEstimatedNumber/EstimatedNumber) %>%
  select(EstimatedNumber, Exp, ExpansionFactor) %>%
  filter(Exp < 2) %>%
  pull(Exp) %>%
  hist(breaks = 50)


# Model assumption of catch expansion factor
# Use alternative values to change data weighting of CWT (re-adjust numbers accordingly)
cwtExp <- 1

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  cwtrelease = cwt_rel_annual$rel,
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatPT = array(round(cwt_pt/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatT = array(round(cwt_t/cwtExp), c(Ldyr, Nages, 1)),
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  bmatt = mat,
  mobase = M_CTC,
  hatchsurv = 0.9,
  gamma = 0.8,
  ssum = 0.4,
  fec = fec_Sarita,
  obsescape = esc_sarita$spawners,
  propwildspawn = round(esc_sarita$p_spawn, 2),
  hatchrelease = c(rel_sarita$Total, rel_sarita$Total[length(rel_sarita$Total)]),
  s_enroute = 1,
  finitPT = 0.8,
  finitT = 0.8,
  cwtExp = cwtExp
)


# Fix these parameters
map <- list()

# Fix maturity
#map$sd_matt <- factor(rep(NA, Nages-2)) # Not estimating year-specific maturity
#map$logit_matt <- factor(rep(NA, Ldyr * (Nages - 2)))

# Fix additional age-1 M
#map$moadd <- factor(NA)

# Fix age-1 density-independent M deviates
#map$wto <- factor(rep(NA, Ldyr))
#map$wto_sd <- factor(NA)

# Fix density dependent egg-smolt M deviates
#map$wt <- factor(rep(NA, Ldyr))
#map$wt_sd <- factor(NA)

# Fix observation error of Sarita escapement (needed, otherwise model can't separate process from obs error)
map$lnE_sd <- factor(NA)

start <- list(log_so = log(3 * max(d$obsescape)))

# Fit with sampling rate = 1
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE)
#dat <- salmonMSE:::get_CMdata(fit)
samp <- sample_CM(fit, chains = 3, cores = 3, iter = 10000, thin = 5)
saveRDS(samp, file = "CM/Sarita_RBT_CM_09.08.25.rds")

samp <- readRDS("CM/Sarita_RBT_CM_09.08.25.rds")
salmonMSE::report_CM(samp, dir = "CM", filename = "Sarita_09.08", year = full_year$BroodYear, name = "Sarita (RBT CWT)")

#shinystan::launch_shinystan(samp)
