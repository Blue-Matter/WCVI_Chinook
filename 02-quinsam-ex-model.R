
library(tidyverse)
library(readxl)

# Quinsam CWT recovery
rec <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Estimated"
) %>%
  mutate(is_catch = TotCatch > 0, is_esc = Escape > 0)

#table(rec$is_catch, rec$is_esc)
#filter(rec, is_catch & is_esc) %>% View()

# Catch area names
#cnames <- colnames(rec)
#data.frame(
#  Name = cnames[grepl("1|2|3", cnames)]
#) %>%
#  readr::write_csv(file = "data/Quinsam/Quinsam_fisheries.csv")


# CWT recoveries
#cwt <- summarise(rec,
#                 n_catch = sum(TotCatch),
#                 n_esc = sum(Escape),
#                 .by = c(Age, BROOD_YEAR))
#
#g <- ggplot(cwt, aes(BROOD_YEAR, n_catch)) +
#  facet_wrap(vars(Age), scales = "free_y") +
#  geom_line() +
#  geom_point() +
#  expand_limits(y = 0)
#g
#
#g <- ggplot(cwt, aes(BROOD_YEAR, n_esc)) +
#  facet_wrap(vars(Age), scales = "free_y") +
#  geom_line() +
#  geom_point() +
#  expand_limits(y = 0)
#g

# CWT by release strategy
cwt_rs <- rec %>%
  mutate(RS = ifelse(RELEASE_STAGE_NAME == "Fed Fry", "Fed Fry", "Seapen/Smolt 0+")) %>%
  summarise(
    n_catch = sum(TotCatch),
    n_esc = sum(Escape),
    .by = c(Age, RELEASE_YEAR, RS)
  )

g <- cwt_rs %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age), BROOD_YEAR = RELEASE_YEAR - 1) %>%
  ggplot(aes(BROOD_YEAR, n_catch, colour = RS)) +
  facet_grid(vars(RS), vars(Age), scales = "free_y") +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT catch") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_catch.png", g, height = 3.5, width = 6)

g <- cwt_rs %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age), BROOD_YEAR = RELEASE_YEAR - 1) %>%
  ggplot(aes(BROOD_YEAR, n_esc)) +
  facet_grid(vars(RS), vars(Age), scales = "free_y") +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT escapement", colour = "Release Strategy") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_esc.png", g, height = 3.5, width = 6)


# Quinsam - CWT releases
rel <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Releases"
)

rel_rs <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Fed Fry", "Smolt 0+", "Seapen 0+")) %>%
  mutate(RS = ifelse(RELEASE_STAGE_NAME == "Fed Fry", "Fed Fry", "Seapen/Smolt 0+")) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(RELEASE_YEAR, RS))

g <- rel_rs %>%
  mutate(BROOD_YEAR = RELEASE_YEAR - 1) %>%
  ggplot(aes(BROOD_YEAR, n_CWT)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(RS), ncol = 1, scales = "free_y") +
  labs(x = "Brood year", y = "CWT releases") +
  theme(legend.position = "bottom") +
  expand_limits(y = 0)
ggsave("figures/Quinsam_CWT_rel.png", g, height = 5, width = 4)


# Fit model ----
full_table <- expand.grid(
  RELEASE_YEAR = seq(min(cwt_rs$RELEASE_YEAR), 2024), # 2006 - 2023
  Age = seq(1, 5),
  RS = c("Fed Fry", "Seapen/Smolt 0+")
) %>%
  left_join(cwt_rs)

cwt_catch <- reshape2::acast(full_table, list("RELEASE_YEAR", "Age", "RS"), value.var = "n_catch", fill = 0)
cwt_esc <- reshape2::acast(full_table, list("RELEASE_YEAR", "Age", "RS"), value.var = "n_esc", fill = 0)

cwt_rel <- left_join(
  full_table %>% filter(Age == 1) %>% select(RELEASE_YEAR, RS),
  rel_rs
) %>%
  reshape2::acast(list("RELEASE_YEAR", "RS"), fill = 0)

# Sarita esc
esc <- readr::read_csv("data/R-OUT_infilled_indicators_escapement_timeseries.csv") %>%
  filter(river == "sarita_river") %>%
  arrange(year) %>%
  right_join(
    full_table %>% filter(Age == 1, RS == "Fed Fry") %>% select(RELEASE_YEAR),
    by = c("year" = "RELEASE_YEAR")
  )


# Data object for model
Ldyr <- dim(cwt_esc)[1]
Nages <- 5

mat <- c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1)
vulT <- rep(0, Nages)

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1))

#fec_Cowichan <- c(0, 87, 1153, 2780, 2700)
fec_Sarita <- c(0, 3000, 3000, 3600, 4600) # See Res Doc
d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 2,
  cwtrelease = cwt_rel,
  cwtesc = round(cwt_esc),
  cwtcatPT = round(cwt_catch),
  cwtcatT = NULL,
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  mobase = M_CTC,
  bmatt = mat,
  hatchsurv = 1,
  gamma = 0.8,
  ssum = 0.4,
  fec = fec_Sarita,
  r_matt = 2,
  obsescape = esc$escapement,
  propwildspawn = rep(1, Ldyr),
  hatchrelease = rep(0, Ldyr + 1),
  finitPT = 0.8,
  finitT = 0.8,
  cwtExp = 1
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
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 2,
                  control=list(adapt_delta = 0.95,
                               #stepsize = 0.01,
                               max_treedepth = 15))
saveRDS(samp, file = "CM/Quinsam_CM_01.05.26.rds")

samp <- readRDS(file = "CM/Quinsam_CM_01.05.26.rds")
#report <- salmonMSE:::get_report(samp)
#d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#shinystan::launch_shinystan(samp)

rs_names <- c("Fed Fry", "Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Sarita (Quinsam CWT)", year = unique(full_table$RELEASE_YEAR),
  dir = "CM", filename = "Quinsam_01.05"
)
