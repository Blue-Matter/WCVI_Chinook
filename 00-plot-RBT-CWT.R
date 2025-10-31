
#### Plot Robertson Creek CWT ----
library(tidyverse)

dat <- readr::read_csv("data/RBT_data_wfisheries.csv") %>%
  mutate(TagCode2 = as.numeric(TagCode))
#problems(dat)
dat[c(1013, 2273), ]

# CWT releases
rel <- dat %>%
  summarise(CWTMark1Count = unique(CWTMark1Count), .by = c(TagCode, BroodYear)) %>%
  arrange(BroodYear, TagCode)

#rel2 <- readr::read_csv("data/releases_indicator_RBT.csv") %>%
#  select(BroodYear, TagCode, CWTMark1Count) %>%
#  arrange(BroodYear, TagCode)

g <- rel %>%
  summarise(rel = sum(CWTMark1Count), .by = BroodYear) %>%
  ggplot(aes(BroodYear, rel)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  labs(x = "Brood Year", y = "CWT releases")
ggsave("figures/rel.png", g, height = 3, width = 4)

# Release strategies
cwt_rs <- readxl::read_excel('data/2025-03-26_Robertson_CWT_Data.xlsx')
g <- cwt_rs %>%
  summarise(Rel = sum(`Total Released`), .by = c(BROOD_YEAR, `Release Strategy`)) %>%
  mutate(p = Rel/sum(Rel), .by = BROOD_YEAR) %>%
  ggplot(aes(BROOD_YEAR, Rel, colour = `Release Strategy`)) +
  geom_point() +
  geom_line()


#### Identify sources of CWT
# Preterminal = US and Canada
# Terminal = Canada
table(dat$fishery_type, dat$country)
table(dat$fishery_type, dat$ERA_fishery_name)

CWT_summary <- dat %>%
  #filter(fishery_type == "terminal") %>%
  summarise(n_records = n(),
            Catch = sum(AdjustedEstimatedNumber) %>% round(),
            .by = c(fishery_type, country, ERA_fishery_name, Coarse_description)) %>%
  arrange(desc(fishery_type), desc(Catch))
readr::write_csv(CWT_summary, "tables/RBT_CWT_catch_summary.csv")

unique(dat$TagCode) %>% length()
range(dat$ExpansionFactor %>% as.numeric(), na.rm = TRUE)
dat$ExpansionFactor %>% as.numeric() %>% is.na() %>% mean()

# Escapement by run year
esc_run <- dat %>%
  filter(fishery_type == "escapement") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(RunYear, Age))

#g <- esc_run %>%
#  ggplot(aes(Age, n)) +
#  facet_wrap(vars(RunYear)) +
#  geom_point() +
#  geom_line() +
#  ggtitle("Escapement at age by run year")

g <- esc_run %>%
  ggplot(aes(RunYear, n)) +
  geom_col(width = 1, colour = "black", fill = "grey80", linewidth = 0.1) +
  facet_wrap(vars(Age)) +
  ggtitle("Escapement at age by run year")
ggsave("figures/esc.png", g, width = 6, height = 4)

esc_annual <- esc_run %>%
  summarise(n = sum(n), .by = RunYear)
g <- esc_annual %>%
  ggplot(aes(RunYear, n)) +
  geom_line() +
  geom_point() +
  labs(y = "CWT escapement")
ggsave("figures/esc_total.png", g, width = 4, height = 3)


# Escapement by brood year
esc_brood <- dat %>%
  filter(fishery_type == "escapement") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(BroodYear, Age))

g <- esc_brood %>%
  ggplot(aes(BroodYear, n)) +
  geom_col(width = 1, colour = "black", fill = "grey80") +
  facet_wrap(vars(Age)) +
  ggtitle("Escapement at age")

#g <- esc_brood %>%
#  ggplot(aes(Age, n)) +
#  facet_wrap(vars(BroodYear)) +
#  geom_point() +
#  geom_line() +
#  ggtitle("Escapement at age by brood year")

g <- esc_brood %>%
  ggplot(aes(BroodYear, n)) +
  geom_col(width = 1, colour = "black", fill = "grey80", linewidth = 0.1) +
  facet_wrap(vars(Age)) +
  ggtitle("CWT Escapement at age")

# Preterminal catch - by run year and country
pt_country <- dat %>%
  filter(fishery_type == "pre-terminal") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(RunYear, Age, country))

g <- pt_country %>%
  ggplot(aes(RunYear, n)) +
  geom_col(width = 1, colour = "black", aes(fill = country), linewidth = 0.1) +
  #facet_wrap(vars(Age)) +
  facet_wrap(vars(Age), scales = "free_y") +
  ggtitle("CWT preterminal catch at age") +
  labs(fill = NULL) +
  #coord_trans(y = "sqrt") +
  theme(legend.position = "bottom")
ggsave("figures/pt_fishery.png", g, width = 6, height = 5.5)


# Terminal catch
term <- dat %>%
  filter(fishery_type == "terminal") %>%
  summarise(n = sum(AdjustedEstimatedNumber), .by = c(RunYear, Age))

g <- term %>%
  ggplot(aes(RunYear, n)) +
  geom_col(width = 1, colour = "black", fill = "grey80", linewidth = 0.1) +
  facet_wrap(vars(Age), scales = "free_y") +
  ggtitle("CWT terminal catch at age")
ggsave("figures/term_fishery.png", g, width = 6, height = 4)
