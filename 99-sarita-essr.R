
library(tidyverse)

# Releases
rel <- readxl::read_excel(
  file.path("data", "Sarita", "2025-05-08_Sarita_ReleaseRep_2000-2024.xlsx"),
  sheet = "Actual Release"
)

# Release strategy numbers
rs <- summarise(
  rel,
  TaggedNum = sum(TaggedNum),
  TotalRelease = sum(TotalRelease),
  .by = c(BROOD_YEAR, RELEASE_STAGE_NAME)
) %>%
  reshape2::melt(id.var = c("BROOD_YEAR", "RELEASE_STAGE_NAME"))
g <- ggplot(rs, aes(BROOD_YEAR, value)) +
  geom_point() +
  geom_line() +
  facet_grid(vars(variable), vars(RELEASE_STAGE_NAME), scales = "free_y") +
  expand_limits(y = 0)

# Total Sarita ESSR
return <- readr::read_csv(
  file.path("data", "ESSR", "River_Returns_2025_04_28_Sarita.csv")
)

ret <- return %>%
  summarise(N4 = sum(`04Excess To Spawning Req (ESSR)`, na.rm = TRUE),
            N3 = sum(`03Given To First Nations (FSC)`, na.rm = TRUE),
            .by = `Recovery Year`) %>%
  mutate(ESSR_catch = N4 + N3)

plot(ESSR_catch ~ `Recovery Year`, ret, typ = 'o')


# Clipped fish - per Lian's email
clip <- tibble::tibble(
  `Recovery Year` = seq(2020, 2024),
  N_noclip = c(26, 5, 0, 38, 3)
) %>% left_join(ret, by = "Recovery Year") %>%
  mutate(N_clip = ESSR_catch - N_noclip)

# CWT sampling (2023 not available)
samp_2020 <- readr::read_csv(
  file.path("data", "ESSR", "Sarita_ESSR_CWT_Batch_BatchId_7414_2020.csv")
)
samp_2021 <- readr::read_csv(
  file.path("data", "ESSR", "Sarita_ESSR_CWT_Batch_BatchId_7641_2021.csv")
)
samp_2022 <- readr::read_csv(
  file.path("data", "ESSR", "Sarita_ESSR_CWT_Batch_BatchId_7951_2022a.csv")
)
samp_2022_hfn <- readr::read_csv(
  file.path("data", "ESSR", "Sarita_HFN_ESSR_CWT_Batch_BatchId_7887_2022b.csv")
)
samp_2024 <- readr::read_csv(
  file.path("data", "ESSR", "Sarita_ESSR_CWT_Batch_BatchId_2024.csv")
)

samp <- rbind(
  samp_2020,
  samp_2021,
  samp_2022,
  samp_2022_hfn,
  samp_2024
)

#unique(samp$Status)
#unique(samp$TagCondition)
#table(samp$Status, samp$TagCondition)



# ESSR summary with total catch, number of clipped fish, and CWTs
cwt <- samp %>%
  summarise(CWT_samp = n(),
            n_CWT = sum(TagCondition == "Normal"),
            .by = `RecoveryYear`)

essr <- ret %>%
  left_join(clip) %>%
  left_join(cwt, by = c("Recovery Year" = "RecoveryYear")) %>%
  filter(`Recovery Year` >= 2020)

#essr %>%
#  filter(`Recovery Year` >= 2020) %>%
#  readr::write_csv(file = "tables/Sarita_ESSR_CWT.csv")


# Merge CWT catch with release
rs_tagcode <- summarise(
  rel,
  TaggedNum = sum(TaggedNum),
  TotalRelease = sum(TotalRelease),
  .by = c(BROOD_YEAR, MRP_TAGCODE, RELEASE_STAGE_NAME)
)
samp_expand <- samp %>%
  select(LabelId, Status, TagCode, RecoveryYear) %>%
  filter(Status == "Tag read OK") %>%
  left_join(rs_tagcode, by = c("TagCode" = "MRP_TAGCODE")) %>%
  left_join(essr, by = c("RecoveryYear" = "Recovery Year")) %>%
  #rename(Catch_noclip = N_noclip, Catch_clip = N_clip,
  #       Nsamp_CWT = CWT_samp) %>%
  select(!N4 & !N3) %>%
  mutate(Age = RecoveryYear - BROOD_YEAR)
readr::write_csv(samp_expand, "data/ESSR/Sarita_ESSR_CWT_recoveries.csv")

# Nominal recoveries, no catch expansion
samp_summary <- samp_expand %>%
  summarise(n = n(), .by = c(BROOD_YEAR, Age, RELEASE_STAGE_NAME)) %>%
  arrange(RELEASE_STAGE_NAME, BROOD_YEAR, Age) %>%
  mutate(A = paste("Age", Age))

g <- samp_summary %>%
  #filter(!is.na(Age)) %>%
  ggplot(aes(BROOD_YEAR, n, colour = RELEASE_STAGE_NAME)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(A)) +
  labs(x = "Brood Year", y = "Nominal CWT recoveries (ESSR)") +
  theme(legend.position = "bottom")
ggsave("figures/Sarita_ESSR_CWT_nominal.png", g, height = 4, width = 6)


# Total releases
sarita_sep <- readxl::read_excel(
  "data/Sarita/WCVI Chinook life cycle master Nov 2024.xlsx",
  sheet = "Sarita-dataSEP",
  range = "A20:AG83"
)
