

# Calculate decline in fecundity based on 5% reduction in mean size at age per decade

# Plot data
d <- readxl::read_excel(
  file.path("data", "2023-11-07 R-OUT_Fecundity-Data_all-years.xlsx")
) %>%
  filter(!is.na(`Resolved Age GR`)) %>%
  mutate(Age = substr(`Resolved Age GR`, 1, 1) |> as.integer(),
         Year = as.integer(Year))

g <- ggplot(d, aes(factor(Age), `Est Fecundity`)) +
  geom_point() +
  geom_violin() +
  facet_wrap(vars(Site))

g <- ggplot(d, aes(factor(Age), `Est Fecundity`)) +
  geom_point() +
  geom_violin() +
  facet_grid(vars(Year), vars(Site))

summarise(d, n = n()/3, m = mean(`Est Fecundity`), .by = c(Age, Site)) %>%
  arrange(Site, Age)


# Fecundity at length model
fit <- readRDS("data/WCVI_hatchery_Chinook_fecundity_predictive_model.rds")
summary(fit$`fecundity model`)

# Prediction grid
g <- fit$`prediction data example` %>%
  ggplot(aes(resolved_age_gr, poh_length_mm, fill = fit)) +
  geom_tile(height = 3.83, width = 1) +
  facet_wrap(vars(site_year)) +
  scale_fill_viridis_c()

fit$`fecundity model`$model

g <- fit$`prediction data example` %>%
  ggplot(aes(poh_center, poh_length_mm)) +
  geom_point() +
  facet_wrap(vars(site_year))


#poh_center = poh_length_mm - 667.425

# Get average length at age for RCH_2023
d <- fit$`fecundity model`$model %>%
  filter(site_year == "RCH_2021") %>%
  mutate(poh_length_mm = poh_center + 667.425)

g <- ggplot(d, aes(resolved_age_gr, poh_length_mm)) +
  geom_jitter()

len <- summarise(d, mlen = mean(poh_length_mm), .by = resolved_age_gr)

y <- 1:30

slope <- data.frame(
  resolved_age_gr = c('31', '41', '51'),
  slope = c(-4.49, -4.30, -5.35)/100*0.1
)


fec_decline <- lapply(y, function(yy) {
  newdata <- len %>%
    mutate(year = yy) %>%
    left_join(slope, by = "resolved_age_gr") %>%
    mutate(poh_length_mm = mlen * (1 + slope)^year) %>%
    mutate(poh_center = matrix(poh_length_mm - 667.425, ncol = 1), site_year = "RCH_2021")

  newdata$fec <- predict(fit$`fecundity model`, newdata = newdata)
  newdata

}) %>%
  bind_rows() %>%
  reshape2::dcast(list("resolved_age_gr", "year"), value.var = "fec")
write_csv(fec_decline, file = "tables/fec_decline.csv")
matplot(fec_decline[, -1] %>% t(), type = 'l', ylim = c(0, 5000))
