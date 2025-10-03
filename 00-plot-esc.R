

library(tidyverse)

# Plot WCVI escapement ----
esc <- readr::read_csv("data/R-OUT_infilled_indicators_escapement_timeseries.csv")

g <- esc %>%
  mutate(river = gsub("_", " ", river) %>% stringr::str_to_title()) %>%
  mutate(`Infill?` = ifelse(is.na(infill), "No", "Yes")) %>%
  ggplot(aes(year, escapement)) +
  geom_line() +
  geom_point(aes(shape = `Infill?`)) +
  labs(x = "Year", y = "Total escapement") +
  scale_shape_manual(values = c(16, 1)) +
  facet_wrap(vars(river), scales = "free_y", ncol = 3) +
  expand_limits(y = 0) +
  theme(legend.position = 'bottom')
ggsave("figures/total_escapement.png", g, height = 9, width = 7)

# Compare LGL estimates of Sarita escapement ----
esc_lgl <- readr::read_csv("data/Sarita/SaritaCN_annualVAR.csv")
esc_sarita <- esc %>%
  filter(river == "sarita_river") %>%
  arrange(year)

g <- esc_lgl %>%
  mutate(type = "LGL") %>%
  ggplot(aes(Year, MedianEscapement, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = perc15, ymax = perc85), alpha = 0.1) +
  geom_point() +
  geom_line() +
  geom_line(data = mutate(esc_sarita, type = "Brown et al."), aes(x = year, y = escapement)) +
  expand_limits(y = 0) +
  labs(y = "Sarita escapement", colour = NULL, fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/Sarita_escapement_comparison.png", g, height = 4, width = 6)
