
library(salmonMSE)
library(ggpubr)

# These plots are all available in the markdown report
p_esc <- 0.75

# Get MCMC values
samp <- readRDS(paste0("CM/Sarita_RBT_CM_01.05.26_pesc", p_esc, ".rds"))
report <- get_report(samp)

# Get data
d <- get_CMdata(samp@.MISC$CMfit)

# Plot maturity
g <- CM_maturity(report, d, 1979, brood = FALSE, annual = TRUE) +
  ggtitle("Calendar year maturity")
ggsave("figures/Sarita_maturity_calendar.png", g, height = 8, width = 6)

g <- CM_maturity(report, d, 1979, brood = TRUE, annual = TRUE) +
  ggtitle("Brood year maturity")
ggsave("figures/Sarita_maturity_brood.png", g, height = 8, width = 6)

g <- CM_maturity(report, d, 1979, brood = TRUE, annual = FALSE)
ggsave("figures/Sarita_maturity_summary.png", g, height = 3, width = 5)

g <- salmonMSE:::CM_ts_origin(report, 1979, var = "Spawners", xlab = "Return Year")
ggsave("figures/Sarita_spawners.png", g, height = 3, width = 5)


# Plot F
g1 <- CM_F(report, year1 = 1979)
g2 <- CM_F(report, PT = FALSE, year1 = 1979)
g <- ggpubr::ggarrange(g1, g2, ncol = 1)
ggsave("figures/Sarita_F.png", g, height = 6, width = 4)

g1$data %>% filter(Year %in% 2013:2018) %>% mutate(u = 1 - exp(-`50%`)) %>% pull(u) %>% mean() # Preterminal exploitation rate
g2$data %>% filter(Year %in% 2013:2018) %>% mutate(u = 1 - exp(-`50%`)) %>% pull(u) %>% mean() # Terminal exploitation rate


# Exploitation rate
g3 <- CM_ER(report, brood = FALSE, type = "PT", year1 = 1979, at_age = FALSE)
g4 <- CM_ER(report, brood = FALSE, type = "T", year1 = 1979, at_age = FALSE)

g5 <- CM_ER(report, brood = TRUE, type = "PT", year1 = 1979, at_age = FALSE)
g6 <- CM_ER(report, brood = TRUE, type = "T", year1 = 1979, at_age = FALSE)

g7 <- CM_ER(report, brood = FALSE, type = "all", year1 = 1979, at_age = FALSE)
g8 <- CM_ER(report, brood = TRUE, type = "all", year1 = 1979, at_age = FALSE)

#g <- ggpubr::ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 2, nrow = 4,
#                       labels = paste0("(", 1:8, ")"))
#ggsave("figures/Sarita_ex.png", g, height = 9, width = 6)

g <- ggpubr::ggarrange(g3, g5, g4, g6, g7, g8, ncol = 2, nrow = 3,
                       labels = paste0("(", LETTERS[1:6], ")"))
ggsave("figures/Sarita_ex.png", g, height = 7, width = 6)

g <- ggpubr::ggarrange(g3, g4, ncol = 1, nrow = 2)
ggsave("figures/Sarita_CYER.png", g, height = 5, width = 4)

# Preterminal CYER
ER_summary <- full_join(
  g3$data %>% select(Year, `50%`) %>% rename(`Preterminal CYER` = `50%`),
  g4$data %>% select(Year, `50%`) %>% rename(`Terminal CYER` = `50%`)
) %>%
  round(2) %>%
  filter(Year %in% c(2013:2018))
mean(ER_summary$`Preterminal CYER`)
mean(ER_summary$`Terminal CYER`)


# Compare exploitation rates to CTC model
ctc <- readr::read_csv("data/rbt_er.csv")[, -1] %>%
  rename(Year = brood) %>%
  mutate(Age = paste("Age", age))

g <- salmonMSE:::CM_ER(report, type = "PT", year1 = 1979, at_age = TRUE) +
  facet_wrap(vars(Age), scales = "free_y") +
  geom_line(data = ctc, aes(y = preterminal_er), colour = "red") +
  guides(colour = guide_legend(title = "Model"))
g$data %>% filter(Year %in% c(2013:2018), Age == "Age 5") %>% pull(`50%`) %>% mean() # 0.43
ggsave(paste0("figures/Sarita_ex_PT_CTC_pesc", p_esc, ".png"), g, height = 4, width = 6)

g <- salmonMSE:::CM_ER(report, type = "T", year1 = 1979, at_age = TRUE) +
  facet_wrap(vars(Age), scales = "free_y") +
  geom_line(data = ctc, aes(y = terminal_er), colour = "red")
g$data %>% filter(Year %in% c(2013:2018), Age == "Age 5") %>% pull(`50%`) %>% mean() # 0.24
ggsave(paste0("figures/Sarita_ex_T_CTC_pesc", p_esc, ".png"), g, height = 4, width = 6)

#g <- salmonMSE:::CM_ER(report, brood = FALSE, type = "all", year1 = 1979, at_age = FALSE)












report_MPD <- samp@.MISC$CMfit$obj$report()
plot(1:d$Nages, report_MPD$vulPT, typ = "o")
plot(1:d$Nages, report_MPD$vulT, typ = "o")

matplot(report_MPD$mo, typ = "o")

full_year <- seq(1979, 2023)

SAMtool::plot_composition(
  full_year,
  obs = d$cwtesc[, , 1],
  fit = report_MPD$ebrood[, , 1],
  annual_yscale = 'raw'
)

SAMtool::plot_composition(
  full_year,
  obs = d$cwtcatPT[, , 1],
  fit = report_MPD$cbroodPT[, , 1],
  annual_yscale = 'raw'
)

SAMtool::plot_composition(
  full_year,
  obs = report_MPD$matt[, , 1],
  annual_yscale = 'raw'
)

SAMtool::plot_composition(
  full_year,
  obs = d$cwtcatT[, , 1],
  fit = report_MPD$cbroodT[, , 1],
  annual_yscale = 'raw'
)

par(mar = c(5, 4, 1, 1), oma = rep(0, 4), mfrow = c(1, 1))
plot(full_year, d$obsescape, typ = 'o')
lines(full_year, exp(report_MPD$logpredesc), col = 2, typ = 'o')

plot(full_year, report_MPD$egg, typ = 'o')
plot(full_year, report_MPD$N[1:d$Ldyr, 1, 1], typ = 'o', ylab = "Smolt production")

plot(full_year, report_MPD$FPT, typ = 'o')
plot(full_year, report_MPD$FT, typ = 'o')

# SRR
plot(report_MPD$egg, report_MPD$N[1:d$Ldyr + 1, 1, 1])

memin <- report_MPD$memin
memax <- report_MPD$memax
log_cr <- memax - memin
exp(log_cr)

salmonMSE:::CM_Megg(report, 1979, surv = FALSE)
salmonMSE:::CM_Megg(report, 1979, surv = TRUE)

