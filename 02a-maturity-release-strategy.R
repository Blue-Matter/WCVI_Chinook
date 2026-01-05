
library(salmonMSE)

samp <- readRDS(file = "CM/Quinsam_CM_01.05.26.rds")
report <- get_report(samp)
d <- get_CMdata(samp@.MISC$CMfit)

# Difference in maturity between release strategies
g1 <- CM_maturity(report, d, year1 = 2005, r = 1) +
  ggtitle("Quinsam Fed Fry")
g2 <- CM_maturity(report, d, year1 = 2005, r = 2) +
  ggtitle("Quinsam Smolt 0+")
g <- ggpubr::ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")
ggsave("figures/Quinsam_maturity_RS.png", width = 6, height = 3)

mat <- rbind(
  g1$data %>% mutate(RS = "Fed Fry"),
  g2$data %>% mutate(RS = "Smolt 0+")
) %>%
  select(!Var2) %>%
  rename(`BroodYear` = Year)
#write_csv(mat, file = "tables/Quinsam_maturity.csv")

# Fit to cwt catch
rs_names <- c("Fed Fry", "Smolt 0+")
g <- salmonMSE:::CM_fit_CWTcatch(report, d, year1 = 2005, rs_names = rs_names)
ggsave("figures/Quinsam_CWTfit_catch.png", g, height = 4, width = 6)

g <- salmonMSE:::CM_fit_CWTesc(report, d, year1 = 2005, rs_names = rs_names)
ggsave("figures/Quinsam_CWTfit_esc.png", g, height = 4, width = 6)

# Exploitation rate and releases
g1 <- salmonMSE:::CM_CWTrel(d$cwtrelease, year1 = 2005, rs_names = rs_names)
g2 <- salmonMSE:::CM_CWT_ER(report, type = "all", year1 = 2005, rs_names = rs_names)
g <- ggpubr::ggarrange(g1, g2, ncol = 1)
ggsave("figures/Quinsam_CWT_ER.png", g, height = 6, width = 4)


# Predict proportion of escapement
# 2017
g <- sapply(report, function(i) {
  apply(i$ebrood[13, , ], 2, function(x) x/sum(x))
}, simplify = "array") %>%
  apply(1:2, median) %>%
  reshape2::melt() %>%
  mutate(Age = factor(Var1, 5:1), `Release Strategy` = factor(rs_names[Var2], rs_names)) %>%
  ggplot(aes(`Release Strategy`, value, fill = Age)) +
  geom_col(colour = "grey20") +
  labs(y = "Proportion escapement")
ggsave("figures/Quinsam_CWT_esc_prop.png", g, height = 4, width = 4)

# Compare maturity from Quinsam and RBT
samp_RBT <- readRDS(file = "CM/Sarita_RBT_CM_01.05.26_pesc0.75.rds")
report_RBT <- get_report(samp_RBT)

g_RBT <- salmonMSE:::CM_maturity(report_RBT, d = salmonMSE:::get_CMdata(samp_RBT@.MISC$CMfit),
                                 year1 = 1979, rs_names = "RBT Smolt 0+", annual = TRUE) %>%
  getElement("data") %>%
  filter(Year >= 2005)

rs_names <- c("Quinsam Fed Fry", "Quinsam Smolt 0+")
col_pal <- scales::pal_hue()(2) %>% c("black") %>% structure(names = c(rs_names, "RBT Smolt 0+"))
g <- salmonMSE:::CM_maturity(report, d, year1 = 2005, rs_names = rs_names, annual = TRUE) +
  geom_point(data = g_RBT) +
  geom_line(data = g_RBT) +
  geom_ribbon(data = g_RBT, aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.25, colour = NA) +
  scale_colour_manual(values = col_pal) +
  scale_fill_manual(values = col_pal)
ggsave("figures/Quinsam_CWT_maturity.png", g, height = 5.5, width = 6)

# Apply difference in Quinsam maturity to RBT to model Sarita
#x <- report[[1]]
#y <- report_RBT[[1]]

calc_Sarita_fedfry_matt <- function(x, y, type = c("logit", "ratio")) { # x = report Quinsam, y = report RBT
  type <- match.arg(type)

  # Fed fry/traditional (< 1 means fed fry mature later)
  if (type == "logit") {
    matt_dev <- qlogis(x$matt[, , 1]/x$matt[, , 2])
  } else {
    matt_dev <- x$matt[, , 1]/x$matt[, , 2]
  }
  i <- seq(1979, 2023) %in% seq(2006, 2023)

  matt_new <- array(0, dim(x$matt)) # 1 = Fed fry, 2 = Smolt 0+
  matt_new[, , 2] <- y$matt[i, , 1] # Traditionals

  # Fed fry
  if (type == "logit") {
    matt_new[, , 1] <- plogis(qlogis(y$matt[i, , 1]) + matt_dev)
  } else {
    matt_new[, , 1] <- y$matt[i, , 1] * matt_dev
  }
  matt_new[is.na(matt_new)] <- 0

  return(list(matt = matt_new))
}

matt_Sarita <- Map(calc_Sarita_fedfry_matt, x = report[1:1000], y = report_RBT[1:1000], type = "ratio")
saveRDS(matt_Sarita, file = "CM/Sarita_maturity.rds")
matt_Sarita <- readRDS("CM/Sarita_maturity.rds")

rs_names <- c("Sarita Fed Fry (hypothesized)", "RBT Smolt 0+")
col_pal <- c("magenta", "black") %>% structure(names = rs_names)
g <- salmonMSE:::CM_maturity(
  matt_Sarita, list(n_r = 2, bmatt = rep(0, 5), Nages = 5),
  year1 = 2005, rs_names = rs_names, annual = TRUE) +
  scale_colour_manual(values = col_pal) +
  scale_fill_manual(values = col_pal)
ggsave("figures/Sarita_maturity.png", g, height = 5.5, width = 6)

# Take average from last 6 complete brood years
matt_avg <- sapply(matt_Sarita, function(x) {
  apply(x$matt[seq(9, 14), , ], 2:3, mean)
}, simplify = "array")

# Plot posterior density
g <- reshape2::melt(matt_avg) %>%
  rename(Maturity = value) %>%
  filter(Var1 %in% 2:4) %>%
  mutate(Age = paste("Age", Var1), `Release Strategy` = rs_names[Var2]) %>%
  ggplot(aes(Maturity, fill = `Release Strategy`, colour = `Release Strategy`)) +
  facet_wrap(vars(Age), ncol = 2, scales = "free") +
  geom_density(alpha = 0.25) +
  scale_colour_manual(values = col_pal) +
  scale_fill_manual(values = col_pal) +
  labs(y = "Posterior density") +
  theme(legend.position = "bottom")

#matplot(matt_avg[, 1, ], typ = 'l')  # Fed Fry
#matlines(matt_avg[, 2, ], typ = 'l') # Traditionals

