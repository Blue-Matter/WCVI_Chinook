

.ts_fn <- function(SMSE, name, var, all_sims = FALSE) {
  require(salmonMSE)

  if (all_sims) {
    s <- 1

    if (var == "Brood") {
      res <- SMSE@NOB[, s, ] + SMSE@HOB[, s, ]
    } else if (var == "Egg") {
      res <- SMSE@Egg_NOS[, s, ] + SMSE@Egg_HOS[, s, ]
    } else {
      res <- plot_statevar_ts(SMSE, var, figure = FALSE, quant = FALSE)
    }

  } else {

    if (var == "Brood") {
      out <- apply(SMSE@NOB + SMSE@HOB, 3, quantile, c(0.025, 0.5, 0.975))
    } else if (var == "Egg") {
      out <- apply(SMSE@Egg_NOS + SMSE@Egg_HOS, 3, quantile, c(0.025, 0.5, 0.975))
    } else {
      out <- plot_statevar_ts(SMSE, var, figure = FALSE, quant = TRUE)
    }

    res <- reshape2::melt(out) %>%
      rename(Year = Var2) %>%
      mutate(name = name) %>%
      reshape2::dcast(Year + name ~ Var1)
  }

  return(res)

}

ts_fn <- function(SMSE_list, name, var) {
  d <- Map(.ts_fn, SMSE = SMSE_list, name = name, MoreArgs = list(var = var)) %>%
    bind_rows() %>%
    filter(!is.na(`50%`), `50%` > 0) %>%
    mutate(var = .env$var)

  g <- ggplot(d, aes(Year, `50%`, colour = name, fill = name)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.25, colour = NA, linetype = 2) +
    geom_line() +
    labs(x = "Year", y = var, colour = "Scenario", fill = "Scenario")
  g
}

plot_dotplot <- function(val_sim) {
  g <- val_sim %>%
    ggplot(aes(scenario, median, ymin = lwr, ymax = upr, shape = factor(ER), colour = factor(pNOB_target))) +
    facet_wrap(vars(variable), scales = "free_x", strip.position = "top") +
    geom_point() +
    geom_linerange() +
    coord_flip()

  g
}


plot_table <- function(df, padding = 0.52) {

  d <- df %>%
    mutate(txt = signif(median, 3)) %>%
    mutate(txt = ifelse(txt < 0.01 & txt > 0, "<0.01", txt)) %>%
    mutate(val_rel = median/max(median, na.rm = TRUE),
           val_0_1 = (median - min(median, na.rm = TRUE)) / diff(range(median, na.rm = TRUE)),
           .by = variable)

  g <- ggplot(d, aes(variable, scenario)) +
    geom_tile(aes(fill = val_rel), alpha = 0.6, color = "white") +
    geom_text(aes(label = txt), size = ggplot2::rel(3)) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL) +
    coord_cartesian(
      expand = FALSE,
      xlim = range(as.numeric(d$variable)) + c(-padding, padding),
      ylim = range(as.numeric(d$scenario)) + c(-padding - 0.01, padding + 0.01)
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "grey10", angle = 90),
      strip.placement = "outside",
      strip.background = element_blank()
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels = levels(d$scenario)) +
    scale_fill_gradient2(low = "deeppink", high = "green4", mid = "white", limits = c(0, 1), midpoint = 0.5)
  g
}

decision_table_grid <- function(x, title = "PNI") {

  vars <- unique(x$Scenario)

  glist <- lapply(1:length(vars), function(i) {
    df <- x %>% dplyr::filter(Scenario == vars[i])
    g <- salmonMSE::plot_decision_table(
      x = df$ER,
      y = df$pNOB_target,
      z = df$value,
      xlab = "ER",
      ylab = "pNOB target"
    ) +
      labs(title = ifelse(i == 1, title, ""), subtitle = vars[i])
    g
  })

  g <- ggpubr::ggarrange(plotlist = glist, ncol = 2, nrow = 2)
  g
}

tradeoff_grid <- function(val_sim, xname = "Total Spawners", yname = "PNI", xlab = xname, ylab = yname,
                          xlim = NULL, ylim = NULL) {
  vars <- unique(val_sim$Scenario)

  glist <- lapply(1:length(vars), function(i) {

    d <- val_sim %>% dplyr::filter(Scenario == vars[i])
    d[is.na(d)] <- 0
    g <- salmonMSE::plot_tradeoff(
      d %>% filter(variable == xname) %>% select(lwr, median, upr) %>% as.matrix(),
      d %>% filter(variable == yname) %>% select(lwr, median, upr) %>% as.matrix(),
      d %>% filter(variable == xname) %>% pull(pNOB_target) %>% factor(),
      d %>% filter(variable == xname) %>% pull(ER) %>% factor(),
      xlab = xname,
      ylab = yname,
      x1lab = "pNOB target",
      x2lab = "ER"
    ) +
      scale_shape_manual(values = c(1, 4, 16)) +
      labs(subtitle = vars[i]) +
      coord_cartesian(xlim = xlim, ylim = ylim)

    g

  })
  g <- ggpubr::ggarrange(plotlist = glist, ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE)
  g
}



plot_histogram <- function(val, var = "PNI", binwidth = 0.025, scales = "free_y") {
  g <- val %>%
    ggplot(aes(value)) +
    geom_histogram(binwidth = binwidth, linewidth = 0.1, fill = "grey80", colour = "black") +
    facet_wrap(vars(scenario), scales = scales) +
    labs(x = var, y = "Frequency") +
    theme(strip.background = element_blank())
  g
}

