

.ts_fn <- function(SMSE, name, var, all_sims = FALSE) {
  require(salmonMSE)

  if (all_sims) {
    s <- 1

    if (var == "Brood") {
      res <- SMSE@NOB[, s, ] + SMSE@HOB[, s, ]
    } else if (var == "Egg") {
      res <- SMSE@Egg_NOS[, s, ] + SMSE@Egg_HOS[, s, ]
    } else if (var == "Mean age") {
      Sp <- SMSE@NOS[, 1, , ] + SMSE@HOS[, 1, , ]
      res <- apply(Sp, c(1, 3), function(w) weighted.mean(x = 1:5, w = w))
    } else if (var == "IR_Return") {
      res <- apply(SMSE@Escapement_NOS[, s, , ] + SMSE@Escapement_HOS[, s, , ], c(1, 3), sum)
    } else if (var == "IR_Catch") {
      res <- SMSE@Misc$inriver_catch
    } else {
      res <- plot_statevar_ts(SMSE, var, figure = FALSE, quant = FALSE)
    }
    dimnames(res) <- NULL

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

  lev <- levels(df$variable)

  df$variable <- sub(" ", "\n", df$variable)
  lev <- sub(" ", "\n", lev)
  df$variable <- factor(df$variable, lev)

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

decision_table_grid <- function(x, title = "PNI", ncol = 2,
                                fill_scheme = scale_fill_gradient2(low = "deeppink", high = "green4", mid = "white", midpoint = 0.5)) {

  g <- salmonMSE::plot_decision_table(
    x = x$ER,
    y = x$pNOB_target,
    z = x$value,
    scenario = x$Scenario,
    title = title,
    xlab = "ER",
    ylab = "pNOB target",
    ncol = ncol
  ) +
    fill_scheme

  #vars <- unique(x$Scenario)
#
  #dt <- data.frame(x = x$ER, y = x$pNOB_target, z = x$value, Scenario = x$Scenario)
  #dt$txt <- format(round(dt$z, 2))
#
  #g <- ggplot(dt, aes(factor(.data$x), factor(.data$y), fill = .data$z)) +
  #  geom_tile(colour = "grey20", alpha = 0.6) +
  #  geom_text(aes(label = .data$txt)) +
  #  coord_cartesian(expand = FALSE) +
  #  guides(fill = "none") +
  #  facet_wrap(vars(.data$Scenario), ncol = 2) +
  #  theme_bw() +
  #  theme(panel.grid.major = element_blank(),
  #        strip.background = element_blank(),
  #        panel.grid.minor = element_blank()) +
  #  fill_scheme +
  #  labs(x = "ER", y = "pNOB target", title = title)
  g
}

tradeoff_grid <- function(val_sim, xname = "Total Spawners", yname = "PNI", xlab = xname, ylab = yname,
                          xlim = NULL, ylim = NULL, ncol = 2) {

  g <- salmonMSE::plot_tradeoff(
    val_sim %>% filter(variable == xname) %>% select(lwr, median, upr) %>% as.matrix(),
    val_sim %>% filter(variable == yname) %>% select(lwr, median, upr) %>% as.matrix(),
    val_sim %>% filter(variable == xname) %>% pull(pNOB_target) %>% factor(),
    val_sim %>% filter(variable == xname) %>% pull(ER) %>% factor(),
    xlab = xlab,
    ylab = ylab,
    x1lab = "pNOB target",
    x2lab = "ER",
    scenario = val_sim %>% filter(variable == xname) %>% pull(Scenario),
    ncol = ncol
  ) +
    scale_shape_manual(values = c(1, 4, 16)) +
    coord_cartesian(xlim = xlim, ylim = ylim)

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


plot_spaghetti <- function(x, sims, OM_name = NULL, MP_name = NULL, alpha = 0.4, by_origin = FALSE) {

  if (by_origin) {
    require(ggborderline)

    meds <- summarise(x, value = median(value), .by = c(Year, var_name, Scenario, scenario, Origin))

    g <- x %>%
      mutate(gr = paste(Simulation, Origin)) %>%
      ggplot(aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line(alpha = alpha, aes(colour = Origin, group = factor(gr))) +
      #geom_line(data = meds, colour = "black", aes(group = Origin), linewidth = 1.5) +
      ggborderline::geom_borderline(data = meds, aes(colour = Origin), bordercolour = "grey40", linewidth = 1) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, colour = NULL) +
      ggtitle(OM_name, subtitle = MP_name)

  } else if (missing(sims)) { # All simulations

    meds <- summarise(x, value = median(value), .by = c(Year, var_name, Scenario, scenario))

    g <- ggplot(x, aes(Year, value)) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line(alpha = alpha, colour = "grey40", aes(group = factor(Simulation))) +
      geom_line(data = meds, colour = "blue", linewidth = 1) +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL) +
      ggtitle(OM_name, subtitle = MP_name)

  } else {

    val_plot <- dplyr::filter(x, Simulation %in% sims)
    g <- ggplot(val_plot, aes(Year, value, colour = factor(Simulation))) +
      facet_wrap(vars(var_name), scales = "free_y") +
      geom_line() +
      expand_limits(y = 0) +
      theme(strip.background = element_blank(), legend.position = "bottom") +
      labs(y = NULL, colour = "Simulation") +
      scale_colour_brewer(palette = "Dark2") +
      ggtitle(OM_name, subtitle = MP_name)
  }
  g
}

# Function that parses text for management options and typesets the ER (pass along to ggplot)
bold_scenario <- function(x) {
  xx <- strsplit(x, ",")
  xout <- sapply(xx, function(i) {
    if (grepl("= 1", i[1])) {
      paste0("italic(underline(\"", i[1], "\"))~\"", i[2], "\"")
    } else if (grepl("= 0.75", i[1])) {
      paste0("plain(\"", i[1], "\")~\"", i[2], "\"")
    } else if (grepl("= 0.5", i[1])) {
      paste0("bold(\"", i[1], "\")~\"", i[2], "\"")
    }
  })
  parse(text = xout)
}

