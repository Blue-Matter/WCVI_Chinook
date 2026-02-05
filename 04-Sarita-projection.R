

library(tidyverse)

# Create grid of scenarios once and save ----
if (FALSE) {
  g <- expand.grid(
    MSF_T = c(FALSE, TRUE),
    MM = c(FALSE, TRUE)
  ) %>%
    mutate(Number = 1:n())

  g2 <- lapply(1:nrow(g), function(i) {
    data.frame(fs = c(0.1, 0.15, 0.2)) %>%
      mutate(Letter = LETTERS[1:3]) %>%
      cbind(g[i, ])
  }) %>%
    bind_rows() %>%
    arrange(fs, MM, MSF_T)

  g3 <- lapply(1:nrow(g2), function(i) {
    IRER <- c(0, 0.25, 0.75)
    pNOB <- c(0.5, 1)

    g_i <- g2[i, ]

    if (!g_i$MM) {
      pNOB <- NA_real_
      #IRER <- 0
    }

    expand.grid(IRER = IRER, pNOB_target = pNOB) %>%
      cbind(g_i) %>%
      mutate(Option_number = 1:n())

  }) %>%
    bind_rows() %>%
    mutate(Scenario = paste0(Letter, Number)) %>%
    mutate(Option_name = paste0("(", Option_number, ") IRER = ", IRER, ifelse(is.na(pNOB_target), "", paste0(", pNOB = ", pNOB_target)))) %>%
    mutate(Scenario_name = paste0(Scenario, ". ", fs, " fs, ",
                                  ifelse(MM, "MM, ", "no MM, "),
                                  ifelse(MSF_T, "MSF_T", "no MSF_T"))) %>%
    mutate(n = 1:n())

  readr::write_csv(g3, file.path("tables", "Sarita_scenarios.csv"))
}


g <- readr::read_csv(file.path("tables", "Sarita_scenarios.csv"))
nOM <- nrow(g)


# Define brood function
# Code assumes zero strays and mark rate = 1!
f_brood <- function(NO, HO, stray, m = 1, ptarget_NOB = 0.5, pmax_NOB = 0.5, pmax_esc = 1/3, min_esc = 600) {

  if (!m %in% c(0, 1)) stop("Brood function assumes mark rate of either zero or one.")

  NOB <- array(0, dim(NO))
  HOB_marked <- HOB_unmarked <- array(0, dim(HO))
  HOB_stray <- array(0, dim(stray))

  # Rule 1: no brood if fewer than 600 returns
  if (sum(NO, HO) > min_esc) {
    # Rule 2: Brood <= 33% of in-river return
    max_brood <- pmax_esc * sum(NO, HO)

    if (m == 1) {
      # Rule 3: total brood cannot exceed twice the natural brood available.
      # This means when natural fish are scarce, total brood scales down.
      brood_total_cap <- min(max_brood, 2 * sum(NO))

      # Take as many natural-origin fish as possible (up to the brood cap)
      NOB_total <- min(pmax_NOB * sum(NO), brood_total_cap)

      # Fill remaining brood with hatchery fish to reach ptarget_NOB
      # ptarget = NOB/(NOB + HOB) --> HOB = (NOB - ptarget * NOB)/ptarget
      HOB_total <- min((NOB_total - ptarget_NOB * NOB_total)/ptarget_NOB, brood_total_cap - NOB_total)

      # Safety: ensure hatchery never exceeds natural due to numerical jitter
      #pNOB <- NOB_total/(NOB_total + HOB_total)
      #if (pNOB < ptarget_NOB) {
      #  excess <- brood_hatchery - brood_natural
      #  brood_hatchery <- brood_hatchery - excess
      #}

      if (sum(NO)) {
        ptake_NOB <- NOB_total/sum(NO)
        NOB[] <- ptake_NOB * NO
      }

      if (sum(HO)) {
        ptake_HOB <- HOB_total/sum(HO)
        HOB_marked[] <- ptake_HOB * HO
      }
    } else if (m == 0) {

      # Rule 3 does not apply (use Rule 2)
      brood_total_cap <- max_brood

      pHOS <- sum(HO)/sum(NO, HO)
      pNOB <- 1 - pHOS

      NOB_total <- pNOB * brood_total_cap
      HOB_total <- (1 - pNOB) * brood_total_cap

      if (sum(NO)) {
        ptake_NOB <- NOB_total/sum(NO)
        NOB[] <- ptake_NOB * NO
      }

      if (sum(HO)) {
        ptake_HOB <- HOB_total/sum(HO)
        HOB_unmarked[] <- ptake_HOB * HO
      }
    } else {
      stop("Brood rule only accommodates mark rate of zero or 1")
    }
  }

  list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
}

# Define catch rule
# Code assumes mark rate = 1!
f_catch <- function(NO, HO, m = 1, targetER = 0.3, min_esc = 1300, output = c("HO", "NO")) {

  output <- match.arg(output)
  total_esc <- sum(NO, HO)

  if (total_esc < min_esc) {
    ER <- 0
  } else {

    if (m == 1) {

      if (output == "HO") {
        HO_notavail <- max(0, min_esc - sum(NO))
        HO_avail <- sum(HO) - HO_notavail

        Catch <- targetER * HO_avail
        ER <- Catch/sum(HO)
      } else {
        ER <- 0
      }

    } else if (m == 0) {

      total_catch <- targetER * (total_esc - min_esc)
      pNO <- sum(NO)/sum(NO, HO)

      if (output == "HO") {
        Catch_HO <- (1 - pNO) * total_catch
        ER <- sum(Catch_HO)/sum(HO)
      } else {
        Catch_NO <- pNO * total_catch
        ER <- sum(Catch_NO)/sum(NO)
      }

    } else {
      stop("Catch rule only accommodates mark rate of zero or 1")
    }
  }
  return(ER)
}


### Do simulation
library(snowfall)
library(pbapply)
sfInit(TRUE, 12)

sfExport(list = c("f_catch", "f_brood"))

sim_function <- function(i, g) {
  require(salmonMSE)

  # Freshwater survival
  if (g$fs[i] == 0.1) {
    SOM <- readRDS(file.path("SOM", "SOM_lowsurv.rds"))
  } else if (g$fs[i] == 0.15) {
    SOM <- readRDS(file.path("SOM", "SOM_mediumsurv.rds"))
  } else {
    SOM <- readRDS(file.path("SOM", "SOM_highsurv.rds"))
  }

  # Brood rule and in-river catch rule
  SOM@Hatchery@f_brood <- f_brood
  SOM@Hatchery@premove_HOS <- f_catch

  # Set in-river ER in catch rule and pNOB target in brood rule
  formals(SOM@Hatchery@premove_HOS)$targetER <- g$IRER[i]
  formals(SOM@Hatchery@f_brood)$ptarget_NOB <- g$pNOB_target[i]

  # Mass marking of Sarita fish
  # If no MM, then mark rate = 0 and add in-river fishery for natural-origin
  if (!g$MM[i]) {
    SOM@Hatchery@m <- 0
    formals(SOM@Hatchery@f_brood)$m <- 0

    SOM@Hatchery@premove_NOS <- f_catch
    formals(SOM@Hatchery@premove_NOS)$targetER <- g$IRER[i]
    formals(SOM@Hatchery@premove_NOS)$output <- "NO"
  }

  # Mark-selective terminal marine fishery
  # It is possible to have MSF with no mass-marking, Sarita ER = 0 for NO, HO
  if (g$MSF_T[i]) {
    if (g$MM[i]) {
      SOM@Harvest@MSF_T <- TRUE
      SOM@Harvest@release_mort <- c(0.05, 0.05)
    } else { # If MSF and no MM then just specify the harvest rate
      SOM@Harvest@u_terminal <- 0.05
    }
  }

  # No strays for now
  SOM@Hatchery@stray_external[] <- 0

  # Run projection
  SMSE <- salmonMSE(SOM)

  # Report in-river catch and ER
  SMSE@Misc$inriver_catch <- list(
    HOS = salmonMSE:::get_salmonMSE_var(salmonMSE_env$H, "HOS_remove"),
    NOS = salmonMSE:::get_salmonMSE_var(salmonMSE_env$N, "NOS_remove")
  )
  SMSE@Misc$ER_inriver <- list(
    HOS = SMSE@Misc$inriver_catch$HOS/(apply(SMSE@Escapement_HOS[, 1, , 1:29], c(1, 3), sum) - SMSE@HOB[, 1, 1:29]),
    NOS = SMSE@Misc$inriver_catch$NOS/(apply(SMSE@Escapement_NOS[, 1, , 1:29], c(1, 3), sum) - SMSE@NOB[, 1, 1:29])
  )
  SMSE@Misc$ER_inriver$HOS[is.na(SMSE@Misc$ER_inriver$HOS)] <- 0
  SMSE@Misc$ER_inriver$NOS[is.na(SMSE@Misc$ER_inriver$NOS)] <- 0

  # Save object
  saveRDS(SMSE, file = file.path("SMSE", paste0("Sarita", i, ".rds")))

  invisible()
}

# Run simulation
pblapply(7:nOM, sim_function, g = g, cl = snowfall::sfGetCluster())
sfStop()
