

# First, create grid of pNOB target and in-river exploitation rate target
# Default
library(tidyverse)
g_init <- expand.grid(
  ER = c(0.5, 0.75, 0.999),
  pNOB_target = c(0.5, 0.75, 0.99),
  ocean_ER_scalar = 1,
  surv = c("avg", "high", "bootstrap"),
  fec = "constant"
)

g <- rbind(
  g_init,
  g_init %>% filter(surv == "avg") %>% mutate(ocean_ER_scalar = 0.75),
  g_init %>% filter(surv == "avg") %>% mutate(ocean_ER_scalar = 1, fec = "decline")
) %>%
  mutate(n = 1:n())

nOM <- nrow(g)


# Define brood function
# Code assumes zero strays and mark rate = 1!
f_brood <- function(NO, HO, stray, m = 1, ptarget_NOB = 0.5, pmax_NOB = 0.5, pmax_esc = 1/3, min_esc = 600) {
  NOB <- array(0, dim(NO))
  HOB_marked <- HOB_unmarked <- array(0, dim(HO))
  HOB_stray <- array(0, dim(stray))

  # Rule 1: no brood if fewer than 600 returns
  if (sum(NO, HO) > min_esc) {
    # Rule 2: Brood <= 33% of return
    max_brood <- pmax_esc * sum(NO, HO)

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

  }

  list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
}

# Define catch rule
# Code assumes mark rate = 1!
f_catch <- function(NO, HO, m, targetER = 0.3, min_esc = 1300) {
  total_esc <- sum(NO, HO)
  if (total_esc < min_esc) {
    ER <- 0
  } else {
    HO_notavail <- max(0, min_esc - sum(NO))
    HO_avail <- sum(HO) - HO_notavail

    Catch <- targetER * HO_avail
    ER <- Catch/sum(HO)
  }
  return(ER)
}


### Start simulation
library(snowfall)
sfInit(TRUE, 9)

sfExport(list = c("f_catch", "f_brood"))

# Run in batches of 9
for (j in 1:5) {
  runs <- 9 * (j - 1) + seq(1, 9)
  print(runs)

  tictoc::tic()
  SMSE_list <- sfLapply(runs, function(i, g) {
    require(salmonMSE)

    if (g$surv[i] == "avg") {
      SOM <- readRDS(file.path("SOM", "SOM_base.rds"))
    } else if (g$surv[i] == "high") {
      SOM <- readRDS(file.path("SOM", "SOM_highsurv.rds"))
    } else if (g$surv[i] == "bootstrap") {
      SOM <- readRDS(file.path("SOM", "SOM_surv_bootstrap.rds"))
    }

    if (g$fec[i] == "decline") {
      SOM_fecdecline <- readRDS(file.path("SOM", "SOM_declinematfec.rds"))
      SOM@Bio@fec <- SOM@Hatchery@fec_brood <- SOM_fecdecline@Bio@fec

      SOM@Bio@p_mature <- SOM_fecdecline@Bio@p_mature
      SOM@Hatchery@p_mature_HOS <- SOM_fecdecline@Hatchery@p_mature_HOS
    }

    SOM@Hatchery@f_brood <- f_brood
    SOM@Hatchery@premove_HOS <- f_catch

    formals(SOM@Hatchery@f_brood)$ptarget_NOB <- g$pNOB_target[i]
    formals(SOM@Hatchery@premove_HOS)$targetER <- g$ER[i]

    SOM@Hatchery@stray_external[] <- 0

    SOM@Harvest@u_terminal <- g$ocean_ER_scalar[i] * SOM@Harvest@u_terminal
    SOM@Harvest@u_preterminal <- g$ocean_ER_scalar[i] * SOM@Harvest@u_preterminal

    SMSE <- salmonMSE(SOM)

    SMSE@Misc$inriver_catch <- salmonMSE:::get_salmonMSE_var(salmonMSE_env$H, "HOS_remove")
    SMSE@Misc$ER_inriver <- SMSE@Misc$inriver_catch/(apply(SMSE@Escapement_HOS[, 1, , 1:29], c(1, 3), sum) - SMSE@HOB[, 1, 1:29])
    SMSE@Misc$ER_inriver[is.na(SMSE@Misc$ER_inriver)] <- 0

    saveRDS(SMSE, file = file.path("SMSE", paste0("Sarita", i, ".rds")))

    invisible()

  }, g = g)
  tictoc::toc()
}


#tictoc::tic()
#SMSE_list <- sfLapply(1:nrow(gadd), function(i, g) {
#  require(salmonMSE)
#
#  if (grepl("HighSurv", g$name[i])) {
#    SOM <- readRDS(file.path("SOM", "SOM_highsurv.rds"))
#  } else {
#    SOM <- readRDS(file.path("SOM", "SOM_base.rds"))
#  }
#
#  SOM@Harvest@u_terminal <- g$scalar[i] * SOM@Harvest@u_terminal
#  SOM@Harvest@u_preterminal <- g$scalar[i] * SOM@Harvest@u_preterminal
#
#  #SOM@Hatchery@stray_external[] <- 0
#
#  if (g$name[i] == "Traditionals") {
#    SOM@Hatchery@n_yearling <- c(0.1, 0.9) * 500000
#  } else if (g$name[i] == "Fed Fry") {
#    SOM@Hatchery@n_yearling <- c(0.9, 0.1) * 500000
#  } else if (grepl("NoHarvestNoHatchery", g$name[i])) {
#    SOM@Hatchery@n_yearling[] <- 0
#    SOM@Hatchery@stray_external[] <- 0
#
#    SOM@Harvest@u_preterminal <-
#      SOM@Harvest@u_terminal <- 1e-8
#
#  } else if (g$name[i] == "DensityDependence") {
#
#    SOM@Habitat@prespawn_rel <- "HS"
#    SOM@Habitat@prespawn_prod <- 1
#    SOM@Habitat@prespawn_capacity <- 1300
#
#  } else if (g$name[i] == "DeclineMaturity") {
#    SOM <- readRDS(file.path("SOM", "SOM_declinematfec.rds"))
#  }
#
#  SMSE <- salmonMSE(SOM)
#  saveRDS(SMSE, file = file.path("SMSE", paste0("Sarita", g$OM[i], ".rds")))
#
#  invisible()
#}, g = gadd)
#
#
#sfStop()
#tictoc::toc()

#SOM <- readRDS("SOM/SOM_base_3sim.rds")
#SMSE <- salmonMSE(SOM)
#saveRDS(SMSE, file = "SMSE/Sarita_base2.rds")

#lhg_name <- c("Early Smalls", "Late Larges")
#rs_name <- c("Fed Fry", "Smolt 0+")
#report(SMSE, dir = "SMSE", filename = "Sarita_base2")
