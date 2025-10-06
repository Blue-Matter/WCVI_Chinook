
#### Calculate productivity in the absence of fishery harvest
library(salmonMSE)
SOM <- readRDS("SOM/SOM_base.rds") # Make this object from running script 3
SOM <- check_SOM(SOM)

s <- 1
x <- 1
y <- 1

p_female <- SOM@Bio[[s]]@p_female
maxage <- SOM@Bio[[s]]@maxage
n_g <- SOM@Bio[[s]]@n_g
p_LHG <- SOM@Bio[[s]]@p_LHG

vulPT <- SOM@Harvest[[s]]@vulPT[x, ]
vulT <- SOM@Harvest[[s]]@vulT[x, ]

fec <- matrix(SOM@Bio[[s]]@fec[x, , y], maxage, n_g)
Mjuv <- matrix(SOM@Bio[[s]]@Mjuv_NOS[x, , y, ], maxage, n_g)
p_mature <- matrix(SOM@Bio[[s]]@p_mature[x, , y], maxage, n_g)

#### Egg per juvenile (phi) ----
surv_juv <- sapply(1:n_g, function(g) p_LHG[g] * salmonMSE:::calc_survival(Mjuv[, g], p_mature[, g])) # First semester due to exploitation
surv_return <- surv_juv * p_mature
surv_esc <- surv_return #* exp(-FT)
surv_spawn <- surv_esc #* s_enroute

EPJ <- sum(surv_spawn * p_female * fec) # Egg per juvenile

#### Juvenile per egg ----
Egg <- 1
surv_egg_fry <- mean(SOM@Habitat[[s]]@fry_sdev[x, ]) # 0.09
JPE <- Egg * surv_egg_fry

prod <- JPE * EPJ

