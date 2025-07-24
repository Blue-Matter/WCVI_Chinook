

# List of state variables, parameters for each MCMC sample from Script #1
samp <- readRDS("CM/Sarita_RBT_CM_05.20.25.rds")

report <- salmonMSE:::get_report(samp)
fit <- samp@.MISC$CMfit
d <- salmonMSE:::get_CMdata(fit)

source("https://raw.githubusercontent.com/Pacific-salmon-assess/samEst/refs/heads/main/R/RP_functions.R")

# 0.08 - 0.14 seconds for 4,000 MCMC samples
tictoc::tic()
ref_pt2 <- sapply(1:length(report), function(i) {

  alpha <- report[[i]]$epro * report[[i]]$alpha # Recruits/spawner
  so <- report[[i]]$so # Unfished spawners
  Smax <- so/log(alpha)

  umsy <- umsyCalc(log(alpha))
  smsy <- smsyCalc(log(alpha), 1/Smax)
  sgen <- sgenCalcDirect(log(alpha), 1/Smax)

  ref_lambert <- structure(c(alpha, umsy, smsy, sgen), names = c("alpha", "UMSY", "SMSY", "Sgen"))
  return(ref_lambert)
})
tictoc::toc()

hist(ref_pt2["SMSY", ])
hist(ref_pt2["alpha", ])

hist(ref_pt2["UMSY", ])
hist(ref_pt2["Sgen", ])


