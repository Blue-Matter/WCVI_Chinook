
# Create csv files in script 5
library(Slick) # install.packages("Slick")
library(tidyverse)
library(reshape2)

# See article at: https://slick.bluematterscience.com/articles/DevelopersGuide.html


#### Create Slick object ----
slick <- Slick()

# Add metadata
slick@Title <- "salmonMSE results"
slick@Subtitle <- "Demonstration of visualization app for Sarita Chinook"
slick@Date <- Sys.Date()
slick@Author <- "Q. Huynh"
slick@Email <-  "[quang@bluematterscience.com](mailto:quang@bluematterscience.com)"
slick@Institution <- "Blue Matter Science"
slick@Introduction <- "

Demonstration of results to develop visualization app for salmonMSE.

"


# Load results ----
ts <- readr::read_csv("tables/Sarita_outcomes_sim_year.csv") %>%
  mutate(value = ifelse(variable == "Egg", value/1e6, value),
         value = ifelse(variable == "Smolt_Rel", value/1e5, value))
perf_sim <- readr::read_csv("tables/Sarita_outcomes_sim.csv")
perf_prob <- readr::read_csv("tables/Sarita_outcomes_prob.csv")


# Describe MPs (management options)
MP_vector <- unique(ts$scenario)

MP_code <- MP_vector %>% strsplit(" ") %>%
  sapply(function(x) {
    xx <- x[c(2, 4, 5, 7)]
    do.call(paste, as.list(xx))
  })

mps <- MPs(
  Code = MP_code,
  #Code = MP_vector,
  Label = MP_vector,
  Description = MP_vector,
  Preset = list(
    "All" = 1:length(MP_vector),
    "ER = 0.5" = c(1, 4, 7),
    "ER = 0.75" = c(2, 5, 8),
    "ER = 1" = c(3, 6, 9),
    "pNOB target = 0.5" = 1:3,
    "pNOB target = 0.75" = 4:6,
    "pNOB target = 1" = 7:9
  )
)
MPs(slick) <- mps # Add to slick object

# Describe operating models (scenarios)
OM_names <- unique(ts$Scenario)

oms <- OMs(
  Factors = data.frame(
    Factor = "Scenario",
    Level = OM_names,
    Description = c(
      "Base with 2013-2018 ocean exploitation rates (ER), average historical freshwater egg-fry survival",
      "Same as A with above average egg-fry survival",
      "Same as A with 25% lower ocean ER",
      "25% lower ocean ER (same as C) + earlier maturity/declining fecundity"
    )
  ),
  Design = data.frame(Scenario = OM_names) %>% `rownames<-`(paste("OM", LETTERS[1:4])),
  Preset = list(
    "All" = list(1:length(OM_names))
  )
)
OMs(slick) <- oms # Add to slick object

# Timeseries
ts_names <- unique(ts$variable)

ts_value <- ts %>%
  mutate(
    Scenario = factor(Scenario, OM_names),
    scenario = factor(scenario, MP_vector),
    variable = factor(variable, ts_names)
  ) %>%
  reshape2::acast(list("Simulation", "Scenario", "scenario", "variable", "Year"), value.var = "value")

ts_Description <- c(
  "Proportionate natural influence (0 - 1)",
  "Total Spawners (absolute abundance, > 0)",
  "Proportion natural broodtake (weighted by age-class fecundity, 0 - 1)",
  "Proportion effective hatchery spawners (0 - 1)",
  "Proportion wild spawners (0 - 1)",
  "Mean age of spawners (2-5 years)",
  "Brood (absolute abundance, > 0)",
  "In-river catch (absolute abundance, > 0)",
  "Escapement from marine fisheries (absolute abundance, > 0)",
  "Egg production in natural environment (millions of eggs)",
  "Hatchery releases (hundreds of thousands)"
)

timeseries <- Timeseries(
  Code = ts_names,
  Label = ts_names,
  Description = ts_Description,
  Time = 1:dim(ts_value)[5],
  TimeNow = 1,
  TimeLab = "Year",
  Value = ts_value,
  Preset = list(
    "All" = 1:length(ts_names)
  ),
  Target = rep(NA_real_, length(ts_names)),
  Limit = rep(NA_real_, length(ts_names))
)
Timeseries(slick) <- timeseries

# Performance metrics at end of simulation
PMs <- unique(perf_sim$variable)

PM_value <- perf_sim %>%
  mutate(
    Scenario = factor(Scenario, OM_names),
    scenario = factor(scenario, MP_vector),
    variable = factor(variable, PMs)
  ) %>%
  reshape2::acast(list("Simulation", "Scenario", "scenario", "variable"), value.var = "value")

PM_Description <- c(
  "Year 30: Proportionate natural influence (0 - 1)",
  "Year 30: Total Spawners (absolute abundance, > 0)",
  "Year 30: Proportion natural broodtake (weighted by age-class fecundity, 0 - 1)",
  "Year 30: Proportion effective hatchery spawners (0 - 1)",
  "Year 30: Proportion wild spawners (0 - 1)",
  "Year 30: Mean age of spawners (2-5 years)",
  "Year 30: Brood (absolute abundance, > 0)",
  "Year 30: In-river catch (absolute abundance, > 0)",
  "Year 30: Escapement from marine fisheries (absolute abundance, > 0)",
  "Year 30: Egg production in natural environment (millions of eggs)",
  "Year 30: Hatchery releases (hundreds of thousands)"
)

# Boxplot
boxplot <- Boxplot(
  Code = PMs,
  Label = PMs,
  Description = PM_Description,
  Value = PM_value,
  Preset = list(
    All = 1:length(PMs)
  )
)
Boxplot(slick) <- boxplot

# Quilt (performance table)
perf_quilt <- rbind(
  perf_sim %>% summarise(value = median(value, na.rm = TRUE), .by = c(Scenario, scenario, variable)),
  perf_prob %>% select(value, Scenario, scenario, variable)
)

quilt_code <- unique(perf_quilt$variable)
quilt_code <- quilt_code[c(1:2, 12:16, 3:11)]

quilt_value <- perf_quilt %>%
  mutate(
    Scenario = factor(Scenario, OM_names),
    scenario = factor(scenario, MP_vector),
    variable = factor(variable, quilt_code)
  ) %>%
  reshape2::acast(list("Scenario", "scenario", "variable"), value.var = "value")
Quilt_Description <- c(
  "Year 30 median: Proportionate natural influence (0 - 1)",
  "Year 30 median: Total Spawners (absolute abundance, > 0)",

  "Year 30 probability natural origin spawners exceed lower benchmark (250)",
  "Year 30 probability natural origin spawners exceed upper benchmark (476)",
  "Year 30 probability total spawners exceed lower benchmark (250)",
  "Year 30 probability total spawners exceed lower benchmark (476)",
  "Year 30 probability total spawners exceed 1300",

  "Year 30 median: Proportion natural broodtake (weighted by age-class fecundity, 0 - 1)",
  "Year 30 median: Proportion effective hatchery spawners (0 - 1)",
  "Year 30 median: Proportion wild spawners (0 - 1)",
  "Year 30 median: Mean age of spawners (2-5 years)",
  "Year 30 median: Brood (absolute abundance, > 0)",
  "Year 30 median: In-river catch (absolute abundance, > 0)",
  "Year 30 median: Escapement from marine fisheries (absolute abundance, > 0)",
  "Year 30 median: Egg production in natural environment (millions of eggs)",
  "Year 30 median: Hatchery releases (hundreds of thousands)"
)


quilt <- Quilt(
  Code = quilt_code,
  Label = quilt_code,
  Description = Quilt_Description,
  Value = quilt_value,
  Preset = list(All = 1:length(quilt_code))
)
Quilt(slick) <- quilt

# Spider (metric must be between 0 - 1!)
i <- c(1, 3:7, 8:10)
spider_Code <- spider_Label <- quilt_code[i]

spider_Description <- Quilt_Description[i]

spider_value <- quilt_value[, , i, drop = FALSE]
spider_value[is.na(spider_value)] <- 0

spider <- Spider(
  Code = spider_Code,
  Label = spider_Label,
  Description = spider_Description,
  Value = spider_value,
  Preset = list(All = 1:length(spider_Code))
)
Spider(slick) <- spider

# Tradeoff plot
tradeoff <- Tradeoff(
  Code = quilt_code,
  Label = quilt_code,
  Description = Quilt_Description,
  Value = quilt_value,
  Preset = list(All = 1:length(quilt_code))
)
Tradeoff(slick) <- tradeoff

# Save object
saveRDS(slick, file = "Slick/Sarita_09.29.2025.slick")

# Open App
#slick <- readRDS(file = "Slick/Sarita_09.29.2025.slick")
#Slick::App(slick = slick)
