####### Libraries #######

library(RulesTools)
library(arules)
library(tidyr)




####### Analysis Parameters #######

NUM_TRANSACTIONS <- 126

covariates <-c(
  "eDNAConc",
  "eFishCatch",
  "AirTemp",
  "WaterTemp",
  "pH",
  "DissolvedOxygen",
  "Conductivity",
  "VolumeFiltered"
)

labels <- c(
  "eDNA Concentration (Copies/µL)",
  "Electrofish Catch",
  "Air Temperature (°C)",
  "Water Temperature (°C)",
  "pH",
  "Dissolved Oxygen (mg/L)",
  "Conductivity (µS)",
  "Volume Filtered (mL)"
)

discretizations <- c(
  13.3,                      # for eDNA concentrations (50% LOD)
  0,                         # for electrofish catch (absence/presence)
  15.03,                     # for air temp (avg temp for sample site in Sept)
  15,                        # for water temp (optimal temp for brook trout)
  7.75,                      # for pH (CCME guideline)
  9.5,                       # for dissolved oxygen (CCME guideline)
  1047,                      # for conductivity (median)***
  1.04                       # for volume filtered (median)
)

consequents <- c(
  "eDNAConc=high", 
  "eDNAConc=low", 
  "eFishCatch=present", 
  "eFishCatch=absent"
)




####### Data Pre-Processing #######

# init empty df
discretized_df <- data.frame(matrix(nrow=NUM_TRANSACTIONS, ncol=0))

# add discretized values to df
for(i in seq_along(discretizations)){
  
  if(covariates[i] == "eFishCatch") {
    temp <- dtize_col(
      column = BrookTrout[[covariates[i]]],
      cutoff = discretizations[i],
      infinity = TRUE,
      labels = c("absent", "present")
    )
  }
  
  else {
    temp <- dtize_col(
      column = BrookTrout[[covariates[i]]],
      cutoff = discretizations[i],
      infinity = TRUE
    )
  }
  
  discretized_df[[covariates[i]]] <- temp
  
}

# backpack and site already discrete, so they can be added to df directly
discretized_df[["Backpack"]] <- as.factor(BrookTrout$Backpack)
discretized_df[["Site"]] <- as.factor(BrookTrout$Site)

# convert df into type transactions
transactions <- as(discretized_df, "transactions")




####### Rule Mining #######

rules <- list()

for (con in consequents) {
  
  # mine rules
  raw_rules <- apriori(
    transactions,
    parameter = list(support = 1/NUM_TRANSACTIONS, confidence = 1/NUM_TRANSACTIONS),
    appearance = list(rhs = con)
  )
  
  # subset non-redundant rules
  pruned_rules <- raw_rules[!is.redundant(raw_rules, measure = "confidence")]
  
  # subset significant rules (unadjusted + multiple comparison corrections)
  unadj_rules = pruned_rules[is.significant(pruned_rules, alpha = 0.05)] %>%
    sort(by = c("lift", "confidence", "support"))
  
  BH_rules = pruned_rules[is.significant(pruned_rules, alpha = 0.05, adjust = "BH")] %>%
    sort(by = c("lift", "confidence", "support"))
  
  BF_rules = pruned_rules[is.significant(pruned_rules, alpha = 0.05, adjust = "bonferroni")] %>%
    sort(by = c("lift", "confidence", "support"))
  
  # final rules list
  rules[[con]] <- list(
    raw = raw_rules,
    pruned = pruned_rules,
    unadj = unadj_rules,
    BH = BH_rules,
    BF = BF_rules
  )
  
}
