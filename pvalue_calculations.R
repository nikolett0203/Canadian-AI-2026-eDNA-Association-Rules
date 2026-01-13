####### Libraries #######

library(RulesTools)
library(arules)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)





####### P-Value Function #######

# fix for edge cases, non-integer contingency cells, etc
get_pvalues <- function(rules, transactions) {
  
  # stores support, confidence, lift, count, etc. because we want count
  quality_df <- quality(rules)
  
  n <- length(transactions)
  
  p_values <- sapply(1:length(rules), function(i) {
    
    rule <- rules[i]
    
    lhs_items <- lhs(rule)
    rhs_items <- rhs(rule)
    
    # figure out how many transactions have rhs and lhs together
    a <- quality_df$count[i]
    
    # figure out how many had lhs only 
    # support returns fraction of itemsets with rule so we *n to get raw count
    lhs_count <- support(lhs_items, transactions) * n
    b <- lhs_count - a
    
    # rhs only
    rhs_count <- support(rhs_items, transactions) * n
    c <- rhs_count - a
    
    # count of transactions with neither lhs nor rhs
    d <- n - a - b - c
    
    # build contingency table 
    cont_table <- matrix(c(a, b, c, d), nrow = 2)
    
    # do the test
    fisher.test(cont_table, alternative = "greater")$p.value
  })
  
  return(p_values)
}






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

# convert df into transactions
transactions <- as(discretized_df, "transactions")





####### Rule Mining #######

rules <- list()

# loop for each targeted consequent
for (con in consequents) {
  
  # mine rules
  raw_rules <- apriori(
    transactions,
    parameter = list(support = 1/NUM_TRANSACTIONS, confidence = 1/NUM_TRANSACTIONS),
    appearance = list(rhs = con)
  )
  
  # subset non-redundant rules
  pruned_rules <- raw_rules[!is.redundant(raw_rules, measure = "confidence")]
  
  # isolate p-values
  pvalues = get_pvalues(pruned_rules, transactions)
  
  # attach to rule quality
  quality(pruned_rules)$pvalue <- pvalues
  quality(pruned_rules)$p_BH <- p.adjust(pvalues, method = "BH")
  quality(pruned_rules)$p_BF <- p.adjust(pvalues, method = "bonferroni")
  
  # subset significant rules (unadjusted + multiple comparison corrections)
  unadj_rules = sort(pruned_rules[quality(pruned_rules)$pvalue <= 0.05],
                     by = c("lift","confidence","support"))
  
  BH_rules <- sort(pruned_rules[quality(pruned_rules)$p_BH <= 0.05],
                   by = c("lift","confidence","support"))
  
  BF_rules <- sort(pruned_rules[quality(pruned_rules)$p_BF <= 0.05],
                   by = c("lift","confidence","support"))
  
  # final rules list
  rules[[con]] <- list(
    raw = raw_rules,
    pruned = pruned_rules,
    unadj = unadj_rules,
    BH = BH_rules,
    BF = BF_rules
  )
  
}





####### P-Value Scatterplots #######

scatterplots <- list()

for (i in 1:4){
  
  subset <- quality(rules[[i]]$BH)
  
  scatterplots[[i]] <- 
    ggplot(subset, aes(x=confidence, y=p_BH)) +
    geom_point(aes(size=support, color=lift)) +
    scale_color_viridis(option="D") +
    labs(x="Confidence", y="BH-Adjusted P-Value") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    )
  
}

# final graphic
(scatterplots[[1]] + scatterplots[[2]] + scatterplots[[3]] + scatterplots[[4]]) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")



