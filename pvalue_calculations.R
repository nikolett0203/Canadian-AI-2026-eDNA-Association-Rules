####### Libraries #######

library(RulesTools)
library(arules)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(dplyr)





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





####### Scatterplot Functions #######

scatter_pval <- function(df, x_var, y_var, size_var, color_var, x_lab, y_lab, alpha, title){
  
  df$neg_log_p <- -log10(df[[y_var]])
  neg_log_a <- -log10(alpha)
  
  ggplot(df, aes(x = .data[[x_var]], y = neg_log_p)) +
    geom_point(aes(size = .data[[size_var]], color = .data[[color_var]])) + 
    scale_color_viridis(option="D") +
    labs(x = x_lab, y = y_lab, title = title) +
    geom_hline(yintercept = neg_log_a, linetype = "dashed", color = "red") +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    ) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.01)
    ) +
    theme(
      plot.title = element_text(
        size = 12,
        hjust = 0.5    # centered title
      )
    ) +
    theme_bw()
  
}





combined_plots <- function(scatter_plots, heatmaps){
  
  for (i in 1:4){
    
    combined_plot <- scatterplots[[i]] | heatmaps[[i]]
    combined_plot <- combined_plot + plot_layout(widths = c(2, 1))
    
    print(combined_plot)
    
  }
  
}





####### DF Functions #######

create_df <- function(rules, method, pruned = FALSE) {
  
  # make sure we get the right pvalues for the adjustment type
  p_type <- switch(
    method,
    unadj = "pvalue",
    BH = "p_BH",
    BF = "p_BF"
  )
  
  if (pruned == TRUE) {
    subset = "pruned"
  } else {
    subset = method
  }
  
  do.call(rbind, lapply(names(rules), function(target) {
    
    ruleset <- rules[[target]][[subset]]
    qty <- quality(ruleset)
    
    data.frame(
      consequent = target,
      support = qty$support,
      confidence = qty$confidence,
      lift = qty$lift,
      p = qty[[p_type]],
      length = qty$rule_length,
      stringsAsFactors = FALSE
    )
  }))
}

summary_stats <- function(x) {
  c(
    min = min(x),
    med = median(x),
    max = max(x)
  )
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

plot_titles <- list(
  "{eDNAConc=high}",
  "{eDNAConc=low}",
  "{eFishCatch=present}",
  "{eFishCatch=absent}"
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
  pvalues <- get_pvalues(pruned_rules, transactions)
  
  # attach to rule quality
  quality(pruned_rules)$pvalue <- pvalues
  quality(pruned_rules)$p_BH <- p.adjust(pvalues, method = "BH")
  quality(pruned_rules)$p_BF <- p.adjust(pvalues, method = "bonferroni")
  quality(pruned_rules)$rule_length <- size(pruned_rules)
  quality(pruned_rules)$lhs <- labels(lhs(pruned_rules))
  
  # subset significant rules (unadjusted + multiple comparison corrections)
  unadj_rules <- sort(pruned_rules[is.significant(pruned_rules, alpha = 0.05, adjust = "none")],
                     by = c("lift","confidence","support"))
  
  BH_rules <- sort(pruned_rules[is.significant(pruned_rules, alpha = 0.05, adjust = "BH")],
                   by = c("lift","confidence","support"))
  
  BF_rules <- sort(pruned_rules[is.significant(pruned_rules, alpha = 0.05, adjust = "bonferroni")],
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

####### Undajusted Scatterplots #######

unadj_plots <- list()


for (i in 1:4){
  
  unadj_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$unadj), 
      x_var = "confidence", 
      y_var = "pvalue", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = "Unadjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(unadj_plots))





####### Benjamini-Hochberg Scatterplots #######

BH_plots <- list()


for (i in 1:4){

  BH_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BH), 
      x_var = "confidence", 
      y_var = "p_BH", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = "BH-Adjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(BH_plots))





####### Bonferroni Scatterplots #######

BF_plots <- list()

for (i in 1:4){
  
  BF_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BF), 
      x_var = "confidence", 
      y_var = "p_BF", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = "Bonferroni-Adjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(BF_plots))





####### Rule Length Scatterplots #######

BH_length <- list()

for (i in 1:4){
  
  BH_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BH), 
      x_var = "rule_length", 
      y_var = "p_BH", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = "BH-Adjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(BH_length))


BF_length <- list()

for (i in 1:4){
  
  BF_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BF), 
      x_var = "rule_length", 
      y_var = "p_BF", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = "Bonferroni-Adjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(BF_length))


unadj_length <- list()

for (i in 1:4){
  
  unadj_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$unadj), 
      x_var = "rule_length", 
      y_var = "pvalue", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = "Unadjusted P-Value", 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}

print(plots_2x2(unadj_length))





####### P-Value Heatmaps #######

BH_heatmaps <- list()

for (i in 1:4){
  BH_heatmaps[[i]] <-
    ggplot(
      quality(rules[[i]]$BH), 
      aes(x = "BH-Adjusted P-Value", y = reorder(lhs, -log10(p_BH)), fill = -log10(p_BH))
    ) +
    geom_tile(color = "white") +
    scale_fill_viridis(option = "magma", direction = -1, name = expression(-log[10](P[adj]))) +
    labs(x = "", y = "Antecedent") +
    theme(axis.text.y = element_text(size = 4))
  
}

print(BH_heatmaps)


BF_heatmaps <- list()

for (i in 1:4){
  BF_heatmaps[[i]] <-
    ggplot(
      quality(rules[[i]]$BH), 
      aes(x = "BH-Adjusted P-Value", y = reorder(lhs, -log10(p_BF)), fill = -log10(p_BH))
    ) +
    geom_tile(color = "white") +
    scale_fill_viridis(option = "magma", direction = -1, name = expression(-log[10](P[adj]))) +
    labs(x = "", y = "Antecedent") +
    theme(axis.text.y = element_text(size = 4))
  
}

print(BF_heatmaps)





####### Other Stats #######

unadj_sig <- create_df(rules, "unadj")
bh_sig    <- create_df(rules, "BH")
bf_sig    <- create_df(rules, "BF")

unadj_sig$method <- "unadj"
bh_sig$method    <- "BH"
bf_sig$method    <- "BF"

sig_rules <- rbind(unadj_sig, bh_sig, bf_sig)

summary_table <- aggregate(
  cbind(support, confidence, lift, p) ~ method,
  sig_rules,
  summary_stats
)

summary_table[ , -1] <- round(summary_table[ , -1], 3)

spearman_sig <- sig_rules %>%
  group_by(method) %>%
  summarise(across(
    all_of(c("support", "confidence", "lift", "length")),
    ~ cor(p, .x, method = "spearman"),
    .names = "rho_{.col}"
  ))

spearman_sig

spearman_sig_by_con <- sig_rules %>%
  group_by(method, consequent) %>%
  summarise(
    across(
      all_of(c("support", "confidence", "lift", "length")),
      ~ cor(p, .x, method = "spearman"),
      .names = "rho_{.col}"
    ),
    .groups = "drop"
  )

spearman_sig_by_con

####### Spearman on Pre-Significance Filtered Sets #######

unadj_nonredund <- create_df(rules, "unadj", TRUE)
bh_nonredund    <- create_df(rules, "BH", TRUE)
bf_nonredund    <- create_df(rules, "BF", TRUE)

unadj_nonredund$method <- "unadj"
bh_nonredund$method    <- "BH"
bf_nonredund$method    <- "BF"

nonredund_rules <- rbind(unadj_nonredund, bh_nonredund, bf_nonredund)

spearman_nonredund <- nonredund_rules %>%
  group_by(method) %>%
  summarise(across(
    all_of(c("support", "confidence", "lift", "length")),
    ~ cor(p, .x, method = "spearman"),
    .names = "rho_{.col}"
  ))

spearman_nonredund


spearman_nonred_by_con <- nonredund_rules %>%
  group_by(method, consequent) %>%
  summarise(
    across(
      all_of(c("support", "confidence", "lift", "length")),
      ~ cor(p, .x, method = "spearman"),
      .names = "rho_{.col}"
    ),
    .groups = "drop"
  )

spearman_nonred_by_con
