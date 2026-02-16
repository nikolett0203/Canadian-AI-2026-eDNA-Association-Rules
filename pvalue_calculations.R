####### Libraries #######

library(RulesTools)




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
  "eFishCatch=present" 
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
















####### Libraries #######


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





####### Plot Functions #######

scatter_pval <- function(df, x_var, y_var, size_var, color_var, x_lab, y_lab, alpha, title){
  
  df$neg_log_p <- -log10(df[[y_var]])
  neg_log_a <- -log10(alpha)
  
  ggplot(df, aes(x = .data[[x_var]], y = neg_log_p)) +
    geom_point(aes(size = .data[[size_var]], color = .data[[color_var]])) + 
    scale_color_viridis(option="D", labels = scales::label_number(accuracy = 0.01)) +
    scale_size_continuous(labels = scales::label_number(accuracy = 0.01)) +
    labs(x = x_lab, y = y_lab, title = title) +
    geom_hline(yintercept = neg_log_a, linetype = "dashed", color = "red") +
    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    ) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 0.01)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        size = 12,
      )
    )
  
}





combined_plots <- function(scatter_plots, heatmaps){
  
  for (i in 1:2){
    
    combined_plot <- scatter_plots[[i]] | heatmaps[[i]]
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

spearman_table <- function(df) {
  
  df %>%
    pivot_longer(
      cols = c(support, confidence, lift, length),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(method, metric) %>%
    summarise(
      {
        ct <- cor.test(p, value, method = "spearman")
        tibble(
          rho  = unname(ct$estimate),
          praw = ct$p.value,
          n    = sum(complete.cases(p, value))
        )
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_BH = p.adjust(praw, method = "BH"),
      p_BF = p.adjust(praw, method = "bonferroni")
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
  #"eDNAConc=low", 
  "eFishCatch=present" 
  #"eFishCatch=absent"
)

plot_titles <- list(
  "{eDNAConc=high}",
  #"{eDNAConc=low}",
  "{eFishCatch=present}"
  #"{eFishCatch=absent}"
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





####### Undajusted Plots #######

unadj_plots <- list()

for (i in 1:2){
  
  unadj_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$unadj), 
      x_var = "confidence", 
      y_var = "pvalue", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = expression(-log[10]("Uncorrected P-value")),
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}






####### Benjamini-Hochberg Plots #######

BH_plots <- list()

for (i in 1:2){
  
  BH_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BH), 
      x_var = "confidence", 
      y_var = "p_BH", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = expression(-log[10]("BH-Adjusted P-value")),
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}





####### Bonferroni Plots #######

BF_plots <- list()

for (i in 1:2){
  
  BF_plots[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BF), 
      x_var = "confidence", 
      y_var = "p_BF", 
      size_var = "support", 
      color_var = "lift", 
      x_lab = "Confidence", 
      y_lab = expression(-log[10]("BF-Adjusted P-value")),
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}





####### Rule Length Scatterplots #######

BH_length <- list()

for (i in 1:2){
  
  BH_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BH), 
      x_var = "rule_length", 
      y_var = "p_BH", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = expression(-log[10]("BH-Adjusted P-value")), 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}


BF_length <- list()

for (i in 1:2){
  
  BF_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$BF), 
      x_var = "rule_length", 
      y_var = "p_BF", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = expression(-log[10]("BF-Adjusted P-value")), 
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}


unadj_length <- list()

for (i in 1:2){
  
  unadj_length[[i]] <-
    scatter_pval(
      df = quality(rules[[i]]$unadj), 
      x_var = "rule_length", 
      y_var = "pvalue", 
      size_var = "lift", 
      color_var = "confidence", 
      x_lab = "Rule Length", 
      y_lab = expression(-log[10]("Uncorrected P-value")),
      alpha = 0.05,
      title = plot_titles[[i]]
    )
}





####### Other Stats #######

unadj_sig <- create_df(rules, "unadj")
bh_sig    <- create_df(rules, "BH")
bf_sig    <- create_df(rules, "BF")

unadj_sig$method <- "unadj"
bh_sig$method    <- "BH"
bf_sig$method    <- "BF"

sig_rules <- rbind(unadj_sig, bh_sig, bf_sig)

spearman_sig <- spearman_table(sig_rules)
spearman_sig

unadj_nonredund <- create_df(rules, "unadj", TRUE)
bh_nonredund    <- create_df(rules, "BH", TRUE)
bf_nonredund    <- create_df(rules, "BF", TRUE)

unadj_nonredund$method <- "unadj"
bh_nonredund$method    <- "BH"
bf_nonredund$method    <- "BF"

nonredund_rules <- rbind(unadj_nonredund, bh_nonredund, bf_nonredund)

spearman_nonredund <- spearman_table(nonredund_rules)
spearman_nonredund

spearman_table_by_con <- function(df) {
  
  df %>%
    pivot_longer(
      cols = c(support, confidence, lift, length),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(method, consequent, metric) %>%
    summarise(
      {
        ct <- cor.test(p, value, method = "spearman")
        tibble(
          rho  = unname(ct$estimate),
          praw = ct$p.value,
          n    = sum(complete.cases(p, value))
        )
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_BH = p.adjust(praw, method = "BH"),
      p_BF = p.adjust(praw, method = "bonferroni")
    )
}

spearman_sig_by_con <- spearman_table_by_con(sig_rules)
spearman_nonred_by_con <- spearman_table_by_con(nonredund_rules)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

write.csv(spearman_sig,
          file = file.path(out_dir, "spearman_sig.csv"),
          row.names = FALSE)

write.csv(spearman_nonredund,
          file = file.path(out_dir, "spearman_nonredund.csv"),
          row.names = FALSE)

write.csv(spearman_sig_by_con,
          file = file.path(out_dir, "spearman_sig_by_con.csv"),
          row.names = FALSE)

write.csv(spearman_nonred_by_con,
          file = file.path(out_dir, "spearman_nonred_by_con.csv"),
          row.names = FALSE)

