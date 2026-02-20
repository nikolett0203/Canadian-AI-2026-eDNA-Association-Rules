####### Libraries #######

library(arules)
library(tidyverse)
library(ggplot2)
library(patchwork)



####### Analysis Parameters #######

# Number of transactions in the BrookTrout dataset
NUM_TRANSACTIONS <- 126

# Continuous variables from the BrookTrout dataset 
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

# Discretization thresholds for converting continuous variables into {low, high}
# (and converting eFishCatch into {absent, present})
discretizations <- c(     # Thresholds chosen:
  13.3,                   # eDNAConc: 50% limit of detection from the molecular assay
  0,                      # eFishCatch: absence/presence threshold 
  15.03,                  # AirTemp: Sept. '19 average daily temperature at Hanlon Creek
  15,                     # WaterTemp: optimal temperature for brook trout
  7.75,                   # pH: CCME longterm freshwater guideline
  9.5,                    # DissolvedOxygen: CCME longterm freshwater guideline
  1047,                   # Conductivity: median from this dataset
  1.04                    # VolumeFiltered: median from this dataset
)

# Targeted consequents, selected for their relevance to domain experts in
# improving species detection success
consequents <- c(
  "eDNAConc=high", 
  "eFishCatch=present" 
)




####### Helper Functions #######

# Simple discretization helper for continuous variables
dtize_col <- function(column, cutoff, labels = c("low", "high")) {
  stopifnot(is.numeric(column), length(cutoff) == 1, is.numeric(cutoff))
  factor(ifelse(column <= cutoff, labels[1], labels[2]), levels = labels)
}

# Compute one-sided Fisher exact-test p-values for a rule set
# This mimics is.significant() from arules but returns the p-values for analysis
get_p <- function(rules, transactions) {
  
  if (length(rules) == 0) 
    return(numeric(0))
  
  quality_df <- quality(rules)
  n <- length(transactions)
  
  # Cache unique LHS/RHS itemsets and supports to avoid repeated support() calls
  lhs <- unique(lhs(rules))
  rhs <- unique(rhs(rules))
  
  lhs_supp <- support(lhs, transactions) * n
  rhs_supp <- support(rhs, transactions) * n
  
  # Create lookup vectors keyed by itemset label for O(1) access inside the loop
  lhs_labels <- labels(lhs)
  rhs_labels <- labels(rhs)
  names(lhs_supp) <- lhs_labels
  names(rhs_supp) <- rhs_labels
  
  # Construct a 2x2 contingency table for each rule, X -> Y:
  #   [a b]
  #   [c d]
  # where a = count(X & Y), b = count(X & !Y), c = count(!X & Y), d = count(!X & !Y)
  # Use one-sided Fisher test for positive association: P(Y|X) > P(Y)
  pvalues <- vapply(1:length(rules), function(i) {
    
    a <- quality_df$count[i]
    
    lhs_label <- labels(lhs(rules[i]))
    b <- lhs_supp[lhs_label] - a
    
    rhs_label <- labels(rhs(rules[i]))
    c <- rhs_supp[rhs_label] - a
    
    d <- n - a - b - c
    
    fisher.test(matrix(c(a, b, c, d), nrow = 2), 
                alternative = "greater")$p.value
    
  }, FUN.VALUE = numeric(1))
  
  return(pvalues)

}

# Filter rules at alpha = 0.05 then sort remaining rules
# by lift, confidence, and support (descending)
sig_sort <- function(rules, method) {
  sort(rules[is.significant(rules, alpha = 0.05, adjust = method)],
       by = c("lift", "confidence", "support"))
}

# Prepare longform dataframe with raw/BH/BF p-values stacked into one column
# Calculate -log10(p) for plotting
prep_for_plots <- function(rules, con) {
  rules[[con]]$pruned %>%
    quality() %>%
    pivot_longer(
      cols      = c(p_raw, p_BH, p_BF),
      names_to  = "p_type",
      values_to = "p_value" 
    ) %>%
    mutate(
      neg_log_p = -log10(p_value),
      p_type    = factor(p_type,
                         levels = c("p_raw", "p_BH", "p_BF"),
                         labels = c("Uncorrected", "Benjamini-Hochberg", "Bonferroni"))
    )
}

# Generate scatterplots of -log10(p) vs selected rule metrics for each correction
make_fig <- function(data, x_vars, x_labs, title) {
  
  plots <- Map(function(x_var, x_lab) {
    
    ggplot(data, aes_string(x_var, "neg_log_p")) +
      facet_wrap(~p_type, nrow = 1) + 
      geom_point(alpha = 0.6, size = 1.2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
      labs(
        x = x_lab,
        y = expression(-log[10](italic(p)))
      ) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.spacing.x = grid::unit(0.8, "lines")
      )
  }, x_vars, x_labs)
  
  wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

# Create a single dataframe of rule metrics for correlation analysis
prep_for_stats <- function(rules, method) {
  
  p_type <- switch(
    method,
    uncorr = "p_raw",
    BH    = "p_BH",
    BF    = "p_BF"
  )
  
  do.call(rbind, lapply(names(rules), function(target) {
    qty <- quality(rules[[target]]$pruned)
    data.frame(
      consequent = target,
      support    = qty$support,
      confidence = qty$confidence,
      lift       = qty$lift,
      p          = qty[[p_type]],
      length     = qty$len,
      stringsAsFactors = FALSE
    )
  }))
}

# Compute Spearman correlations between p-values and rule interestingness metrics
# Correlation p-values are BF- and BH-adjusted, but only BH-adjusted values are
# used in the final analysis
calc_coeffs <- function(data, by_con) {
  
  group_vars <- if(by_con) c("method", "consequent", "metric")
                else c("method", "metric")
  
  data %>%
    pivot_longer(
      cols = c(support, confidence, lift, length),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise({
      ct <- cor.test(p, value, method = "spearman")
      tibble(
        rho  = unname(ct$estimate),
        praw = ct$p.value,
        n    = sum(complete.cases(p, value))
      )
    }, .groups = "drop") %>%
    mutate(
      p_BH = p.adjust(praw, method = "BH"),
      p_BF = p.adjust(praw, method = "bonferroni")
    )
}




####### Data Pre-Processing #######

discretized_df <- data.frame(matrix(nrow=NUM_TRANSACTIONS, ncol=0))

# Load Brook Trout dataset from local CSV
BrookTrout <- read.csv("data/BrookTrout.csv", stringsAsFactors = FALSE)

# Discretize continuous variables into low/high categories
# or (absent/present for eletrofishing) 
for(i in seq_along(discretizations)){
  
  if(covariates[i] == "eFishCatch") {
    temp <- dtize_col(
      column = BrookTrout[[covariates[i]]],
      cutoff = discretizations[i],
      labels = c("absent", "present")
    )
  }
  
  else {
    temp <- dtize_col(
      column = BrookTrout[[covariates[i]]],
      cutoff = discretizations[i]
    )
  }
  
  discretized_df[[covariates[i]]] <- temp
  
}

# Backpack (i.e. sampler type) and site are already discrete
discretized_df[["Backpack"]] <- as.factor(BrookTrout$Backpack)
discretized_df[["Site"]]     <- as.factor(BrookTrout$Site)

# Convert discretized data into transaction format required by arules apriori() 
transactions <- as(discretized_df, "transactions")




####### Rule Mining #######

rules <- list()

for (con in consequents) {
  
  # Mine rules for each target: {eDNAConc=high} and {eFishCatch=present}
  # Thresholds of 1/n were chosen to capture rare rules for exploratory analysis +
  # to generate a large candidate rule set to demonstrate significance pruning
  raw_rules <- apriori(
    transactions,
    parameter = list(support = 1/NUM_TRANSACTIONS, confidence = 1/NUM_TRANSACTIONS),
    appearance = list(rhs = con)
  )
  
  # Remove redundant rules
  pruned_rules <- raw_rules[!is.redundant(raw_rules, measure = "confidence")]
  
  # Compute raw Fisher p-values (mirroring arules is.significant())
  p_values <- get_p(pruned_rules, transactions)
  
  # Save the uncorrected, Benjamini-Hochberg, and Bonferroni corrected rules
  # Also save rule length and antecedent itemsets
  quality(pruned_rules)$p_raw   <- p_values
  quality(pruned_rules)$p_BH    <- p.adjust(p_values, method = "BH")
  quality(pruned_rules)$p_BF    <- p.adjust(p_values, method = "bonferroni")
  quality(pruned_rules)$len     <- size(pruned_rules)
  quality(pruned_rules)$lhs     <- labels(lhs(pruned_rules))
  
  # Save raw rules, non-redundant rules, and significance-filtered subsets
  rules[[con]] <- list(
    raw    = raw_rules,
    pruned = pruned_rules,
    uncorr  = sig_sort(pruned_rules, "none"),
    BH     = sig_sort(pruned_rules, "BH"),
    BF     = sig_sort(pruned_rules, "bonferroni")
  )
  
}




####### Scatterplots #######

# For each consequent, visualize rule metrics vs significance under
# uncorrected/BH/Bonferroni p-values
for (con in consequents) {
  
  data <- prep_for_plots(rules, con)
  
  fig_con_sup <- make_fig(
    data,
    x_vars = c("confidence", "lift"),
    x_labs = c("Confidence", "Lift"),
    title  = paste0("{", con, "}")
  )
  
  fig_lift_len <- make_fig(
    data,
    x_vars = c("support", "len"),
    x_labs = c("Support", "Rule Length"),
    title  = paste0("{", con, "}")
  )
  
  print(fig_con_sup)
  print(fig_lift_len)
  
}




####### Correlation Analysis #######

# Compute Spearman correlations between p-values and interestingness measures
# on the non-redundant rule sets per correction method
methods <- c("uncorr", "BH", "BF")

nonredund_rules <- do.call(rbind, lapply(methods, function(m) {
  df <- prep_for_stats(rules, m)
  df$method <- m
  df
}))

spearman_nonredund     <- calc_coeffs(nonredund_rules, FALSE)
spearman_nonred_by_con <- calc_coeffs(nonredund_rules, TRUE)




####### Save Results #######

# Write mined rules for reproducibility and inspection
save_dir <- "results"

if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

for (con in consequents) {
  
  con_label <- gsub("=", "_", con)
  con_dir <- file.path(save_dir, con_label)
  
  if (!dir.exists(con_dir)) dir.create(con_dir, recursive = TRUE)
  
  # All mined rules
  write.csv(
    as(rules[[con]]$raw, "data.frame"),
    file.path(con_dir, "mined.csv"),
    row.names = FALSE
  )
  
  # Non-redundant rules with all p-values
  write.csv(
    as(rules[[con]]$pruned, "data.frame"),
    file.path(con_dir, "nonredundant.csv"),
    row.names = FALSE
  )
  
  # Significance-filtered subsets
  for (method in c("uncorr", "BH", "BF")) {
    write.csv(
      as(rules[[con]][[method]], "data.frame"),
      file.path(con_dir, paste0("sig_", method, ".csv")),
      row.names = FALSE
    )
  }
}

####### Save Correlation Results #######

corr_dir <- file.path(save_dir, "correlations")
if (!dir.exists(corr_dir)) dir.create(corr_dir, recursive = TRUE)

# Spearman correlations overall
write.csv(
  spearman_nonredund,
  file.path(corr_dir, "spearman_overall.csv"),
  row.names = FALSE
)

# Spearman correlations by consequent
write.csv(
  spearman_nonred_by_con,
  file.path(corr_dir, "spearman_by_consequent.csv"),
  row.names = FALSE
)

