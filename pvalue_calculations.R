####### Libraries #######

library(RulesTools)
library(arules)
library(tidyverse)
library(ggplot2)
library(patchwork)



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
  13.3,                      # For eDNA concentrations (50% LOD)
  0,                         # For electrofishing catch (absence/presence)
  15.03,                     # For air temp (avg temp for sample site in Sept)
  15,                        # For water temp (optimal temp for brook trout)
  7.75,                      # For pH (CCME guideline)
  9.5,                       # For dissolved oxygen (CCME guideline)
  1047,                      # For conductivity (median)
  1.04                       # For volume filtered (median)
)

consequents <- c(
  "eDNAConc=high", 
  "eFishCatch=present" 
)




####### Helper Functions #######

get_p <- function(rules, transactions) {
  
  if (length(rules) == 0) 
    return(numeric(0))
  
  quality_df <- quality(rules)
  n <- length(transactions)                       # could also use NUM_TRANSACTIONS but this is more general
  
  lhs <- unique(lhs(rules))
  rhs <- unique(rhs(rules))
  
  lhs_supp <- support(lhs, transactions) * n      # support() returns fraction, multiply by n to get count
  rhs_supp <- support(rhs, transactions) * n
  
  # create support lookup table
  lhs_labels <- labels(lhs)
  rhs_labels <- labels(rhs)
  names(lhs_supp) <- lhs_labels
  names(rhs_supp) <- rhs_labels
  
  # vectorized p-value calcs
  pvalues <- vapply(1:length(rules), function(i) {
    
    # get count of antecedent&consequent
    a <- quality_df$count[i]
    
    # lookup count of antecedent
    lhs_label <- labels(lhs(rules[i]))
    b <- lhs_supp[lhs_label] - a
    
    # lookup count of consequent
    rhs_label <- labels(rhs(rules[i]))
    c <- rhs_supp[rhs_label] - a
    
    d <- n - a - b - c
    
    fisher.test(matrix(c(a, b, c, d), nrow = 2), 
                alternative = "greater")$p.value
    
  }, FUN.VALUE = numeric(1))
  
  return(pvalues)

}

sig_sort <- function(rules, method) {
  sort(rules[is.significant(rules, alpha = 0.05, adjust = method)],
       by = c("lift", "confidence", "support"))
}

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
discretized_df[["Site"]]     <- as.factor(BrookTrout$Site)

# convert df into transactions
transactions <- as(discretized_df, "transactions")




####### Rule Mining #######

rules <- list()

for (con in consequents) {
  
  # mine rules for each target: eDNA=high and eFishCatch=present
  raw_rules <- apriori(
    transactions,
    parameter = list(support = 1/NUM_TRANSACTIONS, confidence = 1/NUM_TRANSACTIONS),
    appearance = list(rhs = con)
  )
  
  # subset non-redundant rules
  pruned_rules <- raw_rules[!is.redundant(raw_rules, measure = "confidence")]
  
  # isolate ps
  p_values <- get_p(pruned_rules, transactions)
  
  # perform adjustments
  quality(pruned_rules)$p_raw   <- p_values
  quality(pruned_rules)$p_BH    <- p.adjust(p_values, method = "BH")
  quality(pruned_rules)$p_BF    <- p.adjust(p_values, method = "bonferroni")
  quality(pruned_rules)$len     <- size(pruned_rules)
  quality(pruned_rules)$lhs     <- labels(lhs(pruned_rules))
  
  # final rules list
  rules[[con]] <- list(
    raw    = raw_rules,
    pruned = pruned_rules,
    uncorr  = sig_sort(pruned_rules, "none"),
    BH     = sig_sort(pruned_rules, "BH"),
    BF     = sig_sort(pruned_rules, "bonferroni")
  )
  
}




####### Scatterplots #######

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

methods <- c("uncorr", "BH", "BF")

nonredund_rules <- do.call(rbind, lapply(methods, function(m) {
  df <- prep_for_stats(rules, m)
  df$method <- m
  df
}))

spearman_nonredund     <- calc_coeffs(nonredund_rules, FALSE)
spearman_nonred_by_con <- calc_coeffs(nonredund_rules, TRUE)




####### Save Results #######

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




