# Association Rule Mining for eDNA Detection Datasets

This repository contains the R script accompanying the paper:

> *To correct or not to correct?: Assessing the multiple comparisons problem for association rule mining of environmental DNA (eDNA) detection survey datasets*  

---

## Overview

In this script, we apply association rule mining (ARM) to a published brook trout eDNA dataset to (1) demonstrate ARM's utility as an exploratory screening tool for small environmental datasets and (2) evaluate the effect of Bonferroni and Benjamini–Hochberg multiple testing corrections on the quality and statistical reliability of mined rules.

---

## Requirements
### R Packages

| Package | Purpose |
|---|---|
| `arules` | Rule mining, redundancy pruning, significance testing |
| `tidyverse` | Data manipulation |
| `ggplot2` | Plotting |
| `patchwork` | Figure composition |

Install from CRAN:
```r
install.packages(c("arules", "tidyverse", "ggplot2", "patchwork"))
```

### Data

This analysis uses `BrookTrout.csv`, included in this repository. The data are derived from:

> Nolan et al. (2023). *Detection of native brook trout (Salvelinus fontinalis) 
> using environmental DNA and electrofishing.*

It is loaded in the script with:
```r
BrookTrout <- read.csv("data/BrookTrout.csv", stringsAsFactors = FALSE)
```

The dataset comprises 126 transactions across 10 variables (eDNA concentrations, 
electrofishing counts, and eight physicochemical metadata variables) with no missing values.

## Usage

Run the script end-to-end in a single R session:

```r
source("arm_analysis.R")
```

The script executes the following steps in order:
1. **Data pre-processing** — discretizes continuous variables into categorical form
2. **Rule mining** — mines association rules using the Apriori algorithm via `arules`
3. **Redundancy pruning** — removes structurally redundant rules
4. **Significance testing** — filters rules at α = 0.05 using Fisher's Exact Tests
5. **Multiple testing correction** — applies Bonferroni and Benjamini–Hochberg corrections
6. **Scatterplots** — generates figures of rule metrics vs. −log₁₀(p)
7. **Correlation analysis** — computes Spearman correlations between p-values and interestingness measures
8. **Save results** — writes all rule sets and correlation statistics to `results/`


### Output Structure

```
results/
├── eDNAConc_high/
│   ├── mined.csv              # All mined rules (unfiltered)
│   ├── nonredundant.csv       # Non-redundant rules with all p-values
│   ├── sig_uncorr.csv         # Significant rules, no correction
│   ├── sig_BH.csv             # Significant rules, Benjamini-Hochberg
│   └── sig_BF.csv             # Significant rules, Bonferroni
├── eFishCatch_present/
│   ├── mined.csv
│   ├── nonredundant.csv
│   ├── sig_uncorr.csv
│   ├── sig_BH.csv
│   └── sig_BF.csv
└── correlations/
    ├── spearman_overall.csv        # Spearman correlations pooled across consequents
    └── spearman_by_consequent.csv  # Spearman correlations by consequent
```

