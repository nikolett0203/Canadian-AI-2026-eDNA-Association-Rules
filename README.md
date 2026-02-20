# Association Rule Mining for eDNA Detection Datasets

This repository contains the R script accompanying the paper:

> *To correct or not to correct?: Assessing the multiple comparisons problem for association rule mining of environmental DNA (eDNA) detection survey datasets*  
> Canadian AI 2026

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
BrookTrout <- read.csv("BrookTrout.csv")
```

The dataset comprises 126 transactions across 10 variables (eDNA concentrations, 
electrofishing counts, and eight physicochemical metadata variables) with no missing values.

