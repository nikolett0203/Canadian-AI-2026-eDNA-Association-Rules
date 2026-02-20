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
| `RulesTools` | Discretization utilities and the BrookTrout dataset |
| `tidyverse` | Data manipulation |
| `ggplot2` | Plotting |
| `patchwork` | Figure composition |
