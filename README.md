# CIFU NRR Curve Analysis

R scripts for processing, aligning, and comparing OCT and fundus neuroretinal rim (NRR) curves.

## Files

- `CIFU_OCT_vs_Fundus_NRR.R` — cleans and normalizes NRR data, performs elastic curve alignment and fusion, and clusters fused curves using `funFEM`.
- `landmark_registration_nonglaucoma.R` — aligns fundus curves to OCT curves using landmark registration and creates fused curves.
- `Pairwise_functional_cca.R` — compares OCT, fundus, original fused, and landmark-fused curves using functional canonical correlation analysis.

## Requirements

Main R packages:

```r
install.packages(c("fda", "fdasrvf", "funFEM", "pracma", "tidyverse"))
