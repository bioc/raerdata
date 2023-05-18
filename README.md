# raerdata

raerdata provides external  data sets for calling rna-editing sites with [raer](https://rnabioco.github.io/raer).

## Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("raerdata")
```

## To use `raerdata`

``` r
library(ExperimentHub)
eh <- ExperimentHub()

## query
refs <- query(eh, "raerdata")
refs

```
