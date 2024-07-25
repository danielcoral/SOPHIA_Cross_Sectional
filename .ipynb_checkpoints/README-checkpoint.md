# Profiles of phenotypic discordance for a given BMI

Overall analysis pipeline for discovery and replication:

![](./plots/AnalysisPlot.png)

How we used UMAP for clustering:

![](./plots/Clusmethod.png)

## Contents

An overall description of the files in the [scripts](./scripts) folder, which we used to run our analysis is found [here](./scripts_description.md).

## System requirements

We run the analysis in the following R environment and package versions:

```{r}
R version 4.3.3 (2024-02-29)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 9.2 (Plow)

Matrix products: default
BLAS/LAPACK: /gpfs/gpfs0/Home/daniel_c/miniforge3/envs/RLang/lib/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

time zone: Europe/Stockholm
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] glmnet_4.1-8           ggdensity_1.0.0        GGally_2.2.1
 [4] ClustOfVar_1.1         kernlab_0.9-32         Rtsne_0.17
 [7] dcurves_0.5.0          MGMM_1.0.1.1           broom.mixed_0.2.9.4
[10] broom_1.0.6            meta_7.0-0             metadat_1.2-0
[13] lmerTest_3.1-3         lme4_1.1-35.3          scales_1.3.0
[16] reshape2_1.4.4         ggflowchart_1.0.0.9007 survival_3.7-0
[19] dbscan_1.1-12          archetypes_2.2-0.1     nnls_1.5
[22] modeltools_0.2-23      mclust_6.0.1           ggraph_2.1.0.9000
[25] patchwork_1.2.0        ggh4x_0.2.8            ggplot2_3.5.0
[28] stringr_1.5.1          rio_1.0.1              mvtnorm_1.2-4
[31] igraph_2.0.3           uwot_0.1.16            Matrix_1.6-5
[34] lubridate_1.9.3        purrr_1.0.2            tidyr_1.3.1
[37] dplyr_1.1.4            readr_2.1.5

loaded via a namespace (and not attached):
 [1] gridExtra_2.3       rlang_1.1.3         magrittr_2.0.3
 [4] furrr_0.3.1         compiler_4.3.3      vctrs_0.6.5
 [7] shape_1.4.6.1       pkgconfig_2.0.3     fastmap_1.2.0
[10] backports_1.4.1     utf8_1.2.4          tzdb_0.4.0
[13] nloptr_2.0.3        cachem_1.1.0        tweenr_2.0.2
[16] parallel_4.3.3      R6_2.5.1            stringi_1.8.3
[19] RColorBrewer_1.1-3  parallelly_1.37.1   boot_1.3-30
[22] numDeriv_2016.8-1.1 iterators_1.0.14    Rcpp_1.0.12
[25] R.utils_2.12.3      splines_4.3.3       timechange_0.3.0
[28] tidyselect_1.2.0    viridis_0.6.5       codetools_0.2-20
[31] metafor_4.4-0       listenv_0.9.1       lattice_0.22-5
[34] tibble_3.2.1        plyr_1.8.9          withr_3.0.0
[37] future_1.33.2       CompQuadForm_1.4.3  ggstats_0.6.0
[40] polyclip_1.10-6     xml2_1.3.6          pillar_1.9.0
[43] foreach_1.5.2       generics_0.1.3      mathjaxr_1.6-0
[46] hms_1.1.3           munsell_0.5.0       minqa_1.2.7
[49] globals_0.16.3      glue_1.7.0          tools_4.3.3
[52] forcats_1.0.0       graphlayouts_1.1.0  tidygraph_1.3.0
[55] grid_4.3.3          colorspace_2.1-0    nlme_3.1-164
[58] ggforce_0.4.1       cli_3.6.2           fansi_1.0.6
[61] viridisLite_0.4.2   gtable_0.3.4        R.methodsS3_1.8.2
[64] digest_0.6.35       ggrepel_0.9.5       farver_2.1.1
[67] memoise_2.0.1       R.oo_1.26.0         lifecycle_1.0.4
[70] MASS_7.3-60
```

## Installation guide

To install R, please follow the instructions found [here](https://www.r-project.org/).

To install the R packages needed for our analysis you can run the following:

```{r}
install.packages(
    c(
      "readr", "dplyr", "tidyr", "purrr", "lubridate",
      "uwot", "igraph", "mvtnorm", "rio", "stringr", 
      "ggh4x", "patchwork", "ggraph", "mclust", 
      "archetypes", "dbscan", "survival", "ggflowchart", 
      "reshape2", "scales", "lme4", "lmerTest", "meta", 
      "broom", "broom.mixed", "MGMM", "dcurves", "Rtsne", 
      "kernlab", "ClustOfVar", "GGally", "ggdensity", "glmnet"
    )
)
```

On a typical desktop computer this takes around 1 hour.

## Demo

A demo showing how to run our clustering analysis on a simulated dataset can be found [here](./demo/Demo.ipynb). This uses the same functions we used in our analysis. The simulated dataset is derived from the parameters of the profiles we identified, and will have the same format that we used in our original analysis. The parameters of the profiles are stored in the R object [validclusmod.RData](./data/validclusmod.RData). We also show in the demo how to calculate the probabilities to have the profiles we identified for a new set of individuals, using the simulated data as an example. This is useful to map new cohorts to these profiles.

## Reproduction

Due to data access restrictions, we cannot include the data necessary to reproduce all the quantitative results in the manuscript. However, it is possible to run the same analysis in a new dataset provided it has the same format as we used. The simulated dataset from the demo shows how this should look like. Additionally, see [01_Clustering.Rmd](./scripts/01_Clustering.Rmd) and [06_ValClusOutcomes.Rmd](./scripts/06_ValClusOutcomes.Rmd) for a detailed explanation of the formats expected and the functions used.

---
