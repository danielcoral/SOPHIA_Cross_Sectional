# Profiles of phenotypic discordance for a given BMI

Overall analysis pipeline for discovery and replication:

![](./plots/AnalysisPlot.png)

How we used UMAP for clustering:

![](./plots/Clusmethod.png)

An overall description of the files in the `scripts` folder, which we used to run our analysis is found [here](./guides/script_description.md)

## System requirements

We run the analysis in the following R environment:

```{r}
Hello
```

## Installation guide

To install R follow the instructions found [here](https://www.r-project.org/).

To install the R packages needed for our analysis you can run:

```{r}
Hello
```

On a typical desktop computer this takes around 

## Demo

A demo showing how to run our clustering analysis on a simulated dataset can be found [here](./demo/demo.ipynb). This uses the same functions we used in our analysis. The simulated dataset is derived from the parameters of the profiles we identified, and will have the same format that we used in our original analysis. The parameters of the profiles are stored in the R object [validclusmod.RData](./data/validclusmod.RData).

## Reproduction

Due to data access restrictions, we cannot include the data necessary to reproduce all the quantitative results in the manuscript. However, it is possible to run the same analysis in a new dataset provided it has the same format as we used. See [01_Clustering.Rmd](./scripts/01_Clustering.Rmd) and [06_ValClusOutcomes.Rmd](./scripts/06_ValClusOutcomes.Rmd) for a detailed explanation of the formats expected and the functions used.

---
