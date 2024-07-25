# Description of script files

| File name                                                     | Description                                               |
|---------------------------------------------------------------|-----------------------------------------------------------|
| [00_DataPrep.Rmd](./scripts/00_DataPrep.Rmd)                  | UKB data preparation                                      |
| [01_Clustering.Rmd](./scripts/01_Clustering.Rmd)              | Clustering analysis applied in all cohorts                |
| [01_cross_sectional_FX.R](./scripts/01_cross_sectional_FX.R)  | Functions used in `01_Clustering.Rmd`                     |
| [02_ValClus.Rmd](./scripts/02_ValClus.Rmd)                    | Profile validation                                        |
| [03_CentBased.Rmd](./scripts/03_CentBased.Rmd)                | Centroid-based clustering analysis                        |
| [04_BoundaryBased.Rmd](./scripts/04_BoundaryBased.Rmd)        | Boundary-based clustering analysis                        |
| [05_DensityBased.R](./scripts/05_DensityBased.R)              | Density-based clustering analysis                         |
| [06_ValClusOutcomes.Rmd](./scripts/06_ValClusOutcomes.Rmd)    | Profile-outcome associations applied in all cohorts       |
| [06_cross_sectional_FX.R](./scripts/06_cross_sectional_FX.R)  | Functions used in `06_ValClusOutcomes.Rmd`                |
| [07_VisualClus.Rmd](./scripts/07_VisualClus.Rmd)              | Visualisation of concordant and discordant profiles       |
| [08_PrevDx.Rmd](./scripts/08_PrevDx.Rmd)                      | Gathering profile-outcome cross-sectional associations    |
| [09_NoDx.Rmd](./scripts/09_NoDx.Rmd)                          | Distribution of profiles in individuals without diagnosis |
| [10_IncDx.Rmd](./scripts/10_IncDx.Rmd)                        | Gathering profile-outcome longitudinal associations       |
| [11_DataEthnic.Rmd](./scripts/11_DataEthnic.Rmd)              | UKB data preparation for ethnicities other than European  |
| [12_EthnicAssoc.Rmd](./scripts/12_EthnicAssoc.Rmd)            | Ethnic-specific profile-outcome associations              |
| [13_LDL.Rmd](./scripts/13_LDL.Rmd)                            | Net benefit of LDL-cholesterol for benchmarking           |
| [14_Rebuttal.Rmd](./scripts/14_Rebuttal.Rmd)                  | Additional calculations for rebuttal                      |
| [15_LassoCox.Rmd](./scripts/15_LassoCox.Rmd)                  | Applying Lasso penalisation to assess overfitting         |
| [15_LassoFx.R](./scripts/15_LassoFx.R)                        | Functions used in `15_LassoCox.Rmd`                       |
| [16_PrevDxAdj.Rmd](./scripts/16_PrevDxAdj.Rmd)                | Adjusted profile-outcome cross-sectional associations     |