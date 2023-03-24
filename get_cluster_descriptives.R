get_cluster_descriptives <- function(dataset, cluster_probs){
  long_summary <- data.frame(cluster = integer(0), variable = character(0), n_allocated = integer(0), perc_allocated = numeric(0), wtd_mean = numeric(0), wtd_std = numeric(0), wtd_q0 = numeric(0), wtd_q25 = numeric(0), wtd_q50 = numeric(0), wtd_q75 = numeric(0), wtd_q100 = numeric(0))
  for (col in 1:ncol(cluster_probs)){
    clus <- colnames(cluster_probs)[col]
    n_allocated <- sum(cluster_probs[,col] > 0.8)
    perc_allocated <- n_allocated/nrow(dataset)
    for (var in 1:ncol(dataset)){
      varia <- colnames(dataset)[var]
      wtd_mean <- wtd.mean(dataset[,var], cluster_probs[,col])
      wtd_std <- wtd.var(dataset[,var], cluster_probs[,col]) %>% sqrt()
      wtd_iqr <- wtd.quantile(dataset[,var], cluster_probs[,col])
      row <- nrow(long_summary) + 1
      long_summary[row,] <- c(clus, varia, n_allocated, perc_allocated, wtd_mean, wtd_std, wtd_iqr)
    }
  }
  long_summary
}

# use the function like this: cluster_descriptives <- map2(strat_dat, cluster_probs, get_cluster_descriptives)
# requires Hmisc package
