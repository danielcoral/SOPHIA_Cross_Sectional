get_cluster_descriptives <- function(dataset, clus_probs){
  dataset <- data.frame(dataset)
  long_summary <- data.frame(cluster = integer(0), n_allocated = integer(0), perc_allocated = numeric(0), variable = character(0),  wtd_mean = numeric(0) , wtd_std = numeric(0), wtd_q0 = numeric(0), wtd_q25 = numeric(0), wtd_q50 = numeric(0), wtd_q75 = numeric(0), wtd_q100 = numeric(0))
  for (col in 1:ncol(clus_probs)){
    clus_name <- colnames(clus_probs)[col]
    probs_col <- clus_probs[,col]
    n_allocated <- sum(probs_col > 0.8)
    perc_allocated <- n_allocated/nrow(dataset)
    prob_sum <- sum(probs_col)
    for (var in 2:ncol(dataset)){
      var_name <- colnames(dataset)[var]
      var_col <- dataset[,var]
      wtd_mean <- sum(probs_col*var_col)/prob_sum
      wtd_std <- sum(probs_col*(var_col - wtd_mean))/prob_sum %>% sqrt()
      var_prob <- data.frame(var_col, probs_col)
      var_prob <- var_prob[order(var_prob$var_col),]
      var_prob$probs_cum <- cumsum(var_prob$probs_col)
      q0 <- var_prob$var_col[1]
      i25 <- which.min(abs(var_prob$probs_cum - 0.25*prob_sum))
      q25 <- var_prob$var_col[i25]
      i50 <- which.min(abs(var_prob$probs_cum - 0.5*prob_sum))
      q50 <- var_prob$var_col[i50]
      i75 <- which.min(abs(var_prob$probs_cum - 0.75*prob_sum))
      q75 <- var_prob$var_col[i75]
      q100 <- max(var_col)
      row <- nrow(long_summary) + 1
      long_summary[row,] <- c(clus_name, n_allocated, perc_allocated, var_name, wtd_mean, wtd_std, q0, q25, q50, q75, q100)
    }
  }
  long_summary
}

# use the function like this: cluster_descriptives <- map2(strat_dat, cluster_probs, get_cluster_descriptives)
# works for both continuous data and binary data
# doesn't give prevalence per 1000, should be easily adaptable
# doesn't require Hmisc anymore
