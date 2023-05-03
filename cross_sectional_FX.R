## Functions to run clustering analysis to detect subgroups with different relationships between BMI and risk

# Declaring variables that are expected to be found in the data to run analysis
IDCOL <- "eid"
SEXVAR <- "sex"
SEXGROUPS <- c("Female", "Male")
BODYSIZEINDEX <- "bmi"
COVARIATES <- c("age", "smoking")
BIOMARKERS <- c("whr", "sbp", "dbp", "alt", "scr", "crp", "hdl", "tg", "ldl", "fg")
PERSEXVARS <- c(BODYSIZEINDEX, COVARIATES, BIOMARKERS)
ALLVARS <- c(IDCOL, SEXVAR, PERSEXVARS)


# Outlier removal
remove_outliers <- function(x, sdunits = 5){
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    upperbound <- m + (sdunits*s)
    lowerbound <- m - (sdunits*s)
    ifelse((x > lowerbound) & (x < upperbound), x, NaN)
}

# Checking data input
checkinput <- function(X){
    stopifnot(is.data.frame(X))
    stopifnot(ALLVARS %in% names(X))
    stopifnot(names(X) %in% ALLVARS)
}

# Variable types
vartype <- function(x){
    if(is.numeric(x)){
        if(length(unique(x)) > 10){
            return("Numeric")
        } else { 
            return("Categorical")
        }
    } else {
        if(is.character(x) | is.factor(x)){
            return("Categorical")
        } else {
            return(NA)
        }
    }
}

# Summary of categorical variables
categorical_sumstats <- function(x){
    N <- length(x)
    N_miss <- sum(is.na(x)) 
    props <- table(x)
    props <- sort(props, decreasing = TRUE)
    props <- data.frame(N = N, N_miss = N_miss,
                        Summary1 = names(props), 
                        Summary2 = as.vector(props))
    props$Summary2 <- paste0(props$Summary2, " (", round(100 * props$Summary2 / sum(props$Summary2), 2), "%)")
    props
}

# Summary of continuous variables
continuous_sumstats <- function(x){
    N <- length(x)
    N_miss <- sum(is.na(x))
    Mn <- round(mean(x, na.rm = TRUE), 2)
    Md <- round(median(x, na.rm = TRUE), 2)
    Sd <- round(sd(x, na.rm = TRUE), 2)
    p25_75 <- round(quantile(x, probs = c(.25, .75), na.rm = TRUE), 2)
    p25_75 <- paste(p25_75, collapse = " - ")
    data.frame(N = N, N_miss = N_miss,
               Summary1 = paste0(Mn, " (", Sd, ")"),
               Summary2 = paste0(Md, " (", p25_75, ")"))
}

# Pick the summary function given variable type
sumstatfx <- function(x){
    vartyp <- vartype(x)
    xres <- switch(
        vartyp,
        Categorical = categorical_sumstats(x),
        Numeric = continuous_sumstats(x),
        stop("A column is not in the right format")
    )
    data.frame(Type = vartyp, xres)
}

# Summary of all variables
gendescfx <- function(X){
    checkinput(X)
    purrr::map_dfr(
        setNames(PERSEXVARS, PERSEXVARS),
        function(VRBL){ sumstatfx(X[[VRBL]]) },
        .id = "Variable"
    )
}

# General descriptives for each sex
get_general_descriptives <- function(sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ gendescfx(sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Generate models for all biomarkers
bmimodfx <- function(X){
    checkinput(X)
    mdtb <- tibble::tibble(Biomarker = BIOMARKERS)
    dplyr::mutate(
        mdtb,
        mod = purrr::map(
            Biomarker,
            function(biomarker){
                lm(reformulate(response = biomarker, termlabels = c(BODYSIZEINDEX, COVARIATES)), data = X, na.action = na.exclude)
            }
        )
    )
}

# Generate models by sex
get_bmimods <- function(sexstrats){
    purrr::map_dfr(setNames(SEXGROUPS, SEXGROUPS), 
                   function(sexgr){ bmimodfx(sexstrats[[sexgr]]) },
                   .id = "sex")
}

# Extract BMI effect on a biomarker
bmieff <- function(MOD){
    modsum <- summary(MOD)
    modcoefs <- modsum$coefficients["bmi", 1:2, drop = FALSE]
    modCI <- confint(MOD)["bmi", 1:2, drop = FALSE]
    modcoefs <- cbind(modcoefs, modCI)
    modcoefs <- round(modcoefs, digits = 5)
    colnames(modcoefs) <- c("Estimate", "SE", "lowerCI", "upperCI")
    data.frame(modcoefs)
}

# Extract BMI effects on all biomarkers in each sex
get_bmicoefs <- function(mdtb){
    mdtb$mod <- purrr::map(mdtb$mod, bmieff)
    tidyr::unnest(mdtb, mod)
}

# Generate table of residuals for each sex
get_residtabs <- function(sexstrats, mdtb){
    mdtb$mod <- purrr::map(mdtb$mod, resid)
    mdtb$mod <- purrr::map(mdtb$mod, scale)
    mdtb$mod <- purrr::map(mdtb$mod, as.vector)
    mdtb <- tidyr::pivot_wider(mdtb, names_from = Biomarker, values_from = mod)
    purrr::map(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){
            residtb <- tidyr::unnest(mdtb[mdtb$sex == sexgr,], -sex)
            residtb$sex <- NULL
            tibble::tibble(eid = sexstrats[[sexgr]]$eid, residtb)
        }
    )
}

# Fitting a Gaussian mixture model for given centers and covariance matrices
gmmfx <- function(DATA, CLUS){
    mus <- lapply(CLUS, "[[", "center")
    covmats <- lapply(CLUS, "[[", "cov")
    pdfs <- mapply(function(mu, covmat){ mvtnorm::dmvnorm(DATA, mu, covmat) }, mus, covmats)
    k <- length(mus)
    ws <- rep(1 / k, k)
    L <- pdfs %*% diag(ws)
    probs <- L / rowSums(L)
    Lclus <- colSums(probs)
    loglik <- sum(log(rowSums(L)))
    Delta <- 1
    iter <- 1
    itermax <- 100
    while(Delta > 1e-10 && iter <= itermax){
        ws <- Lclus / sum(Lclus)
        L <- pdfs %*% diag(ws)
        probs <- L / rowSums(L)
        Lclus <- colSums(probs)
        loglik_current <- sum(log(rowSums(L)))
        Delta <- loglik_current - loglik
        loglik <- loglik_current
        iter <- iter + 1
    }
    colnames(probs) <- names(mus)
    probs <- data.frame(probs)
    names(ws) <- names(mus)
    message("\t\tConvergence reached in ", iter, " iterations.")
    list(probs = probs, weights = ws)
}

# Function to find clusters in residual data from BMI
umapclus_fx <- function(RESIDMAT){
    message("\t1. Setting parameters...")
    eids <- RESIDMAT$eid
    Xmat <- RESIDMAT[,names(RESIDMAT) != "eid"]
    n_total <- nrow(Xmat)
    nn <- max(10, round(10 + 15 * (log10(n_total) - 4)))
    message("\t2. Running UMAP...")
    umapres <- uwot::umap(
        Xmat, n_components = 2, 
        n_neighbors = nn, nn_method = "annoy", n_trees = 100, n_sgd_threads = "auto", 
        init = "pca", n_epochs = 500, approx_pow = TRUE,
        binary_edge_weights = TRUE, dens_scale = 1, 
        ret_extra = "fgraph", verbose = FALSE       
    )
    message("\t3. Extracting graph...")
    umap_graph <- igraph::graph_from_adjacency_matrix(umapres$fgraph, mode = "undirected")
    igraph::V(umap_graph)$name <- eids
    message("\t4. Initializing graph clustering using the leading eigen vector method...")
    initgraphclus <- igraph::cluster_leading_eigen(umap_graph)
    message("\t5. Optimizing graph clustering using the Leiden algorithm...")
    graphclus <- igraph::cluster_leiden(umap_graph, objective_function = "modularity", 
                                        initial_membership = initgraphclus$membership, n_iterations = 500)
    Nclus <- length(igraph::sizes(graphclus))
    Qual <- round(graphclus$quality, 2)
    message("\t   ", Nclus, " clusters found. Modularity = ", Qual, ".")
    message("\t6. Extracting cluster membership...")
    member_idx <- igraph::membership(graphclus)
    member_idx <- data.frame(eid = as.numeric(names(member_idx)), cluster = as.numeric(member_idx))
    member_idx <- split(member_idx, member_idx$cluster)
    message("\t7. Calculating eigen centrality in clusters...")
    cluster_subgraphs <- purrr::map(member_idx, function(clus){ igraph::induced_subgraph(umap_graph, as.character(clus$eid)) })
    eigencent <- purrr::map_dfr(
        cluster_subgraphs,
        function(clus){
            res <- igraph::eigen_centrality(clus, directed = FALSE)$vector
            data.frame(
                eid = as.numeric(names(res)),
                value = as.numeric(res)
            )
        },
        .id = "cluster"
    )
    message("\t8. Calculating Gaussian subdistributions weighted by eigen centralities...")
    eigencent <- dplyr::inner_join(eigencent, RESIDMAT, by = "eid")
    eigencent <- split(eigencent, eigencent$cluster)
    cluster_pars <- purrr::map(
        eigencent, 
        function(clus){
            cluseigen <- clus$value
            clusmat <- clus[,!(names(clus) %in% c("eid", "cluster", "value"))]
            cov.wt(clusmat, wt = cluseigen)
        }
    )
    message("\t9. Adding central/concordant subdistribution...")
    clus0_mu <- rep(0, ncol(Xmat))
    names(clus0_mu) <- colnames(Xmat)
    clus0_cov <- diag(rep(1, ncol(Xmat)))
    dimnames(clus0_cov) <- list(colnames(Xmat), colnames(Xmat))
    clus0 <- list(
        cov = clus0_cov,
        center = clus0_mu
    )
    cluster_pars[["0"]] <- clus0
    message("\t10. Fitting a mixture of Gaussians with subdistributions found...")
    gmmres <- gmmfx(Xmat, cluster_pars)
    message("\t11. Organizing results...")
    for(i in seq_along(cluster_pars)){
        cluster_pars[[i]][["weight"]] <- gmmres$weights[[i]]
    }
    names(cluster_pars) <- paste("cluster", names(cluster_pars), sep = "_")
    probdat <- data.frame(eid = eids, gmmres$probs)
    colnames(probdat) <- gsub("X", "cluster_", colnames(probdat))
    embedding <- data.frame(eid = eids, umapres$embedding)
    colnames(embedding) <- gsub("X", "UMAP", colnames(embedding))
    message("Done!")
    return(list(umap_model = umapres, probs = probdat,
                clusters = cluster_pars, modularity = Qual))
}

# Running clustering on sex-specific residual data
get_cluster_results <- function(RESIDMATS){
    purrr::map(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){
            message(paste("Starting analysis for the", sexgr, "subset:"))
            umapclus_fx(RESIDMATS[[sexgr]])
        }
    )
}

# Summary of clusters
clussumfx <- function(clusdat){
    clusnam <- names(clusdat$clusters)
    purrr::map_dfr(
        setNames(clusnam, clusnam),
        function(clus){
            mu <- clusdat$clusters[[clus]][["center"]]
            covmat <- clusdat$clusters[[clus]][["cov"]]
            VAR <- diag(covmat)
            SD <- sqrt(VAR)
            w <- clusdat$clusters[[clus]][["weight"]]
            tibble::tibble(trait = names(mu), center = mu, SD = SD, weight = w)
        },
        .id = "cluster"
    )
}

# Summary of clusters per sex
get_clustersummary <- function(clusres){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS), 
        function(sexgr){ clussumfx(clusres[[sexgr]]) },
        .id = "sex"
    )
}

# Weighted mean
wtd_mean <- function(x, p){
    sum(x*p)/sum(p)
}

# Weighted variance
wtd_var <- function(x, p){
    xhat <- wtd_mean(x, p)
    sum(((x - xhat)^2)*p)/sum(p)
}

# Weighted quantiles
# Taken from: https://aakinshin.net/posts/weighted-quantiles/
wtd_quantile <- function(x, p, q){
    n <- sum(p)
    indexes <- order(x)
    x <- x[indexes]
    p <- p[indexes]
    p <- p / n
    cdfprobs <- cumsum(c(0, p))
    sapply(
        setNames(q, q), 
        function(QUANTILE) {
            Q <- pbeta(cdfprobs, n * QUANTILE, n * (1 - QUANTILE))
            W <- tail(Q, -1) - head(Q, -1)
            sum(W * x)
        }
    )
}

# Weighted proportions
wtd_table <- function(x, p){
    LEVELS <- unique(x)
    sapply(
        setNames(LEVELS, LEVELS),
        function(LEVEL){
            sum((x == LEVEL) * p) / sum(p)
        }
    )
}

# Weighted summary of categorical variables
wtd_categorical_sumstats <- function(x, p){
    props <- wtd_table(x, p)
    props <- sort(props, decreasing = TRUE)
    props <- data.frame(Summary1 = names(props), 
                        Summary2 = as.vector(props))
    props$Summary2 <- paste(round(100 * props$Summary2, 2), "%")
    return(props)
}

# Weighted summary of continuous variables
wtd_continuous_sumstats <- function(x, p){
    Mn <- round(wtd_mean(x, p), 2)
    Sd <- round(sqrt(wtd_var(x, p)), 2)
    quants <- wtd_quantile(x, p, q = c(.25, .5, .75))
    Md <- round(quants[2], 2)
    p25_75 <- round(quants[c(1, 3)], 2)
    p25_75 <- paste(p25_75, collapse = " - ")
    data.frame(Summary1 = paste0(Mn, " (", Sd, ")"),
               Summary2 = paste0(Md, " (", p25_75, ")"))
}

# Weighted summary of a variable
wtd_sumstatfx <- function(x, p){
    vartyp <- vartype(x)
    xres <- switch(
        vartyp,
        Categorical = wtd_categorical_sumstats(x, p),
        Numeric = wtd_continuous_sumstats(x, p),
        stop("A column is not in the right format")
    )
    data.frame(Type = vartyp, xres)
}

# Weighted summary of all variables
wtdgendescfx <- function(X, p){
    purrr::map_dfr(
        setNames(PERSEXVARS, PERSEXVARS),
        function(VRBL){ wtd_sumstatfx(X[[VRBL]], p) },
        .id = "Variable"
    )
}

clusgendescfx <- function(clusdat, X){
    probs <- clusdat$probs
    probs <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    joindf <- dplyr::inner_join(probs, X, by = "eid")
    joindf <- dplyr::group_by(joindf, cluster)
    joindf <- tidyr::nest(joindf)
    joindf <- dplyr::transmute(
        joindf,
        TotalN = purrr::map_dbl(data, nrow),
        Weighted_N = purrr::map_dbl(data, function(d){ round(sum(d$cluster_prob), 2) }),
        N80Perc = purrr::map_dbl(data, function(d){ sum(d$cluster_prob >= .8) }),
        varsum = purrr::map(data, function(d){ wtdgendescfx(d, d$cluster_prob) })
    )
    joindf <- tidyr::unnest(joindf, varsum)
    dplyr::ungroup(joindf)
}

# Cluster descriptives
get_cluster_descriptives <- function(clusres, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ clusgendescfx(clusres[[sexgr]], sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Models of BMI-biomarker association weighted by cluster probability
cluslinmod <- function(X, p){
    purrr::map_dfr(
        setNames(BIOMARKERS, BIOMARKERS),
        function(biomarker){
            wmod <- lm(reformulate(response = biomarker, termlabels = c(BODYSIZEINDEX, COVARIATES)), data = X, weights = p)
            bmieff(wmod)
        },
        .id = "Biomarker"
    )
}
            

# Cluster-specific BMI-biomarker associations
bmieffclusfx <- function(clusdat, X){
    probs <- clusdat$probs
    probs <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    joindf <- dplyr::inner_join(probs, X, by = "eid")
    joindf <- dplyr::group_by(joindf, cluster)
    joindf <- tidyr::nest(joindf)
    joindf <- dplyr::mutate(
        joindf,
        data = purrr::map(data, function(d){ cluslinmod(d, d$cluster_prob) })
    )
    joindf <- tidyr::unnest(joindf, data)
    dplyr::ungroup(joindf)
}

# Cluster-specific BMI-biomarker associations for each sex
get_bmicoefs_clusters <- function(clusres, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ bmieffclusfx(clusres[[sexgr]], sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Disease prevalence
diseasesumfx <- function(diseasedf, X){
    dxdf <- diseasedf[diseasedf$eid %in% X$eid,]
    diseases <- names(dxdf)
    diseases <- diseases[diseases != "eid"]
    purrr::map_dfr(
        setNames(diseases, diseases),
        function(disease){
            dx <- dxdf[,disease, drop = TRUE]
            dxnomiss <- dx[!is.na(dx)]
            if(length(dxnomiss) > 0){
                data.frame(Disease = disease, categorical_sumstats(dxnomiss))
            } else {
                data.frame(NULL)
            }
        },
        .id = "Disease"
    )
}

# Disease prevalence by sex
get_diseasesummary <- function(diseasedat, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ diseasesumfx(diseasedat, sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Disease prevalence by cluster
clusdiseasesumfx <- function(diseasedf, clusdat){
    probs <- clusdat$probs
    dxdf <- diseasedf[diseasedf$eid %in% probs$eid,]
    probs <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    dxdf <- tidyr::pivot_longer(dxdf, -eid, names_to = "disease", values_to = "status", values_drop_na = TRUE)
    joindf <- dplyr::inner_join(probs, dxdf, by = "eid", relationship = "many-to-many")
    joindf <- dplyr::group_by(joindf, cluster, disease)
    joindf <- tidyr::nest(joindf)
    joindf <- dplyr::transmute(
        joindf,
        varsum = purrr::map(
            data, 
            function(d){ wtd_categorical_sumstats(d$status, d$cluster_prob) }
        )
    )
    joindf <- tidyr::unnest(joindf, varsum)
    dplyr::ungroup(joindf)
}

# Disease prevalence by cluster in each sex
get_diseasesummary_clusters <- function(diseasedat, clusres){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ clusdiseasesumfx(diseasedat, clusres[[sexgr]]) },
        .id = "sex"
    )
}

# Converting probabilities to log probabilities
ptologp <- function(p){
    p <- ifelse(p < .Machine$double.eps, .Machine$double.eps, p)
    p <- ifelse(p > .9999999999999999, .9999999999999999, p)
    log(p)
}

# Extract association estimates of cluster allocation probabilities to a disease
extractclusglmcoefs <- function(mod){
    modsum <- summary(mod)
    modcoefs <- modsum$coefficients
    modcoefs <- modcoefs[grepl("cluster", rownames(modcoefs)), 1:2, drop = FALSE]
    modCI <- confint.default(mod)
    modCI <- modCI[grepl("cluster", rownames(modCI)), 1:2, drop = FALSE]
    modcoefs <- cbind(modcoefs, modCI)
    modcoefs <- round(modcoefs, digits = 5)
    modcoefs <- cbind(rownames(modcoefs), modcoefs)
    colnames(modcoefs) <- c("cluster", "Estimate", "SE", "lowerCI", "upperCI")
    data.frame(modcoefs)
}

# Models of association between cluster allocation probabilities and a disease
adjclusdiseasesumfx <- function(clusdat, diseasedf, medsdf, X){
    probs <- clusdat$probs
    clusnams <- names(probs)
    clusnams <- clusnams[!(clusnams %in% c("eid", "cluster_0"))]
    joindf <- probs[,c("eid", clusnams)]
    for(clus in clusnams){
        joindf[[clus]] <- ptologp(joindf[[clus]])
    }
    basecovardf <- X[,c("eid", BODYSIZEINDEX, COVARIATES)]
    joindf <- dplyr::inner_join(joindf, basecovardf, by = "eid")
    meds <- names(medsdf)
    meds <- meds[meds != "eid"]
    joindf <- dplyr::inner_join(joindf, medsdf, by = "eid")
    dxdf <- tidyr::pivot_longer(diseasedf, -eid, names_to = "disease", values_to = "status", values_drop_na = TRUE)
    joindf <- dplyr::inner_join(joindf, dxdf, by = "eid")
    joindf <- tidyr::drop_na(joindf)
    joindf <- dplyr::group_by(joindf, disease)
    joindf <- tidyr::nest(joindf)
    joindf <- dplyr::mutate(
        joindf,
        data = purrr::map(
            data,
            function(d){
                mod1 <- glm(reformulate(response = "status", termlabels = clusnams), data = d, family = "binomial")
                mod2 <- glm(reformulate(response = "status", termlabels = c(clusnams, BODYSIZEINDEX, COVARIATES, meds)), data = d, family = "binomial")
                dplyr::bind_rows(list(Unadjusted = extractclusglmcoefs(mod1), Adjusted = extractclusglmcoefs(mod2)), .id = "model")
            }
        )
    )
    joindf <- tidyr::unnest(joindf, data)
    dplyr::ungroup(joindf)
}

# Models of association between cluster allocation probabilities and a disease by sex
get_adjdiseasesummary_clusters <- function(clusres, diseasedat, medsdat, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ adjclusdiseasesumfx(clusres[[sexgr]], diseasedat, medsdat, sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Models of association between cluster labels given at a probability threshold and a disease
adjclusdiseasesumfx_thresh <- function(clusdat, diseasedf, medsdf, X, thresh){
    probs <- clusdat$probs
    joindf <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    basecovardf <- X[,c("eid", BODYSIZEINDEX, COVARIATES)]
    joindf <- dplyr::inner_join(joindf, basecovardf, by = "eid")
    meds <- names(medsdf)
    meds <- meds[meds != "eid"]
    joindf <- dplyr::inner_join(joindf, medsdf, by = "eid")
    dxdf <- tidyr::pivot_longer(diseasedf, -eid, names_to = "disease", values_to = "status", values_drop_na = TRUE)
    joindf <- dplyr::inner_join(joindf, dxdf, by = "eid")
    joindf <- tidyr::drop_na(joindf)
    joindf$cluster <- relevel(factor(joindf$cluster), ref = "cluster_0")
    purrr::map_dfr(
        thresh,
        function(TR){
            joindfthresh <- joindf[joindf$cluster_prob >= TR,]
            joindfthresh <- dplyr::group_by(joindfthresh, disease)
            joindfthresh <- tidyr::nest(joindfthresh)
            joindfthresh <- dplyr::transmute(
                joindfthresh,
                data = purrr::map(
                    data,
                    function(d){
                        countdat <- data.frame(thresh = TR, table(d$cluster))
                        names(countdat) <- c("thresh", "cluster", "N")
                        mod1 <- glm(reformulate(response = "status", termlabels = "cluster"), data = d, family = "binomial")
                        mod2 <- glm(reformulate(response = "status", termlabels = c("cluster", BODYSIZEINDEX, COVARIATES, meds)), data = d, family = "binomial")
                        modres <- dplyr::bind_rows(list(Unadjusted = extractclusglmcoefs(mod1), Adjusted = extractclusglmcoefs(mod2)), .id = "model")
                        modres$N0 <- countdat$N[countdat$cluster == "cluster_0"]
                        modres$cluster <- gsub("^cluster", "", modres$cluster)
                        dplyr::inner_join(countdat, modres, by = "cluster")
                    }
                )
            )
            joindfthresh <- tidyr::unnest(joindfthresh, data)
            dplyr::ungroup(joindfthresh)
        }
    )
}

# Models of association between cluster labels given at a probability threshold and a disease by sex
get_adjdiseasesummary_clustersthresh <- function(clusres, diseasedat, medsdat, sexstrats, thresh = c(.6, .7, .8)){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ adjclusdiseasesumfx_thresh(clusres[[sexgr]], diseasedat, medsdat, sexstrats[[sexgr]], thresh) },
        .id = "sex"
    )
}

# Incidence
incidencesumfx <- function(incidencedf, X){
    incdf <- incidencedf[incidencedf$eid %in% X$eid,]
    ntotal <- length(unique(incdf$eid))
    ncases <- sum(incdf$outcome_value, na.rm = TRUE)
    personyrs <- sum(incdf$outcome_timeyrs, na.rm = TRUE)
    avgtime <- mean(incdf$outcome_timeyrs, na.rm = TRUE)
    sdtime <- sd(incdf$outcome_timeyrs, na.rm = TRUE)
    mediqrtime <- quantile(incdf$outcome_timeyrs, probs = c(.25, .5, .75), na.rm = TRUE, names = FALSE)
    data.frame(
        ntotal = ntotal,
        ncases = ncases, 
        personyrs = personyrs,
        casesper1e5py = 1e5 * ncases / personyrs,
        avgtime = avgtime,
        sdtime = sdtime,
        medtime = mediqrtime[2],
        timeq25 = mediqrtime[1],
        timeq75 = mediqrtime[3]
    )
}

# Incidence by sex
get_incidencesummary <- function(incidencedat, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ incidencesumfx(incidencedat, sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Incidence by cluster
clusincidencesumfx <- function(incidencedf, clusdat){
    probs <- clusdat$probs
    joindf <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    joindf <- dplyr::inner_join(joindf, incidencedf, by = "eid")
    joindf <- tidyr::drop_na(joindf)
    joindf <- dplyr::group_by(joindf, cluster)
    dplyr::summarise(
        joindf,
        ncases = sum(outcome_value * cluster_prob),
        personyrs = sum(outcome_timeyrs * cluster_prob),
        casesper1e5py = 1e5 * ncases / personyrs,
        .groups = "drop"
    )
}

# Incidence by cluster in each sex
get_incidencesummary_clusters <- function(incidencedat, clusres){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ clusincidencesumfx(incidencedat, clusres[[sexgr]]) },
        .id = "sex"
    )
}

# Models of association between cluster allocation probabilities and incidence
adjclusincidencesumfx <- function(clusdat, incidencedf, medsdf, X){
    probs <- clusdat$probs
    clusnams <- names(probs)
    clusnams <- clusnams[!(clusnams %in% c("eid", "cluster_0"))]
    joindf <- probs[,c("eid", clusnams)]
    for(clus in clusnams){
        joindf[[clus]] <- ptologp(joindf[[clus]])
    }
    basecovardf <- X[,c("eid", BODYSIZEINDEX, COVARIATES)]
    joindf <- dplyr::inner_join(joindf, basecovardf, by = "eid")
    meds <- names(medsdf)
    meds <- meds[meds != "eid"]
    joindf <- dplyr::inner_join(joindf, medsdf, by = "eid")
    joindf <- dplyr::inner_join(joindf, incidencedf, by = "eid")
    joindf <- tidyr::drop_na(joindf)
    survterm <- "survival::Surv(time = outcome_timeyrs, event = outcome_value)"
    mod1 <- survival::coxph(reformulate(response = survterm, termlabels = clusnams), data = joindf)
    mod2 <- survival::coxph(reformulate(response = survterm, termlabels = c(clusnams, BODYSIZEINDEX, COVARIATES, meds)), data = joindf)
    modres <- dplyr::bind_rows(list(Unadjusted = extractclusglmcoefs(mod1), Adjusted = extractclusglmcoefs(mod2)), .id = "model")
    tibble::tibble(modres)
}

# Models of association between cluster allocation probabilities and incidence by sex
get_adjincidencesummary_clusters <- function(clusres, incidencedat, medsdat, sexstrats){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ adjclusincidencesumfx(clusres[[sexgr]], incidencedat, medsdat, sexstrats[[sexgr]]) },
        .id = "sex"
    )
}

# Models of association between cluster labels given at a probability threshold and incidence
adjclusincidencesumfx_thresh <- function(clusdat, incidencedf, medsdf, X, thresh){
    probs <- clusdat$probs
    joindf <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
    basecovardf <- X[,c("eid", BODYSIZEINDEX, COVARIATES)]
    joindf <- dplyr::inner_join(joindf, basecovardf, by = "eid")
    meds <- names(medsdf)
    meds <- meds[meds != "eid"]
    joindf <- dplyr::inner_join(joindf, medsdf, by = "eid")
    joindf <- dplyr::inner_join(joindf, incidencedf, by = "eid")
    joindf <- tidyr::drop_na(joindf)
    joindf$cluster <- relevel(factor(joindf$cluster), ref = "cluster_0")
    survterm <- "survival::Surv(time = outcome_timeyrs, event = outcome_value)"
    purrr::map_dfr(
        thresh,
        function(TR){
            joindfthresh <- joindf[joindf$cluster_prob >= TR,]
            countdat <- data.frame(thresh = TR, table(joindfthresh$cluster))
            names(countdat) <- c("thresh", "cluster", "N")
            mod1 <- survival::coxph(reformulate(response = survterm, termlabels = "cluster"), data = joindfthresh)
            mod2 <- survival::coxph(reformulate(response = survterm, termlabels = c("cluster", BODYSIZEINDEX, COVARIATES, meds)), data = joindfthresh)
            modres <- dplyr::bind_rows(list(Unadjusted = extractclusglmcoefs(mod1), Adjusted = extractclusglmcoefs(mod2)), .id = "model")
            modres$N0 <- countdat$N[countdat$cluster == "cluster_0"]
            modres$cluster <- gsub("^cluster", "", modres$cluster)
            dplyr::inner_join(countdat, modres, by = "cluster")
        }
    )
}

# Models of association between cluster labels given at a probability threshold and incidence by sex
get_adjincidencesummary_clustersthresh <- function(clusres, incidencedat, medsdat, sexstrats, thresh = c(.6, .7, .8)){
    purrr::map_dfr(
        setNames(SEXGROUPS, SEXGROUPS),
        function(sexgr){ adjclusincidencesumfx_thresh(clusres[[sexgr]], incidencedat, medsdat, sexstrats[[sexgr]], thresh) },
        .id = "sex"
    )
}

## Global survival model
global_survivalmodel <- function(sexstrat, incidencedf){
    purrr::map(
        sexstrat,
        function(X){
            DAT <- dplyr::inner_join(X, incidencedf, by = "eid")
            survterm <- "survival::Surv(time = outcome_timeyrs, event = outcome_value)"
            survival::coxph(reformulate(response = survterm, termlabels = PERSEXVARS), data = DAT)
        }
    )
}

# Weighted ROC curve
wtd_roc <- function(label, preds, weights = rep(1, length(label))){
    ord <- order(preds, decreasing = TRUE)
    y <- label[ord]
    y.hat <- preds[ord]
    cases <- y == 1
    controls <- y == 0
    w.cases <- w.controls <- weights[ord]
    w.cases[controls] <- 0
    w.controls[cases] <- 0
    cum.cases <- cumsum(w.cases)
    cum.controls <- cumsum(w.controls)
    is.end <- c(diff(y.hat) != 0, TRUE)
    total.cases <- cum.cases[length(y)]
    total.controls <- cum.controls[length(y)]
    thresholds <- c(y.hat[is.end], 0)
    TN <- c(0, cum.controls[is.end])
    FN <- c(0, cum.cases[is.end])
    FP <- total.controls - TN
    TP <- total.cases - FN
    data.frame(thresholds, TN, TP, FP, FN)
}

# Weighted ROC AUC variance parameters
wtd_rocaucvar <- function(labels, preds, weights = rep(1, length(labels))) {
    cases <- preds[labels == 1]
    controls <- preds[labels == 0]
    w.cases <- weights[labels == 1]
    w.controls <- weights[labels == 0]
    n <- sum(w.controls)
    m <- sum(w.cases)
    wmat <- outer(w.cases, w.controls, "*")
    equalsurv  <- outer(cases, controls, "==") * 0.5
    lowersurv <- outer(cases, controls, "<") * 1.0
    MW <- (equalsurv + lowersurv) * wmat
    V <- list(
        theta = sum(MW)/(m * n),
        X = rowSums(MW) / n,
        n = n,
        Y = colSums(MW) / m,
        m = m
    )
    return(V)
}

# Global survival metrics at 5 years
global_survmetrics <- function(survmod, sexstrat, incidencedf){
    purrr::map2(
        survmod, sexstrat,
        function(mod, X){
            STATUS_Y5 <- dplyr::transmute(
                incidencedf, 
                eid,
                outcome_value = ifelse(outcome_value == 1 & outcome_timeyrs > 5, 0, outcome_value),
                outcome_timeyrs = 5
            )
            STATUS_Y5 <- dplyr::inner_join(STATUS_Y5, X, by = "eid")
            STATUS_Y5$pred_survprob <- survival:::predict.coxph(object = mod, newdata = STATUS_Y5, type = "survival")
            roctab <- wtd_roc(label = STATUS_Y5$outcome_value, preds = STATUS_Y5$pred_survprob)
            rocaucvar <- wtd_rocaucvar(label = STATUS_Y5$outcome_value, preds = STATUS_Y5$pred_survprob)
            return(list(roctab = roctab, rocaucvar = rocaucvar))
        }
    )
}

# Global survival metrics at 5 years
cluster_survmetrics <- function(survmod, sexstrat, incidencedf, clusres){
    purrr::pmap(
        list(survmod, sexstrat, clusres),
        function(mod, X, clusdat){
            STATUS_Y5 <- dplyr::transmute(
                incidencedf, 
                eid,
                outcome_value = ifelse(outcome_value == 1 & outcome_timeyrs > 5, 0, outcome_value),
                outcome_timeyrs = 5
            )
            STATUS_Y5 <- dplyr::inner_join(STATUS_Y5, X, by = "eid")
            STATUS_Y5$pred_survprob <- survival:::predict.coxph(object = mod, newdata = STATUS_Y5, type = "survival")
            probs <- clusdat$probs
            probs <- tidyr::pivot_longer(probs, -eid, names_to = "cluster", values_to = "cluster_prob")
            STATUS_Y5 <- dplyr::inner_join(STATUS_Y5, probs, by = "eid")
            STATUS_Y5 <- split(STATUS_Y5, ~cluster)
            purrr::map(
                STATUS_Y5,
                function(cl){
                    roctab <- wtd_roc(label = cl$outcome_value, preds = cl$pred_survprob, weights = cl$cluster_prob)
                    rocaucvar <- wtd_rocaucvar(label = cl$outcome_value, preds = cl$pred_survprob, weight = cl$cluster_prob)
                    return(list(roctab = roctab, rocaucvar = rocaucvar))
                }
            )
        }
    )
}

# ROC AUC 
aucfx <- function(roctab){
    roctab$FPR <- roctab$FP / (roctab$FP + roctab$TN)
    roctab$TPR <- roctab$TP / (roctab$TP + roctab$FN)
    right <- roctab[-nrow(roctab),]
    left <- roctab[-1,]
    width <- right$FPR - left$FPR
    rect.area <- left$TPR * width
    triangle.h <- right$TPR - left$TPR
    triangle.area <- triangle.h * width / 2
    sum(rect.area, triangle.area)
}

# Confidence intervals of ROC AUC curves
wtd_ci.auc <- function(Vobj, conf.level){
    with(
        Vobj,
        {
            SX <- sum((X - theta)^2)/(m-1)
            SY <- sum((Y - theta)^2)/(n-1)
            S <- SX/m + SY/n
            qnorm(c((1-conf.level)/2, 1-((1-conf.level)/2)), mean = theta, sd = sqrt(S))
        }
    )
}

# Comparing two *independent* ROC AUC curves
compare_rocs <- function (VR, VS){
    nR <- VR$n
    mR <- VR$m
    nS <- VS$n
    mS <- VS$m
    SRX <- sum((VR$X - VR$theta)^2)/(mR-1)
    SSX <- sum((VS$X - VS$theta)^2)/(mS-1)
    SRY <- sum((VR$Y - VR$theta)^2)/(nR-1)
    SSY <- sum((VS$Y - VS$theta)^2)/(nS-1)
    SR <- SRX/mR + SRY/nR
    SS <- SSX/mS + SSY/nS
    ntotR <- nR + mR
    ntotS <- nS + mS
    SSR <- sqrt((SR) + (SS))
    t <- (VR$theta - VS$theta) / SSR
    df <- ((SR + SS)^2) / ((SR^2 / (ntotR - 1)) + (SS^2 / (ntotS - 1)))
    return(c(t, df))
}
