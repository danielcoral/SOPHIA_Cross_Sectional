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
    pc <- round(quantile(x, probs = c(.025, .25, .75, .975), na.rm = TRUE), 2)
    pc <- paste(pc, collapse = " - ")
    data.frame(N = N, N_miss = N_miss,
               Summary1 = paste0(Mn, " (", Sd, ")"),
               Summary2 = paste0(Md, " (", pc, ")"))
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
    modcoefs <- modsum$coefficients[,1:2, drop = FALSE]
    modCI <- confint(MOD)[,1:2, drop = FALSE]
    modcoefs <- cbind(modcoefs, modCI)
    modcoefs <- round(modcoefs, digits = 5)
    colnames(modcoefs) <- c("Estimate", "SE", "lowerCI", "upperCI")
    trm <- rownames(modcoefs)
    modcoefs <- data.frame(term = trm, modcoefs)
    rownames(modcoefs) <- NULL
    return(modcoefs)
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
        init = "pca", n_epochs = 500,
        binary_edge_weights = TRUE, dens_scale = 1, 
        ret_extra = c("model", "fgraph"), verbose = FALSE       
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

## Get UMAP models
get_umap_mod <- function(clusres){ purrr::map(clusres, "umap_model") }

# Only retain embedding in UMAP model
clean_cluster_results <- function(clusres){
    purrr::map(
        clusres, 
        function(CLUSRES){
            CLUSRES[["umap_model"]] <- list(embedding = CLUSRES[["umap_model"]][["embedding"]])
            CLUSRES
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
    quants <- wtd_quantile(x, p, q = c(.025, .25, .5, .75, .975))
    Md <- round(quants[3], 2)
    pc <- round(quants[-3], 2)
    pc <- paste(pc, collapse = " - ")
    data.frame(Summary1 = paste0(Mn, " (", Sd, ")"),
               Summary2 = paste0(Md, " (", pc, ")"))
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