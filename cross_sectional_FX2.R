## Functions to assess predictive ability of clusters

## Columns needed
BIOMARKERS <- c("whr", "sbp", "dbp", "alt", "scr", "crp", "hdl", "tg", "ldl", "fg")
BODYSIZEINDEX <- "bmi"
COVARIATES <- c("age", "smoking")
DISEASES <- c("HT", "CHD", "Stroke", "PAD", "CKD", "LiverFailure", "RA", "T2D", "T1D")
MEDICATION <- c("Insulin", "AntiDM", "AntiHT", "LipidLower")

# Calculating cluster probabilities
clusterprobcalc <- function(X){
    RES <- dplyr::mutate(
        X,
        # Calculate standardized residuals
        resids = purrr::map2(
            data, residmod,
            function(DATA, RESIDMOD){
                sapply(
                    BIOMARKERS,
                    function(biomarker){
                        PREDICTDAT <- as.matrix(cbind(Intercept = 1, DATA[,c(BODYSIZEINDEX, COVARIATES)]))
                        MODELDAT <- dplyr::filter(RESIDMOD, Biomarker == biomarker)
                        RESIDSD <- MODELDAT$SDRes
                        COEFDAT <- MODELDAT[,c("Intercept", BODYSIZEINDEX, COVARIATES)]
                        COEFDAT <- t(COEFDAT)
                        PREDVAL <- as.vector(PREDICTDAT %*% COEFDAT)
                        RESIDVAL <- DATA[[biomarker]] - PREDVAL
                        RESIDVAL / RESIDSD
                    }
                )
            }
        ),
        # Calculate probabilities to valid clusters
        probs = purrr::map2(
            resids, clusmod,
            function(RESIDDF, CLUSMOD){
                mus <- CLUSMOD$validclus_centers
                mus <- purrr::map(mus, function(mu){ mu[BIOMARKERS] })
                covmats <- CLUSMOD$validclus_covmats
                covmats <- purrr::map(covmats, function(covmat){ covmat[BIOMARKERS, BIOMARKERS] })
                wts <- CLUSMOD$validclus_weights
                pdfs <- mapply(function(mu, covmat){ mvtnorm::dmvnorm(RESIDDF, mu, covmat) }, mus, covmats)
                L <- pdfs %*% diag(wts)
                colnames(L) <- paste0("prob", CLUSMOD$validclus_name)
                data.frame(L / rowSums(L))
            }
        ),
        data = purrr::map2(data, probs, bind_cols)
    )
    dplyr::select(RES, -c(resids, probs))
}

# CLR transformation
clrcalc <- function(X){
    dplyr::mutate(
        X,
        data = purrr::map(
            data,
            ~{
                lower_bound <- 1e-10
                upper_bound <- 1 - lower_bound
                RES <- dplyr::mutate(
                    .x,
                    dplyr::across(
                        dplyr::starts_with("prob"),
                        function(CLUS){
                            CLUS <- ifelse(CLUS < lower_bound, lower_bound, CLUS)
                            CLUS <- ifelse(CLUS > upper_bound, upper_bound, CLUS)
                            return(CLUS)
                        },
                        .names = "clr{.col}"
                    ),
                    Gmean = exp(rowMeans(log(dplyr::pick(dplyr::starts_with("clr"))))),
                    dplyr::across(dplyr::starts_with("clr"), ~log(.x/Gmean))
                )
                RES <- dplyr::rename_with(
                    RES,
                    function(COLN){ gsub("prob", "", COLN) },
                    dplyr::starts_with("clrprob")
                )
                dplyr::select(RES, -clrBC)
            }
        )
    )
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
            sum((x == LEVEL) * p)
        }
    )
}

# Weighted summary of categorical variables
wtd_categorical_sumstats <- function(x, p){
    vals <- wtd_table(x, p)
    vals <- sort(vals, decreasing = TRUE)
    tibble::tibble(
        N = sum(p),
        Summary1 = names(vals), 
        Summary2 = as.character(round(vals, 2))
    )
}

# Weighted summary of continuous variables
wtd_continuous_sumstats <- function(x, p){
    Mn <- round(wtd_mean(x, p), 2)
    Sd <- round(sqrt(wtd_var(x, p)), 2)
    quants <- wtd_quantile(x, p, q = c(.025, .25, .5, .75, .975))
    Md <- round(quants[3], 2)
    pc <- round(quants[-3], 2)
    pc <- paste(pc, collapse = " - ")
    tibble::tibble(
        N = sum(p),
        Summary1 = paste0(Mn, " (", Sd, ")"),
        Summary2 = paste0(Md, " (", pc, ")")
    )
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
markerdistribfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::map(
            data,
            ~{
                RES2 <- tidyr::expand_grid(
                    Variable = c(BIOMARKERS, COVARIATES),
                    Cluster = names(.x)[grepl("^prob", names(.x))]
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    RES3 = purrr::map2(
                        Variable, Cluster,
                        function(COL, PROB, DAT){
                            wtd_sumstatfx(DAT[[COL]], DAT[[PROB]])
                        },
                        DAT = .x
                    )
                )
                tidyr::unnest(RES2, RES3)
            }
        )
    )
    RES <- tidyr::unnest(RES, RES1)
    dplyr::mutate(RES, Cluster = gsub("prob", "", Cluster))
}

# BMI effect on biomarkers specific to each cluster
bmieffmarkerfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::map(
            data,
            ~{
                RES2 <- tidyr::expand_grid(
                    Variable = BIOMARKERS,
                    Cluster = names(.x)[grepl("^prob", names(.x))]
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    RES3 = purrr::map2(
                        Variable, Cluster,
                        function(COL, PROB, DAT, COVARS){
                            D <- dplyr::select(DAT, all_of(c(COL, COVARS)))
                            mod <- lm(reformulate(response = COL, termlabels = "."),
                                      data = D, weights = DAT[[PROB]])
                            modcoefs <- coef(summary(mod))
                            modcoefs <- data.frame(term = rownames(modcoefs), modcoefs)
                            modcoefs <- dplyr::select(
                                modcoefs,
                                term, estimate = Estimate, se = Std..Error
                            )
                            return(modcoefs)
                        },
                        DAT = .x, COVARS = c(COVARIATES, BODYSIZEINDEX)
                    )
                )
                tidyr::unnest(RES2, RES3)
            }
        )
    )
    RES <- tidyr::unnest(RES, RES1)
    dplyr::mutate(RES, Cluster = gsub("prob", "", Cluster))
}

# Counts of individuals with diseases/medication per cluster
countcovarsfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        countres = purrr::map(
            data,
            function(DAT){
                clusnams <- names(DAT)[grepl("^prob", names(DAT))]
                DXandMed <- c(DISEASES, MEDICATION)
                RES1 <- tibble::tibble(
                    Cluster = "Overall",
                    Covariate = DXandMed
                )
                RES1 <- dplyr::mutate(
                    RES1,
                    Nclus = nrow(DAT),
                    NclusDX = purrr::map_dbl(Covariate, ~sum(DAT[[.x]]))
                )
                RES2 <- tidyr::expand_grid(
                    Cluster = clusnams,
                    Covariate = DXandMed
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    Nclus = purrr::map_dbl(Cluster, ~sum(DAT[[.x]])),
                    NclusDX = purrr::map2_dbl(
                        Cluster, Covariate,
                        function(CL, CO, D){
                            sum(D[[CO]] * D[[CL]])
                        },
                        D = DAT
                    )
                )
                RES3 <- dplyr::bind_rows(RES1, RES2)
                dplyr::mutate(RES3, Cluster = gsub("prob", "", Cluster))
            }
        )
    )
    tidyr::unnest(RES, countres)
}

# Counts of individuals with specific diseases and medication combinations per cluster
countspectxfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        countres = purrr::map(
            data,
            function(DAT){
                clusnams <- names(DAT)[grepl("^prob", names(DAT))]
                combs <- tibble::tibble(
                    DX = c("HT", "CHD", "T2D", "T2D"),
                    MED = list("AntiHT", "LipidLower", "Insulin",
                               c("Insulin", "AntiDM"))
                )
                RES1 <- dplyr::mutate(
                    combs,
                    Cluster = "Overall",
                    NclusDXM = purrr::map2_dbl(
                        DX, MED, 
                        ~{
                            disease <- DAT[[.x]] == 1
                            anymed <- rowSums(DAT[,.y]) > 0
                            sum(disease & anymed)
                        }
                    )
                )
                RES2 <- tidyr::expand_grid(
                    combs,
                    Cluster = clusnams
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    NclusDXM = purrr::pmap_dbl(
                        list(Cluster, DX, MED),
                        ~{
                            Prob <- DAT[[..1]]
                            disease <- DAT[[..2]] == 1
                            anymed <- rowSums(DAT[,..3]) > 0
                            sum((disease & anymed) * Prob)
                        }
                    )
                )
                RES3 <- dplyr::bind_rows(RES1, RES2)
                dplyr::mutate(
                    RES3,
                    Cluster = gsub("prob", "", Cluster),
                    MED = purrr::map_chr(MED, paste, collapse = "Or")
                )
            }
        )
    )
    tidyr::unnest(RES, countres)
}

# Logistic regressions with cluster CLRs against current diseases
assocdxfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::map(
            data,
            function(DAT){
                clusnames <- names(DAT)[grepl("^clr", names(DAT))]
                RES2 <- tibble::tibble(DX = DISEASES)
                RES2 <- dplyr::mutate(
                    RES2,
                    modres = purrr::map(
                        DX,
                        function(dx, D, clusnm){
                            if(sum(D[[dx]]) > 5){
                                mod1 <- glm(
                                    formula = reformulate(
                                        termlabels = clusnm, 
                                        response = dx
                                    ), 
                                    family = binomial, data = D
                                )
                                mod2 <- glm(
                                    formula = reformulate(
                                        termlabels = c(clusnm, MEDICATION), 
                                        response = dx
                                    ), 
                                    family = binomial, data = D
                                )
                                bind_rows(
                                    OnlyClusters = dplyr::select(broom::tidy(mod1), term, estimate, se = std.error),
                                    ClustersMed = dplyr::select(broom::tidy(mod2), term, estimate, se = std.error),
                                    .id = "model"
                                )
                            } else {
                                data.frame(term = NA, estimate = NaN)
                            }
                        },
                        D = DAT, clusnm = clusnames
                    )
                )
                tidyr::unnest(RES2, modres)
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

# Adding survival data - MACE
addsurvmacedat <- function(DAT, SURVDATA){
    D <- dplyr::inner_join(DAT, SURVDATA, by = "eid")
    dplyr::filter(
        D,
        CHD == 0,
        Stroke == 0,
        PAD == 0
    )
}

# Adding survival data - Diabetes
addsurvdmdat <- function(DAT, SURVDATA){
    D <- dplyr::inner_join(DAT, SURVDATA, by = "eid")
    dplyr::filter(
        D,
        T2D == 0,
        T1D == 0,
        Insulin == 0,
        AntiDM == 0
    )
}

# Kaplan-Meier estimates
kmestfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        dplyr::across(
            c(macedf, dmdf),
            function(OUTDFs){
                purrr::map(
                    OUTDFs,
                    function(OUTDF){
                        RES1 <- survival:::survfit.formula(
                            survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                            data = OUTDF
                        )
                        RES1 <- survival:::summary.survfit(RES1, times = 10)
                        RES1 <- RES1[c("surv", "upper", "lower")]
                        RES1 <- 1 - data.frame(RES1)
                        colnames(RES1) <- c("risk", "lower", "upper")
                        RES1 <- cbind(Cluster = "Overall", RES1)
                        RES2 <- tibble::tibble(Cluster = names(OUTDF)[grepl("^prob", names(OUTDF))])
                        RES2 <- dplyr::mutate(
                            RES2,
                            RES3 = purrr::map(
                                Cluster,
                                function(CL, D){
                                    RES4 <- survival:::survfit.formula(
                                        survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                                        data = D,
                                        weights = get(CL),
                                        robust = FALSE
                                    )
                                    RES4 <- survival:::summary.survfit(RES4, times = 10)
                                    RES4 <- RES4[c("surv", "upper", "lower")]
                                    RES4 <- 1 - data.frame(RES4)
                                    colnames(RES4) <- c("risk", "lower", "upper")
                                    RES4
                                },
                                D = OUTDF
                            )
                        )
                        RES2 <- tidyr::unnest(RES2, RES3)
                        dplyr::bind_rows(RES1, RES2)
                    }
                )
            }
        ),
        RES4 = purrr::map2(
            macedf, dmdf,
            ~dplyr::bind_rows(MACE = .x, DM = .y, .id = "Outcome")
        )
    )
    RES <- dplyr::select(RES, sex, RES4)
    tidyr::unnest(RES, RES4)
}

# Running Cox models - MACE
coxmodelsmace <- function(X){
    dplyr::transmute(
        X,
        sex, macedf,
        score2 = purrr::map2(
            macedf, sex,
            function(DATA, SEX){
                MODELDAT <- dplyr::transmute(
                    DATA,
                    outcome_value,
                    outcome_timeyrs,
                    age,
                    smoking,
                    sbp,
                    T2D,
                    tchol = hdl + ldl + (tg/2.2),
                    hdl,
                    age_smoking = age*smoking,
                    age_sbp = age*sbp,
                    age_t2d = age*T2D,
                    age_tchol = age*tchol,
                    age_hdl = age*hdl,
                    t2donsetage = ifelse(T2D == 1, T2Dage, 0),
                    fg,
                    scr = scr / 88.42,
                    alpha = dplyr::case_when(SEX == "Female" ~ -0.241, SEX == "Male" ~ -0.302),
                    kappa = dplyr::case_when(SEX == "Female" ~ 0.7, SEX == "Male" ~ 0.9),
                    creat_kappa = scr / kappa,
                    minkappa = pmin(creat_kappa, 1), 
                    maxkappa = pmax(creat_kappa, 1),
                    minkappa_alpha = minkappa^alpha,
                    maxkappa_exp = maxkappa^(-1.200),
                    age_term = 0.9938^age,
                    sex_term = dplyr::case_when(SEX == "Female" ~ 1.012, SEX == "Male" ~ 1), 
                    egfr = 142 * minkappa_alpha * maxkappa_exp * age_term * sex_term,
                    lnegfr = log(egfr),
                    lnegfrsq = lnegfr^2,
                    age_fg = age*fg,
                    age_lnegfr = age*lnegfr,
                    bmi, 
                    whr, alt, crp,
                    HT, CKD, LiverFailure, RA, T1D,
                    Insulin, AntiDM, AntiHT, LipidLower,
                    dplyr::across(dplyr::starts_with("clr"))
                )
                dplyr::select(
                    MODELDAT, 
                    -c(scr, alpha, kappa, creat_kappa, minkappa, maxkappa, 
                       minkappa_alpha, maxkappa_exp, age_term, sex_term, egfr)
                )
            }
        ),
        mod_null = purrr::map(
            score2, 
            ~survival::coxph(
                survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                data = .x
            )
        ),
        mod_score2 = purrr::map(
            score2,
            ~{
                MODELDAT <- dplyr::select(.x, -starts_with("clr"))
                survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = MODELDAT)
            }
        ),
        mod_score2clr = purrr::map(
            score2,
            ~survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = .x)
        )
    )
}

macesurvcoefx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        dplyr::across(
            c(mod_score2, mod_score2clr),
            ~purrr::map(
                .x, 
                function(MOD){ 
                    RES <- coef(summary(MOD))[,c(1,3)]
                    colnames(RES) <- c("estimate", "se")
                    data.frame(term = rownames(RES), RES)
                }
            )
        ),
        RES1 = purrr::map2(
            mod_score2, mod_score2clr,
            ~bind_rows(score2 = .x, score2clr = .y, .id = "model")
        ),
        mod_score2 = NULL,
        mod_score2clr = NULL
    )
    tidyr::unnest(RES, RES1)
}

comparemodsmace <- function(X){
    dplyr::transmute(
        X,
        sex,
        LL0 = purrr::map_dbl(mod_null, logLik),
        dplyr::across(
            c(mod_score2, mod_score2clr),
            ~purrr::map_dbl(.x, logLik),
            .names = "LL{.col}"
        ),
        dplyr::across(
            c(mod_score2, mod_score2clr),
            ~purrr::map_dbl(.x, \(MOD) sum(!is.na(MOD$coefficients))),
            .names = "NV{.col}"
        ),
        LRTstat = 2 * abs(LLmod_score2clr - LLmod_score2),
        LRTdf = NVmod_score2clr - NVmod_score2,
        LRTp = stats::pchisq(
            q = LRTstat,
            df = LRTdf,
            lower.tail = FALSE
        ),
        AdeqInd = (LLmod_score2 - LL0) / (LLmod_score2clr - LL0),
        Cdat = purrr::map2(
            mod_score2, mod_score2clr,
            survival::concordance
        ),
        Cmod_score2 = purrr::map_dbl(Cdat, ~.x$concordance[1]),
        Cmodse_score2 = purrr::map_dbl(Cdat, ~diag(.x$var)[1]),
        Cmod_score2clr = purrr::map_dbl(Cdat, ~.x$concordance[2]),
        Cmodse_score2clr = purrr::map_dbl(Cdat, ~diag(.x$var)[1]),
        Cdiff = purrr::map_dbl(Cdat, ~c(1,-1) %*% coef(.x)),
        Cdiffse = purrr::map_dbl(Cdat, ~sqrt(c(1,-1) %*% vcov(.x) %*% c(1,-1))),
        Cdiffp = 2 * pnorm(-abs(Cdiff/Cdiffse)),
        Cdat = NULL
    )
}

AdeqIndClusMACEFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        dplyr::across(
            c(mod_score2, mod_score2clr),
            ~purrr::pmap(
                list(macedf, score2, .x),
                function(DAT, MODDAT, MOD){
                    clusnams <- names(DAT)[grepl("^prob", names(DAT))]
                    sapply(
                        DAT[clusnams],
                        function(PROBCL){
                            MODDAT$PROB <- PROBCL
                            MODNULLCL <- survival::coxph(
                                survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                                data = MODDAT,
                                weights = PROB,
                                robust = FALSE
                            )
                            LL0CL <- logLik(MODNULLCL)
                            MODCL <- survival::coxph(
                                formula = as.formula(MOD), 
                                data = MODDAT,
                                weights = PROB,
                                robust = FALSE,
                                init = coef(MOD),
                                control = survival::coxph.control(iter.max = 0)
                            )
                            LLclus <- as.numeric(logLik(MODCL))
                            LL0clus <- as.numeric(logLik(MODNULLCL))
                            2 * abs(LLclus - LL0clus)
                        }
                    )
                }
            )
        ),
        RES1 = purrr::map2(
            mod_score2, mod_score2clr, 
            ~pmin(.x/.y, 1)
        ),
        RES1 = purrr::map(
            RES1, ~tibble::tibble(Cluster = names(.x), AdeqInd = .x)
        ),
        mod_score2 = NULL, mod_score2clr = NULL
    )
    tidyr::unnest(RES, RES1)
}

AdeqIndByPreMACEFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(score2, mod_score2, mod_score2clr),
            function(DAT, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DAT, outcome_timeyrs = 10)
                yhat1 <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                sapply(
                    setNames(seq(0, .25, .05), seq(0, .25, .05)),
                    function(thresh){
                        TDAT <- DAT[yhat1 >= thresh,]
                        MODNULL <- survival::coxph(
                            formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1,
                            data = TDAT
                        )
                        NEWMOD1 <- survival::coxph(
                            formula = as.formula(MOD1),
                            data = TDAT,
                            init = coef(MOD1), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        NEWMOD2 <- survival::coxph(
                            formula = as.formula(MOD2),
                            data = TDAT,
                            init = coef(MOD2), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        LL0 <- as.numeric(logLik(MODNULL))
                        LL1 <- as.numeric(logLik(NEWMOD1))
                        LL2 <- as.numeric(logLik(NEWMOD2))
                        LR1 <- 2*abs(LL1 - LL0)
                        LR2 <- 2*abs(LL2 - LL0)
                        pmin(1, LR1/LR2)
                    }
                )
            }
        ),
        RES1 = purrr::map(
            RES1, ~tibble::tibble(threshold = names(.x), AdeqInd = .x)
        )
    )
    tidyr::unnest(RES, RES1)
}

surv_to_risk <- function(mat, time, robust){
    res <- survival:::survfit.formula(
        mat[,1] ~ 1, 
        data = mat,
        weights = mat[,"weights"],
        robust = robust
    )
    res <- broom::tidy(res)
    res <- dplyr::mutate(res, estimate = 1 - estimate)
    res <- dplyr::filter(res, time <= time)
    res <- dplyr::slice_tail(res)
    res$estimate
}

dca <- function(formula, data, thresholds, time, weights = rep(1, nrow(data)), robust = FALSE){
    model_frame <- stats::model.frame(formula, data)
    outcome_name <- names(model_frame)[1]
    preds <- setdiff(names(model_frame), outcome_name)
    model_frame$all <- 1
    model_frame$none <- 0
    model_frame <- dplyr::mutate(
        model_frame,
        dplyr::across(
            -dplyr::all_of(outcome_name),
            ~dplyr::case_when(.x == 0 ~ 0 - .Machine$double.eps,
                              .x == 1 ~ 1 + .Machine$double.eps,
                              TRUE ~ .x)
        )
    )
    model_frame$weights <- weights
    outcome <- model_frame[[outcome_name]]
    RES <- tidyr::expand_grid(
        pred = c(preds, "all", "none"),
        threshold = thresholds
    )
    RES <- dplyr::mutate(
        RES,
        n = sum(weights),
        pos_rate = surv_to_risk(model_frame, time, robust),
        test_pos_rate = purrr::map2_dbl(
            pred, threshold,
            function(P, T, MODF){
                Pr <- MODF[,P]
                W <- MODF[,"weights"]
                TB <- wtd_table(Pr >= T, p = W)
                NPOS <- TB["TRUE"]
                if(is.na(NPOS)){ 
                    return(0)
                } else {
                    return(NPOS / sum(W))
                }
            },
            MODF = model_frame
        ),
        risk_rate_among_test_pos = purrr::map2_dbl(
            pred, threshold,
            function(P, T, MODF, time, robust){
                Pr <- MODF[,P]
                MODFPOS <- MODF[Pr >= T,]
                tryCatch(
                    surv_to_risk(MODFPOS, time, robust),
                    error = function(e) {
                        if(sum(Pr >= T) == 0L){ return(0) }
                        NA_real_
                    }
                ) 
            },
            MODF = model_frame, time = time, robust = robust
        ),
        tp_rate = risk_rate_among_test_pos * test_pos_rate,
        fp_rate = (1 - risk_rate_among_test_pos) * test_pos_rate,
        net_benefit = tp_rate - threshold / (1 - threshold) * fp_rate
    )
    RES <- dplyr::select(RES, pred, n, threshold, pos_rate, tp_rate, fp_rate, net_benefit)
    allnb <- dplyr::filter(RES, pred == "all")
    allnb <- dplyr::select(allnb, threshold, net_benefit_all = net_benefit)
    RES <- dplyr::left_join(RES, allnb, by = "threshold")
    dplyr::mutate(
        RES,
        net_intervention_avoided = (net_benefit - net_benefit_all) /
        (threshold / (1 - threshold)),
        net_benefit_all = NULL
    )
}

DCurvMACEFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(score2, mod_score2, mod_score2clr),
            function(DAT, MOD1, MOD2){
                NEWDAT <- mutate(DAT, outcome_timeyrs = 10)
                DAT$score2_y <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DAT$score2clr_y <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                dca(
                    formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ score2_y + score2clr_y,
                    data = DAT,
                    thresholds = seq(0, .25, .01),
                    time = 10
                )
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

DCurvMACEbyClFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(macedf, score2, mod_score2, mod_score2clr),
            function(DAT, DATMOD, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DATMOD, outcome_timeyrs = 10)
                DATMOD$score2_y <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DATMOD$score2clr_y <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                RES2 <- tibble::tibble(
                    Cluster = names(DAT)[grepl("^prob", names(DAT))]
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    RES3 = purrr::map(
                        Cluster,
                        function(CL){
                            dca(
                                formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ score2_y + score2clr_y,
                                data = DATMOD,
                                thresholds = seq(0, .25, .01),
                                time = 10,
                                weights = DAT[[CL]]
                            )
                        }
                    )
                )
                tidyr::unnest(RES2, RES3)
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

################################################################################

# Running Cox models - DM
coxmodelsdm <- function(X){
    dplyr::transmute(
        X,
        sex, dmdf,
        baseclr = purrr::map(
            dmdf, 
            ~dplyr::select(
                .x,
                outcome_timeyrs, outcome_value,
                all_of(
                    c(
                        BIOMARKERS,
                        BODYSIZEINDEX,
                        COVARIATES,
                        DISEASES[!DISEASES %in% c("T1D", "T2D")],
                        MEDICATION[!MEDICATION %in% c("Insulin", "AntiDM")],
                        names(.x)[grepl("^clr", names(.x))]
                    )
                )
            )
        ),      
        mod_null = purrr::map(
            baseclr, 
            ~survival::coxph(
                survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                data = .x
            )
        ),
        mod_base = purrr::map(
            baseclr,
            ~{
                MODELDAT <- dplyr::select(.x, -starts_with("clr"))
                survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = MODELDAT)
            }
        ),
        mod_baseclr = purrr::map(
            baseclr,
            ~survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = .x)
        )
    )
}

dmsurvcoefx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        dplyr::across(
            c(mod_base, mod_baseclr),
            ~purrr::map(
                .x, 
                function(MOD){ 
                    RES <- coef(summary(MOD))[,c(1,3)]
                    colnames(RES) <- c("estimate", "se")
                    data.frame(term = rownames(RES), RES)
                }
            )
        ),
        RES1 = purrr::map2(
            mod_base, mod_baseclr,
            ~bind_rows(base = .x, baseclr = .y, .id = "model")
        ),
        mod_base = NULL, 
        mod_baseclr = NULL
    )
    tidyr::unnest(RES, RES1)
}

comparemodsdm <- function(X){
    dplyr::transmute(
        X,
        sex,
        LL0 = purrr::map_dbl(mod_null, logLik),
        dplyr::across(
            c(mod_base, mod_baseclr),
            ~purrr::map_dbl(.x, logLik),
            .names = "LL{.col}"
        ),
        dplyr::across(
            c(mod_base, mod_baseclr),
            ~purrr::map_dbl(.x, \(MOD) sum(!is.na(MOD$coefficients))),
            .names = "NV{.col}"
        ),
        LRTstat = 2 * abs(LLmod_baseclr - LLmod_base),
        LRTdf = NVmod_baseclr - NVmod_base,
        LRTp = stats::pchisq(
            q = LRTstat,
            df = LRTdf,
            lower.tail = FALSE
        ),
        AdeqInd = (LLmod_base - LL0) / (LLmod_baseclr - LL0),
        Cdat = purrr::map2(
            mod_base, mod_baseclr,
            survival::concordance
        ),
        Cmod_base = purrr::map_dbl(Cdat, ~.x$concordance[1]),
        Cmodse_base = purrr::map_dbl(Cdat, ~diag(.x$var)[1]),
        Cmod_baseclr = purrr::map_dbl(Cdat, ~.x$concordance[2]),
        Cmodse_baseclr = purrr::map_dbl(Cdat, ~diag(.x$var)[1]),
        Cdiff = purrr::map_dbl(Cdat, ~c(1,-1) %*% coef(.x)),
        Cdiffse = purrr::map_dbl(Cdat, ~sqrt(c(1,-1) %*% vcov(.x) %*% c(1,-1))),
        Cdiffp = 2 * pnorm(-abs(Cdiff/Cdiffse)),
        Cdat = NULL
    )
}

AdeqIndClusDMFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        dplyr::across(
            c(mod_base, mod_baseclr),
            ~purrr::pmap(
                list(dmdf, baseclr, .x),
                function(DAT, DATMOD, MOD){
                    clusnams <- names(DAT)[grepl("^prob", names(DAT))]
                    sapply(
                        DAT[clusnams],
                        function(PROBCL){
                            DATMOD$PROB <- PROBCL
                            MODNULLCL <- survival::coxph(
                                survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                                data = DATMOD,
                                weights = PROB,
                                robust = FALSE
                            )
                            LL0CL <- logLik(MODNULLCL)
                            MODCL <- survival::coxph(
                                formula = as.formula(MOD), 
                                data = DATMOD,
                                weights = PROB,
                                robust = FALSE,
                                init = coef(MOD),
                                control = survival::coxph.control(iter.max = 0)
                            )
                            LLclus <- as.numeric(logLik(MODCL))
                            LL0clus <- as.numeric(logLik(MODNULLCL))
                            2 * abs(LLclus - LL0clus)
                        }
                    )
                }
            )
        ),
        RES1 = purrr::map2(
            mod_base, mod_baseclr, 
            ~pmin(.x/.y, 1)
        ),
        RES1 = purrr::map(
            RES1, ~tibble::tibble(Cluster = names(.x), AdeqInd = .x)
        ),
        mod_base = NULL, mod_baseclr = NULL
    )
    tidyr::unnest(RES, RES1)
}

AdeqIndByPreDMFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(dmdf, mod_base, mod_baseclr),
            function(DAT, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DAT, outcome_timeyrs = 10)
                yhat1 <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                sapply(
                    setNames(seq(0, .25, .05), seq(0, .25, .05)),
                    function(thresh){
                        TDAT <- DAT[yhat1 >= thresh,]
                        MODNULL <- survival::coxph(
                            formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1,
                            data = TDAT
                        )
                        NEWMOD1 <- survival::coxph(
                            formula = as.formula(MOD1),
                            data = TDAT,
                            init = coef(MOD1), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        NEWMOD2 <- survival::coxph(
                            formula = as.formula(MOD2),
                            data = TDAT,
                            init = coef(MOD2), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        LL0 <- as.numeric(logLik(MODNULL))
                        LL1 <- as.numeric(logLik(NEWMOD1))
                        LL2 <- as.numeric(logLik(NEWMOD2))
                        LR1 <- 2*abs(LL1 - LL0)
                        LR2 <- 2*abs(LL2 - LL0)
                        pmin(1, LR1/LR2)
                    }
                )
            }
        ),
        RES1 = purrr::map(
            RES1, ~tibble::tibble(threshold = names(.x), AdeqInd = .x)
        )
    )
    tidyr::unnest(RES, RES1)
}

DCurvDMFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(dmdf, mod_base, mod_baseclr),
            function(DAT, MOD1, MOD2){
                NEWDAT <- mutate(DAT, outcome_timeyrs = 10)
                DAT$base_y <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DAT$baseclr_y <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                dca(
                    formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ base_y + baseclr_y,
                    data = DAT,
                    thresholds = seq(0, .25, .01),
                    time = 10
                )
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

DCurvDMbyClFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex,
        RES1 = purrr::pmap(
            list(dmdf, mod_base, mod_baseclr),
            function(DAT, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DAT, outcome_timeyrs = 10)
                DAT$base_y <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DAT$baseclr_y <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                RES2 <- tibble::tibble(
                    Cluster = names(DAT)[grepl("^prob", names(DAT))]
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    RES3 = purrr::map(
                        Cluster,
                        function(CL){
                            dca(
                                formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ base_y + baseclr_y,
                                data = DAT,
                                thresholds = seq(0, .25, .01),
                                time = 10,
                                weights = DAT[[CL]]
                            )
                        }
                    )
                )
                tidyr::unnest(RES2, RES3)
            }
        )
    )
    tidyr::unnest(RES, RES1)
}
