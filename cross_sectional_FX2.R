## Functions to assess predictive ability of clusters

## Columns needed
BIOMARKERS <- c("whr", "sbp", "dbp", "alt", "scr", "crp", "hdl", "tg", "ldl", "fg")
BODYSIZEINDEX <- "bmi"
COVARIATES <- c("age", "smoking")
DISEASES <- c("HT", "CHD", "Stroke", "PAD", "CKD", "LiverFailure", "RA", "T2D", "T1D")
MEDICATION <- c("Insulin", "AntiDM", "AntiHT", "LipidLower")

# Calculating cluster probabilities
clusterprobcalc <- function(ClusModDf, StratDat){
    alldat <- dplyr::mutate(ClusModDf, data = purrr::map(sex, ~StratDat[[.x]]))
    dplyr::transmute(
        alldat,
        sex,
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
        data = purrr::map2(data, probs, bind_cols),
        resids = NULL, probs = NULL
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
                    Variable = c(BIOMARKERS, COVARIATES, BODYSIZEINDEX),
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

# Adding covariate data
addcovardat <- function(X, CovarDat){
    dplyr::mutate(
        X,
        data = purrr::map(data, inner_join, CovarDat, by = "eid")
    )
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
                RES1 <- tidyr::expand_grid(
                    Cluster = clusnams,
                    Covariate = DXandMed
                )
                dplyr::mutate(
                    RES1,
                    Nclus = purrr::map_dbl(Cluster, ~sum(DAT[[.x]], na.rm = TRUE)),
                    Ncases = purrr::map2_dbl(
                        Cluster, Covariate,
                        function(CL, CO, D){
                            sum(D[[CO]] * D[[CL]], na.rm = TRUE)
                        },
                        D = DAT
                    ),
                    Nnoncases = Nclus - Ncases
                )
            }
        )
    )
    tidyr::unnest(RES, countres)
}

# Counts of individuals with specific diseases and medication combinations per cluster
countdxmedfx <- function(X){
    RES <- dplyr::mutate(
        X,
        data = purrr::map(
            data,
            ~{
                RES1 <- tidyr::pivot_longer(.x, dplyr::starts_with("prob"), names_to = "Cluster", values_to = "prob")
                RES1 <- tidyr::pivot_longer(RES1, dplyr::all_of(DISEASES), names_to = "Dx", values_to = "Dx_value")
                RES1 <- dplyr::group_by(RES1, Dx, Dx_value, Cluster)
                RES1 <- dplyr::summarise(
                    RES1,
                    NoMed = sum((rowSums(dplyr::pick(dplyr::all_of(MEDICATION))) == 0) * prob),
                    dplyr::across(dplyr::all_of(MEDICATION), \(x) round(sum(x * prob))),
                    .groups = "drop"
                )
                RES1 <- tidyr::pivot_longer(RES1, c(NoMed, dplyr::all_of(MEDICATION)), names_to = "Med", values_to = "Med_value")
                RES1 <- dplyr::mutate(RES1, Dx_value = dplyr::case_match(Dx_value, 1 ~ "Ncases", 0 ~ "Nnoncases"))
                tidyr::pivot_wider(RES1, names_from = Dx_value, values_from = Med_value)
            }
        )
    )       
    tidyr::unnest(RES, data)
}

# Logistic regressions with cluster probabilities against current diseases
assocdxfx <- function(X){
    RES <- dplyr::mutate(
        X,
        data = purrr::map(
            data,
            function(DAT){
                DAT <- tidyr::pivot_longer(
                    DAT,
                    dplyr::all_of(DISEASES),
                    names_to = "Dx_name",
                    values_to = "Dx_value"
                )
                DAT <- tidyr::nest(DAT, DxDat = -Dx_name)
                DAT <- dplyr::filter(DAT, purrr::map_dbl(DxDat, \(x) sum(x$Dx_value)) > 5)
                DAT <- dplyr::mutate(
                    DAT,
                    DxDat = purrr::map(
                        DxDat,
                        function(DXD){
                            DXD <- dplyr::transmute(
                                DXD,
                                Dx_value,
                                dplyr::across(dplyr::starts_with("prob"), \(x) log(x/probBC), .names = "alr{.col}"),
                                alrprobBC = NULL,
                                dplyr::across(dplyr::all_of(c(BODYSIZEINDEX, COVARIATES, BIOMARKERS, MEDICATION)))
                            )
                            D1 <- dplyr::select(DXD, Dx_value, dplyr::starts_with("alr"))
                            MOD1 <- glm(Dx_value ~ ., family = binomial, data = D1)
                            RES1 <- tibble::tibble(
                                model = "OnlyClusters",
                                estimates = list(coef(MOD1)),
                                varcovmat = list(vcov(MOD1))
                            )
                            MOD2 <- glm(Dx_value ~ ., family = binomial, data = DXD)
                            RES2 <- tibble::tibble(
                                model = "FullModel",
                                estimates = list(coef(MOD2)),
                                varcovmat = list(vcov(MOD2))
                            )
                            bind_rows(RES1, RES2)
                        }
                    )
                )
                tidyr::unnest(DAT, DxDat)
            }
        )
    )
    tidyr::unnest(RES, data)
}

# Adding survival data
addsurvdat <- function(X, SurvMACEDf, SurvDMDf){
    X1 <- dplyr::transmute(
        X, sex, outcome = "MACE",
        data = purrr::map(data, dplyr::inner_join, SurvMACEDf, by = "eid"),
        data = purrr::map(data, filter, CHD == 0, Stroke == 0, PAD == 0)
    )
    X2 <- dplyr::transmute(
        X, sex, outcome = "DM",
        data = purrr::map(data, dplyr::inner_join, SurvDMDf, by = "eid"),
        data = purrr::map(data, filter, T2D == 0, T1D == 0, Insulin == 0, AntiDM == 0)
    )
    RES <- bind_rows(X1, X2)
    dplyr::mutate(
        RES,
        data = purrr::map(data, dplyr::select, -dplyr::where(\(x) all(x == 0)))
    )
}

# Rates of outcome by cluster
ratesclusfx <- function(X){
    RES <- dplyr::mutate(
        X,
        data = purrr::map(
            data,
            function(DAT){
                DAT <- tidyr::pivot_longer(
                    DAT,
                    dplyr::starts_with("prob"),
                    names_to = "Cluster",
                    values_to = "prob"
                )
                DAT <- dplyr::group_by(DAT, Cluster)
                dplyr::summarise(
                    DAT,
                    Ncases = sum(outcome_value * prob),
                    TPT = sum(outcome_timeyrs * prob)
                )
            }
        )
    )
    tidyr::unnest(RES, data)
}

# Rates of outcome by cluster and medication
ratesclusmedfx <- function(X){
    RES <- dplyr::mutate(
        X,
        data = purrr::map(
            data,
            function(DAT){
                DAT <- dplyr::mutate(
                    DAT, 
                    NoMed = 1 * (rowSums(dplyr::pick(dplyr::any_of(MEDICATION))) == 0)
                )
                DAT <- tidyr::pivot_longer(
                    DAT,
                    dplyr::starts_with("prob"),
                    names_to = "Cluster",
                    values_to = "prob"
                )
                DAT <- tidyr::pivot_longer(
                    DAT,
                    c(NoMed, dplyr::any_of(MEDICATION)),
                    names_to = "Med_name",
                    values_to = "Med_value"
                )
                DAT <- dplyr::group_by(DAT, Cluster, Med_name)
                dplyr::summarise(
                    DAT,
                    Ncases = sum(Med_value * outcome_value * prob),
                    TPT = sum(Med_value * outcome_timeyrs * prob),
                    .groups = "drop"
                )
            }
        )
    )
    tidyr::unnest(RES, data)
}

# Kaplan-Meier estimates
kmestfx <- function(X){
    RES <- dplyr::mutate(
        X,
        data = purrr::map(
            data,
            function(DAT){
                RES1 <- survival:::survfit.formula(
                    survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                    data = DAT
                )
                RES1 <- survival:::summary.survfit(RES1, times = 10)
                RES1 <- RES1[c("surv", "upper", "lower")]
                RES1 <- 1 - data.frame(RES1)
                colnames(RES1) <- c("risk", "lower", "upper")
                RES1 <- cbind(Cluster = "Overall", RES1)
                RES2 <- tibble::tibble(Cluster = names(DAT)[grepl("^prob", names(DAT))])
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
                        D = DAT
                    )
                )
                RES2 <- tidyr::unnest(RES2, RES3)
                dplyr::bind_rows(RES1, RES2)
            }
        )
    )
    tidyr::unnest(RES, data)
}

# Running Cox models
coxmodels <- function(X){
    X1 <- dplyr::filter(X, outcome == "MACE")
    X1 <- dplyr::mutate(
        X1,
        survdf = purrr::map2(
            data, sex,
            function(DATA, SEX){
                MODELDAT <- dplyr::transmute(
                    DATA,
                    eid,
                    outcome_value,
                    outcome_timeyrs,
                    age,
                    smoking,
                    bmi,
                    sbp,
                    tchol = hdl + ldl + (tg/2.2),
                    hdl,
                    age_smoking = age*smoking,
                    age_sbp = age*sbp,
                    age_tchol = age*tchol,
                    age_hdl = age*hdl,
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
                    lnegfr = log(egfr + 0.01),
                    lnegfrsq = lnegfr^2,
                    age_fg = age*fg,
                    age_lnegfr = age*lnegfr,
                    whr, alt, crp,
                    dplyr::across(
                        any_of(
                            c("T2D", "T2Dage", "HT", "CKD", "LiverFailure", "RA", "T1D",
                              "Insulin", "AntiDM", "AntiHT", "LipidLower")
                        )
                    ),
                    dplyr::across(dplyr::starts_with("prob"))
                )
                if(any(names(MODELDAT) == "T2D")){
                    MODELDAT$age_t2d <- MODELDAT$age * MODELDAT$T2D
                }
                dplyr::select(
                    MODELDAT, 
                    -c(scr, alpha, kappa, creat_kappa, minkappa, maxkappa, 
                       minkappa_alpha, maxkappa_exp, age_term, sex_term, egfr)
                )
            }
        )
    )
    X2 <- dplyr::filter(X, outcome == "DM")
    X2 <- dplyr::mutate(
        X2,
        survdf = purrr::map(
            data,
            function(DATA){
                dplyr::transmute(
                    DATA,
                    eid,
                    outcome_timeyrs, outcome_value,
                    dplyr::across(
                        dplyr::any_of(
                            c(
                                BIOMARKERS,
                                BODYSIZEINDEX,
                                COVARIATES,
                                DISEASES[!DISEASES %in% c("T1D", "T2D")],
                                MEDICATION[!MEDICATION %in% c("Insulin", "AntiDM")]
                            )
                        )
                    ),
                    dplyr::across(dplyr::starts_with("prob"))
                )
            }
        )
    )
    NewX <- dplyr::bind_rows(X1, X2)
    dplyr::mutate(
        NewX,
        survdf = purrr::map(
            survdf, 
            function(DATA){
                DATA <- dplyr::mutate(
                    DATA,
                    dplyr::across(
                        dplyr::everything(),
                        function(x){ 
                            x[is.infinite(x)] <- NaN
                            return(x)
                        }
                    )
                )
                DATA[complete.cases(DATA),]
            }
        ),
        data = purrr::map2(
            data, survdf,
            function(D1, D2){
                ID <- dplyr::select(D2, eid)
                dplyr::inner_join(ID, D1, by = "eid")
            }
        ),
        survdf = purrr::map(survdf, select, -eid),
        survdf = purrr::map(
            survdf,
            function(DATA){
                dplyr::mutate(
                    DATA,
                    dplyr::across(dplyr::starts_with("prob"), \(P) log(P/probBC), .names = "alr{.col}"),
                    alrprobBC = NULL,
                    dplyr::across(dplyr::starts_with("prob"), \(P) NULL)
                )
            }
        ),
        NullMod = purrr::map(
            survdf, 
            ~survival::coxph(
                survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1, 
                data = .x
            )
        ),
        mod_base = purrr::map(
            survdf,
            ~{
                MODELDAT <- dplyr::select(.x, -dplyr::starts_with("alrprob"))
                MOD1 <- survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = MODELDAT)
                COEFS <- coef(MOD1)
                COEFSNA <- COEFS[is.na(COEFS)]
                MODELDAT <- dplyr::select(MODELDAT, -all_of(names(COEFSNA)))
                survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = MODELDAT)
            }
        ),
        mod_clus = purrr::map(
            survdf,
            ~{
                MOD1 <- survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = .x)
                COEFS <- coef(MOD1)
                COEFSNA <- COEFS[is.na(COEFS)]
                MODELDAT <- dplyr::select(.x, -all_of(names(COEFSNA)))
                survival::coxph(survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ ., data = MODELDAT)
            }
        )
    )
}

# Extracting coefficients of survival model
survcoefx <- function(X){
    RES <- dplyr::select(X, sex, outcome, starts_with("mod_"))
    RES <- tidyr::pivot_longer(
        RES, 
        starts_with("mod_"), 
        names_to = c(".value", "model"), 
        names_sep = "_"
    )
    RES <- dplyr::mutate(
        RES,
        mod = purrr::map(
            mod, 
            function(MOD){
                tibble::tibble(
                    estimates = list(coef(MOD)),
                    varcovmat = list(vcov(MOD)),
                    Means = list(MOD$means),
                    Afit = list(
                        survival:::agsurv(
                            MOD$y, 
                            model.matrix(MOD) - rep(MOD$means, each = nrow(MOD$y)), 
                            rep(1.0, nrow(MOD$y)), 
                            exp(MOD$linear.predictors), 
                            3, 3
                        )
                    )
                )
            }
        )
    )
    tidyr::unnest(RES, mod)
}

# Obtain cumulative hazard/survival probabilities with summary data generated
PredCumHaz <- function(AFit, Means, Coefs, VCoVMat, NewX){
    newx <- dplyr::select(NewX, -outcome_time)
    newx <- as.matrix(newx)
    newtime <- NewX$outcome_time
    newx <- newx - rep(Means, each = nrow(newx))
    newrisk <- c(exp(newx %*% Coefs))
    j1 <- approx(
        x = AFit$time, y = 1:length(AFit$time), xout = newtime, 
        method = "constant", f = 0, yleft = 0, yright = length(AFit$time)
    )
    j1 <- j1$y
    chaz <-c(0, AFit$cumhaz)[j1+1]
    pred <- chaz * newrisk
    varh <- c(0, cumsum(AFit$varhaz))[j1+1]
    xbar <- rbind(0, AFit$xbar)[j1+1,,drop=F]
    dt <- (chaz * newx) - xbar
    se <- sqrt(varh + rowSums((dt %*% VCoVMat) *dt)) * newrisk
    tibble(fit = drop(pred), se.fit = drop(se))
    # To convert to survival: 
    # pred <- exp(-pred)
    # se <- se * exp(-pred)
}

# Comparing survival models
comparemods <- function(X){
    dplyr::transmute(
        X,
        sex, outcome,
        LL0 = purrr::map_dbl(NullMod, logLik),
        LLBase = purrr::map_dbl(mod_base, logLik),
        NVBase = purrr::map_dbl(mod_base, \(MOD) sum(!is.na(MOD$coefficients))),
        LLBaseCl = purrr::map_dbl(mod_clus, logLik),
        NVBaseCl = purrr::map_dbl(mod_clus, \(MOD) sum(!is.na(MOD$coefficients))),
        LRTstat = 2 * abs(LLBaseCl - LLBase),
        LRTdf = NVBaseCl - NVBase,
        LRTp = stats::pchisq(
            q = LRTstat,
            df = LRTdf,
            lower.tail = FALSE
        ),
        AdeqInd = (LLBase - LL0) / (LLBaseCl - LL0),
        CBase = purrr::map_dbl(mod_base, \(MOD) MOD[["concordance"]][["concordance"]]),
        CseBase = purrr::map_dbl(mod_base, \(MOD) MOD[["concordance"]][["std"]]),
        CBaseCl = purrr::map_dbl(mod_clus, \(MOD) MOD[["concordance"]][["concordance"]]),
        CseBaseCl = purrr::map_dbl(mod_clus, \(MOD) MOD[["concordance"]][["std"]]),
        Cdat = purrr::map2(mod_clus, mod_base, survival::concordance),
        Cdiff = purrr::map_dbl(Cdat, \(MOD) c(1, -1) %*% coef(MOD)),
        Cdiffse = purrr::map_dbl(Cdat, \(MOD) sqrt(c(1, -1) %*% vcov(MOD) %*% c(1, -1))),
        Cdiffp = 2*pnorm(-abs(Cdiff / Cdiffse)),
        Cdat = NULL
    )
}

# Adequacy index by cluster
AdeqIndClusFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex, outcome,
        Cluster = purrr::map(data, ~names(.x)[grepl("^prob", names(.x))]),
        dplyr::across(
            c(mod_base, mod_clus),
            ~purrr::pmap(
                list(data, Cluster, survdf, .x),
                function(DAT, CLUS, MODDAT, MOD){
                    sapply(
                        DAT[CLUS],
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
                            abs(LLclus - LL0clus)
                        }
                    )
                }
            )
        ),        
        AdeqInd = purrr::map2(
            mod_base, mod_clus, 
            ~pmin(.x/.y, 1)
        ),
        mod_base = NULL, 
        mod_clus = NULL
    )
    tidyr::unnest(RES, c(Cluster, AdeqInd))
}

# Adequacy index by pre-test probability
AdeqIndByPreFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex, outcome,
        RES1 = purrr::pmap(
            list(survdf, mod_base, mod_clus),
            function(DAT, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DAT, outcome_timeyrs = 10)
                yhat1 <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                RES2 <- purrr::map(
                    seq(0, .25, .01),
                    function(thresh){
                        TDAT <- DAT[yhat1 >= thresh,]
                        MODNULL <- survival::coxph(
                            formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 1,
                            data = TDAT
                        )
                        LL0 <- as.numeric(logLik(MODNULL))
                        NEWMOD1 <- survival::coxph(
                            formula = as.formula(MOD1),
                            data = TDAT,
                            init = coef(MOD1), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        LLMOD1 <- as.numeric(logLik(NEWMOD1))
                        LRMOD1 <- LLMOD1 - LL0
                        NEWMOD2 <- survival::coxph(
                            formula = as.formula(MOD2),
                            data = TDAT,
                            init = coef(MOD2), 
                            control = survival::coxph.control(iter.max = 0)
                        )
                        LLMOD2 <- as.numeric(logLik(NEWMOD2))
                        LRMOD2 <- LLMOD2 - LL0
                        tibble::tibble(
                            threshold = thresh, 
                            AdeqInd = pmin(1, LRMOD1/LRMOD2),
                        )
                    }
                )
                bind_rows(RES2)
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

# Estimate risk at a given time
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

# Decision curve analysis
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

# DCA
DCurvFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex, outcome,
        RES1 = purrr::pmap(
            list(survdf, mod_base, mod_clus),
            function(DAT, MOD1, MOD2){
                NEWDAT <- mutate(DAT, outcome_timeyrs = 10)
                DAT$base <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DAT$clus <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                dca(
                    formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~ 
                    base + clus,
                    data = DAT,
                    thresholds = seq(0, .25, .01),
                    time = 10
                )
            }
        )
    )
    tidyr::unnest(RES, RES1)
}

# DCA by cluster
DCurvbyClFx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex, outcome,
        RES1 = purrr::pmap(
            list(data, survdf, mod_base, mod_clus),
            function(DAT, DATMOD, MOD1, MOD2){
                NEWDAT <- dplyr::mutate(DATMOD, outcome_timeyrs = 10)
                DATMOD$base <- 1 - survival:::predict.coxph(MOD1, NEWDAT, type = "survival")
                DATMOD$clus <- 1 - survival:::predict.coxph(MOD2, NEWDAT, type = "survival")
                RES2 <- tibble::tibble(
                    Cluster = names(DAT)[grepl("^prob", names(DAT))]
                )
                RES2 <- dplyr::mutate(
                    RES2,
                    RES3 = purrr::map(
                        Cluster,
                        function(CL){
                            dca(
                                formula = survival::Surv(time = outcome_timeyrs, event = outcome_value) ~
                                base + clus,
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

## Test interactions between cluster allocations and medications
interactmodfx <- function(X){
    RES <- dplyr::transmute(
        X,
        sex, outcome,
        intmod = map(
            data,
            function(DAT){
                DAT <- dplyr::transmute(
                    DAT,
                    outcome_timeyrs, outcome_value,
                    dplyr::across(all_of(c(BODYSIZEINDEX, COVARIATES, BIOMARKERS))),
                    dplyr::across(any_of(MEDICATION)),
                    dplyr::across(dplyr::starts_with("prob"), \(P) log(P/probBC), .names = "alr{.col}"),
                    alrprobBC = NULL,
                    dplyr::across(dplyr::starts_with("prob"), \(P) NULL)
                )
                MEDS <- colnames(dplyr::select(DAT, dplyr::any_of(MEDICATION)))
                CLUS <- colnames(dplyr::select(DAT, dplyr::starts_with("alrprob")))
                RES1 <- tibble::tibble(Med_name = MEDS)
                RES1 <- dplyr::mutate(
                    RES1,
                    mods = purrr::map(
                        Med_name,
                        function(MED){
                            modlist <- list(
                                OnlyClusMod = coxph(
                                    reformulate(
                                        response = "Surv(time = outcome_timeyrs, event = outcome_value)",
                                        termlabels = paste(CLUS, MED, sep = "*")
                                    ),
                                    data = DAT
                                ),
                                FullMod = coxph(
                                    reformulate(
                                        response = "Surv(time = outcome_timeyrs, event = outcome_value)",
                                        termlabels = c(
                                            paste(CLUS, MED, sep = "*"),
                                            BODYSIZEINDEX, COVARIATES, BIOMARKERS
                                        )
                                    ),
                                    data = DAT
                                )
                            )
                            modlist <- purrr::map(
                                modlist,
                                ~tibble::tibble(
                                    estimates = list(coef(.x)),
                                    varcovmat = list(vcov(.x)),
                                    Means = list(.x$means),
                                    Afit = list(
                                        survival:::agsurv(
                                            .x$y, 
                                            model.matrix(.x) - rep(.x$means, each = nrow(.x$y)), 
                                            rep(1.0, nrow(.x$y)), 
                                            exp(.x$linear.predictors), 
                                            3, 3
                                        )
                                    )
                                )
                            )
                            dplyr::bind_rows(modlist, .id = "model")    
                        }
                    )
                )
                tidyr::unnest(RES1, mods)
            }
        )
    )
    tidyr::unnest(RES, intmod)
}