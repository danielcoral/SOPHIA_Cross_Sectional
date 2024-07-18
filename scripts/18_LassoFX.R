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


# Adding covariate data
addcovardat <- function(X, CovarDat){
    dplyr::mutate(
        X,
        data = purrr::map(data, inner_join, CovarDat, by = "eid")
    )
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
    RES <- dplyr::bind_rows(X1, X2)
    dplyr::mutate(
        RES,
        data = purrr::map(data, dplyr::select, -dplyr::where(\(x) all(x == 0)))
    )
}

# Follow-up subsets
futsubsetsfx <- function(X){
    RES <- tidyr::expand_grid(
        X,
        fut = c(5, 10)
    )
    RES <- dplyr::mutate(
        RES,
        data = purrr::map2(
            data, fut,
            function(DAT, FUT){
                dplyr::mutate(
                    DAT,
                    outcome_value = outcome_value * (outcome_timeyrs <= FUT),
                    outcome_timeyrs = pmin(outcome_timeyrs, FUT)
                )
            }
        ),
        maxfut = map_dbl(data, ~max(.x$outcome_timeyrs)),
        valid = maxfut > (fut - 1)
    )
    RES <- dplyr::filter(RES, valid)
    dplyr::select(RES, -c(valid, maxfut))
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
                    dplyr::across(dplyr::starts_with("prob"), \(x) ifelse(x < 1e-15, 1e-15, x)),
                    dplyr::across(dplyr::starts_with("prob"), \(x) ifelse(x > (1 - 1e-15), (1 - 1e-15), x)),
                    dplyr::across(dplyr::starts_with("prob"), \(P) log(P/probBC), .names = "alr{.col}"),
                    alrprobBC = NULL,
                    dplyr::across(dplyr::starts_with("prob"), \(P) NULL)
                )
            }
        ),
        mod_clus = purrr::map(
            survdf,
            ~{
                X <- as.matrix(dplyr::select(.x, -c(outcome_value, outcome_timeyrs)))
                Y <- with(.x, survival::Surv(time = outcome_timeyrs, event = outcome_value))
                cvfit <- glmnet::cv.glmnet(x = X, y = Y, alpha = 1, family = "cox")
                coef(cvfit, s = "lambda.min")
            }
        )
    )
}
