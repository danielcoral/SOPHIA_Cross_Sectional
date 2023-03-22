# Functions

# ==============================================================================
#
# DESCRIPTIVE OVERVIEW OF NUMERIC AND FACTOR VARIABLES
# WITH OVERALL TESTING OF DIFFERENCES BETWEEN BMI CATEGORIES
#
# Author: Carl Delfin (CRDF)
# Version: 3
# Date: 2022-10-28
#
# Instructions:
#
# (1) You must have the openxlsx package in order to save output as Excel sheets:
#     install.packages("openxlsx", dependencies = TRUE)
#
# (2) You must have the data loaded in a data frame named 'df' before running
#     the script, and the data must be processed according to the variable 
#     format table provided.
#
# (3) You must set a path where the output will be saved, edit the 'path'
#     variable below so that it matches your system.
#
# (4) If you want to have output by groups, you can add as many grouping
#     variables as you want in the 'group_vars' vector below. You will get
#     one output for each grouping variable.
#
# (5) OPTIONAL: If you want to filter out cases where the smallest group or
#     count is less than 5, set 'filter_5' to TRUE below.
#
# (6) When the above is done, just run the script, e.g. by pressing the
#     'Source' button in RStudio.
#
# ==============================================================================
df <- recoded_dat
df$sex <- as.factor(df$sex)
df$smoking <- as.factor(df$smoking)
# change as necessary:
path <- ""

# change as necessary
filter_5 <- FALSE 

# WOMEN

df <- filter(df, sex == "Female")
# ------------------------------------------------------------------------------
cat("Creating necessary functions...\n")
# ------------------------------------------------------------------------------

# get descriptives from numeric variables:

get_numeric_desc <- function(variable) {
  tmp <- df[, variable]
  out <- data.frame(
    variable = variable,
    n_total = length(tmp),
    n_missing = sum(is.na(tmp)),
    n_complete = sum(!is.na(tmp)),
    perc_complete = round((sum(!is.na(tmp)) / length(tmp)) * 100, 2),
    min = round(min(tmp, na.rm = TRUE), 2),
    max = round(max(tmp, na.rm = TRUE), 2),
    mean = round(mean(tmp, na.rm = TRUE), 2),
    sd = round(sd(tmp, na.rm = TRUE), 2),
    se = round(sd(tmp, na.rm = TRUE) / sqrt(sum(!is.na(tmp))), 2),
    median = round(median(tmp, na.rm = TRUE), 2),
    q25 = round(quantile(tmp, probs = 0.25, na.rm = TRUE), 2),
    q75 = round(quantile(tmp, probs = 0.75, na.rm = TRUE), 2))
  
  rownames(out) <- NULL
  return(out)
}


# descriptives from factor variables:

get_factor_desc <- function(variable) {
  
  res <- NULL
  out <- NULL
  
  for (i in levels(df[, variable])) {
    tmp <- subset(df, eval(parse(text = variable)) == i)
    tmp <- tmp[, variable]
    
    res <- data.frame(
      variable = variable,
      n_complete = sum(!is.na(df[, variable])),
      variable_level = i,
      count = length(tmp))
    
    rownames(res) <- NULL
    out <- rbind(res, out)
  }
  
  out$perc_of_n_complete <- round((out$count / sum(!is.na(df[, variable]))) * 100, 2)
  return(out)
}


# ------------------------------------------------------------------------------
cat("Creating descriptive summaries for numeric variables...\n")
# ------------------------------------------------------------------------------

# all numeric variables in 'df'
numeric_desc <- do.call(rbind.data.frame, 
                        lapply(colnames(df[sapply(df, is.numeric)]), 
                               get_numeric_desc))

openxlsx::write.xlsx(numeric_desc,
                     paste0(path, "numeric_desc_women.xlsx"))


# ------------------------------------------------------------------------------
cat("Creating descriptive summaries for factor variables...\n")
# ------------------------------------------------------------------------------

# all factor variables in 'df'
factor_desc <- do.call(rbind.data.frame, 
                       lapply(colnames(df[sapply(df, is.factor)]), 
                              get_factor_desc))
# filter?
if (filter_5 == TRUE) {
  factor_desc <- factor_desc[factor_desc$count >= 5, ]
}

openxlsx::write.xlsx(factor_desc,
                     paste0(path, "factor_desc_women.xlsx"))


# ------------------------------------------------------------------------------

df <- recoded_dat
df$sex <- as.factor(df$sex)
df$smoking <- as.factor(df$smoking)
# change as necessary:
path <- ""

# change as necessary
filter_5 <- FALSE 

df <- filter(df, sex == "Male")

# get descriptives from numeric variables:

get_numeric_desc <- function(variable) {
  tmp <- df[, variable]
  out <- data.frame(
    variable = variable,
    n_total = length(tmp),
    n_missing = sum(is.na(tmp)),
    n_complete = sum(!is.na(tmp)),
    perc_complete = round((sum(!is.na(tmp)) / length(tmp)) * 100, 2),
    min = round(min(tmp, na.rm = TRUE), 2),
    max = round(max(tmp, na.rm = TRUE), 2),
    mean = round(mean(tmp, na.rm = TRUE), 2),
    sd = round(sd(tmp, na.rm = TRUE), 2),
    se = round(sd(tmp, na.rm = TRUE) / sqrt(sum(!is.na(tmp))), 2),
    median = round(median(tmp, na.rm = TRUE), 2),
    q25 = round(quantile(tmp, probs = 0.25, na.rm = TRUE), 2),
    q75 = round(quantile(tmp, probs = 0.75, na.rm = TRUE), 2))
  
  rownames(out) <- NULL
  return(out)
}


# descriptives from factor variables:

get_factor_desc <- function(variable) {
  
  res <- NULL
  out <- NULL
  
  for (i in levels(df[, variable])) {
    tmp <- subset(df, eval(parse(text = variable)) == i)
    tmp <- tmp[, variable]
    
    res <- data.frame(
      variable = variable,
      n_complete = sum(!is.na(df[, variable])),
      variable_level = i,
      count = length(tmp))
    
    rownames(res) <- NULL
    out <- rbind(res, out)
  }
  
  out$perc_of_n_complete <- round((out$count / sum(!is.na(df[, variable]))) * 100, 2)
  return(out)
}


# ------------------------------------------------------------------------------
cat("Creating descriptive summaries for numeric variables...\n")
# ------------------------------------------------------------------------------

# all numeric variables in 'df'
numeric_desc <- do.call(rbind.data.frame, 
                        lapply(colnames(df[sapply(df, is.numeric)]), 
                               get_numeric_desc))

openxlsx::write.xlsx(numeric_desc,
                     paste0(path, "numeric_desc_men.xlsx"))


# ------------------------------------------------------------------------------
cat("Creating descriptive summaries for factor variables...\n")
# ------------------------------------------------------------------------------

# all factor variables in 'df'
factor_desc <- do.call(rbind.data.frame, 
                       lapply(colnames(df[sapply(df, is.factor)]), 
                              get_factor_desc))
# filter?
if (filter_5 == TRUE) {
  factor_desc <- factor_desc[factor_desc$count >= 5, ]
}

openxlsx::write.xlsx(factor_desc,
                     paste0(path, "factor_desc_men.xlsx"))


# ------------------------------------------------------------------------------

