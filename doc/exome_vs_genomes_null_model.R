## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7
)

## -----------------------------------------------------------------------------
library(ClusterBurden)

controls = collect_gnomad_controls(dataset="exome")
cases = collect_gnomad_controls(dataset="genome")
cases[,aff:=1]

attributes(controls)$provenance = attributes(cases)$provenance = NULL

binpval = BIN_test_WES(cases, controls)

head(binpval)



## -----------------------------------------------------------------------------
# nominal significance 
binpval[,sum(BIN.test_pvalue < 0.05, na.rm=T)/.N]

# Bonferonni significance 
binpval[,sum(BIN.test_pvalue < 0.05/.N)/.N]

## -----------------------------------------------------------------------------
# function to generate lambda inflation metric
lambda = function(p) median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

# lambda across all non-NA p-values
lambda(binpval[!is.na(BIN.test_pvalue), BIN.test_pvalue])

# lambda across p-values not in giant proteins i.e. > 2000 residues
lambda(binpval[!is.na(BIN.test_pvalue) & !grepl("high", pl_flag), BIN.test_pvalue])

## -----------------------------------------------------------------------------
manhattan(binpval, "bin-test", 10, SCALE=0.5)

## -----------------------------------------------------------------------------
attributes(controls)$provenance = attributes(cases)$provenance = "auto"

binpval = BIN_test_WES(cases, controls)

binpval

## -----------------------------------------------------------------------------
# nominal significance 
binpval[,sum(BIN.test_pvalue < 0.05, na.rm=T)/.N]

# Bonferonni significance 
binpval[,sum(BIN.test_pvalue < 0.05/.N, na.rm=T)/.N]

# function to generate lambda inflation metric
lambda = function(p) median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

# lambda across all non-NA p-values
lambda(binpval[!is.na(BIN.test_pvalue), BIN.test_pvalue])

# manhattan
manhattan(binpval, "bin-test", 10, SCALE=0.5)

## -----------------------------------------------------------------------------
binpval = binpval[!grepl("high", cov_flag1) & !grepl("high", cov_flag2) & !grepl("high", pl_flag)]

# nominal significance 
binpval[,sum(BIN.test_pvalue < 0.05, na.rm=T)/.N]

# Bonferonni significance 
binpval[,sum(BIN.test_pvalue < 0.05/.N, na.rm=T)/.N]

# function to generate lambda inflation metric
lambda = function(p) median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

# lambda across all non-NA p-values
lambda(binpval[!is.na(BIN.test_pvalue), BIN.test_pvalue])

# manhattan
manhattan(binpval, "bin-test", 10, SCALE=0.5)

