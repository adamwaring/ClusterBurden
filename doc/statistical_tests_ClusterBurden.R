## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7
)

## -----------------------------------------------------------------------------
 library(ClusterBurden)
 nresidues = 1000 
 probs = rexp(nresidues)^2 

 case_residues = sample(1:nresidues, 100, rep=T, probs)
 control_residues = sample(1:nresidues, 100, rep=T, probs)

## -----------------------------------------------------------------------------
BIN_test(case_residues, control_residues)

## -----------------------------------------------------------------------------
 case_probs = probs * rep(c(1, 3, 1), c(200, 200, 600))
 control_probs = probs * rep(c(1, 0.5, 1), c(200, 200, 600))

 case_residues = sample(1:nresidues, 100, rep=T, case_probs)
 control_residues = sample(1:nresidues, 100, rep=T, control_probs)

 plot_distribs(case_residues, control_residues)

## ----fig.height=7-------------------------------------------------------------
BIN_test(case_residues, control_residues)

## -----------------------------------------------------------------------------
  case_coverage = rep(c(0, 1), c(200, 800))
  control_coverage = rep(c(1, 0.8), c(800, 200))

  case_probs = probs * case_coverage
  control_probs = probs * control_coverage
  
  case_residues = sample(1:nresidues, 100, rep=T, case_probs)
  control_residues = sample(1:nresidues, 100, rep=T, control_probs)
  
  BIN_test(case_residues, control_residues)
 

## -----------------------------------------------------------------------------
  case_coverage = data.table(protein_position=1:nresidues, over_10=case_coverage)
  control_coverage = data.table(protein_position=1:nresidues, over_10=control_coverage)
  
  BIN_test(case_residues, control_residues, case_coverage, control_coverage)

## -----------------------------------------------------------------------------

  binpval = BIN_test(case_residues, control_residues, case_coverage, control_coverage)
  burdp = burden_test(n1 = 80, n2 = 40, ss1 = 1000, ss2 = 1000) # 2X excess in cases
  
  combine_ps(burdp, binpval)

