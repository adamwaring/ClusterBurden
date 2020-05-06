#' ClusterBurden test for burden and position of rare-missense variants
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return A data.table with columns: symbol, ncases, ncontrols, case_ss, control_ss, binp, burdp, clusterburden
#'
#' @param cases case data in format: data.table(aff, symbol, protein_position, ac)
#' @param controls control data in format: data.table(aff, symbol, protein_position, ac)
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cases_ss sample sizes for cases either a scalar for all genes or a data.table with two columns: symbol, ss
#' @param controls_ss sample sizes for controls either a scalar for all genes or a data.table with two columns: symbol, ss
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param covstats include additional coverage information as a column in results?
#' @param messages print messages to the terminal (e.g. how many variants removed by coverage)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Perform ClusterBurden test across an exome set
#'
#' @export


ClusterBurden_WES = function(cases, controls, case_coverage=NULL, control_coverage=NULL, cases_ss=NULL, controls_ss=NULL, cov_threshold=0.5, alternative="greater", covstats=F, messages=T){


  binpval = BIN_test_WES(cases, controls, case_coverage, control_coverage, cov_threshold, covstats, messages)

  burdenpval = burden_test_WES(cases, controls, cases_ss, controls_ss, case_coverage, control_coverage, cov_threshold, alternative, covstats, messages)

  pvals = merge(binpval, burdenpval, by="symbol", all = T)

  pvals[,ClusterBurden:=mapply(function(x, y) combine_ps(x, y), BIN.test_pvalue, burden_pvalue)]

  return(pvals)


}
