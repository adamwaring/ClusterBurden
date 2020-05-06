#' Carry out a Fisher's exact test
#'
#' Burden test with formatting step to make a contingency table
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns an object of class htest - use $p for p-value
#'
#' @param n1 number of carriers cases
#' @param n2 number of carriers controls
#' @param ss1 sample size cases
#' @param ss2 sample size controls
#' @param pval return p-value or full output?
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Fisher.test from base R with extra formatting steps.
#'
#' For a one-sided test i.e. excess in cases, use the default alternative="greater"
#'
#' if pval=F
#' Use $p for p-value
#' Use $est for odds-ratio
#' Use $conf.int[1:2] for 95% CI
#'
#' @export
#' @examples
#' # For the simple case where you have the vector status (binary case or control) and a vector of whether they
#' have a variant 'car' (presumably after filtering) then;
#'
#' N = 100
#' car = sample(c(0, 1), N, rep=T)
#' status = sample(c(0, 1), N, rep=T)
#'
#' burden_test(car, status)
#'
#' burden_test(n1=10, n2=20, ss1=100, ss2=100)

# 1. control for uneven coverage

burden_test = function(n1=NULL, n2=NULL, ss1=NULL, ss2=NULL, pval=T, case_coverage=NULL, control_coverage=NULL, cov_threshold=0.5, alternative="greater"){

  if(!any(!is.null(n1) & !is.null(n2) & !is.null(ss1) & !is.null(ss2))){

    stop("Must supply n1, n2, ss1 and ss2")

  }

  if(!is.null(case_coverage)){

    ss1 = ss1 * mean(case_coverage$over_10)

  }

  if(!is.null(control_coverage)){

    ss2 = ss2 * mean(control_coverage$over_10)

  }


  if(n1 > ss1) n1 = ss1
  if(n2 > ss2) n2 = ss2

  if(n1==0) n1=0.5
  if(n2==0) n2=0.5

  contig = rbind(c(n1, ss1-n1), c(n2, ss2-n2))

  fisher_output = suppressWarnings(fisher.test(contig, alternative = alternative))

  if(pval) return(fisher_output$p) else return(fisher_output)
}
