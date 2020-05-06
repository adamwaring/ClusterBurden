#' Combine p-values
#'
#' Combine p-values using chi-squared or z-score combination method.
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a single unnamed p-value of type num
#'
#' @param burden_pval P-value from a burden test such as Fisher's exact-test
#' @param position_pval P-value from a position test such as Chi-squared two-sample test
#' @param method Method to combine p-values. "chisq" for Fisher's method. "z-score" for Stouffer's method.
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details This is a method to combine p-values using either Fisher's method or Stouffer's method.
#'
#' Usage is simple however you must first calculate burden and position p-values.
#' Burden p-values can be calculated using ClusterBurden::burden_test() or any other RVAT that does not use positional information
#' Position p-values can be calculated using ClusterBurden::BIN_test() or any other position test that is not sensitive to burden
#'
#' If any p-values are zero then the function returns zero.
#'
#' The functions have been adapted from the R package metap but re-written here to reduce dependencies.
#'
#' @export


combine_ps = function(burden_pval, position_pval, method="chisq"){

  pvals = c(position_pval, burden_pval)

  if(sum(is.na(pvals))==2) return(NA_real_)
  if(sum(is.na(pvals))==1) return(pvals[!is.na(pvals)])
  if(any(pvals == 0)) return(0)



  if(method == "chisq"){
    # take the natural log of the pvalues
    # multiply by -2 to give the chisq
    # degrees of freedom is 2 * the number of pvalues
    # use pchisq to get the pvalues with lower.tail = F
    chisq = -2 * sum(log(pvals))
    pchisq(chisq, 4, lower.tail = F)

  } else if(method == "z-score"){
    # get qnorm
    # matrix multiply by weights
    # divide by the sum of squared weights to get z-score
    # calculate the normal p
    z = (qnorm(pvals, lower.tail = FALSE) %*% weights)/sqrt(sum(weights^2))
    pnorm(z, lower.tail = FALSE)

  }

}

