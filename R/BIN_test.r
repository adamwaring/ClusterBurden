#' Split data into bins and carry out a two-sample goodness-of-fit test
#' Calculate a p-value for positional differences of rare missense-variant residue positions between cases and controls.
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns an object of class htest or p.value depending on value of argument "pval=?"
#'
#' @param case_residues vector of case variant residue positions
#' @param control_residues vector of control variant residue positions
#' @param case_coverage optional coverage data for cases in format: data.table(protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(protein_position, over_10)
#' @param cov_threshold threshold at which to exclude a residue from the analysis (choose 0 to keep all residues)
#' @param pval return only p-value or return chi-squared test output?
#' @param method method to bin data either mann-wald or nbins
#' @param nbins number of bins to use if method == "nbins"
#' @param plot_resids should chi-squared residuals be plotted? Defaults to False
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details The function takes a vector of case and control missense-variant residue positions (aggregated over a protein-coding-region) as input and returns a p-value representing the significance of variant clustering within the gene.
#' The linear sequence of the protein is split into 'bins' and the the counts for variant within each bin for each cohort are used to construct a kx2 contigency table where k is the number of bins and 2 is for the two cohorts: cases and controls.
#' The binning method is either:
#' - "mann_wald" where the number of bins k is determined by the total number of observed variants n by the equation k ~ n^(2/5)
#' - "nbins" where the user selects a specific number of bins (reasonable values here would be ~10-20 bins)
#' Setting "plot_resids" to true allows the residuals for each cell in the kx2 contigency table to be plotted - this allows the user to determine which cells (protein regions) contribute towards the significance of the test.
#'
#' When coverage files are supplied then regions with a 10X coverage below "cov_threshold" (default=0.5) are exluded from the analysis. For the remaining regions, cell counts are adjusted by the reciprocal of the mean coverage across the bin.
#'
#' @export
#' @examples
#' # The essential inputs are case_residues and control_residues
#'
#' # Example 1: Simulated NULL data
#' # Bin by the Mann-Wald heuristic; n ^ (2/5) where n = length(case_residues) + length(control_residues)
#' # simulate case-control residue positions from the same distribution
#'
#' nresidues = 1000 # length of the protein
#' probs = rexp(nresidues)^2 # probability of a missense variant at each residue
#'
#' case_residues = sample(1:nresidues, 100, rep=T, probs)
#' control_residues = sample(1:nresidues, 100, rep=T, probs)
#'
#' BIN_test(case_residues, control_residues)
#'
#' # Example 2: Simulated DISEASE data
#' # simulate case-control residue positions from different distributions
#'
#' nresidues = 1000 # length of the protein
#' probs = rexp(nresidues)^2
#' case_probs = probs * rep(c(1, 3, 1), c(200, 200, 600))
#' control_probs = probs * rep(c(2, 1, 2), c(200, 200, 600))
#'
#' case_residues = sample(1:nresidues, 100, rep=T, case_probs)
#' control_residues = sample(1:nresidues, 100, rep=T, control_probs)
#'
#' plot_distribs(case_residues, control_residues)
#'
#' BIN_test(case_residues, control_residues, plot_resids = T)

BIN_test = function(case_residues, control_residues, case_coverage=NULL, control_coverage=NULL, cov_threshold=0.5, pval=T, method="mann-wald", nbins=NULL, plot_resids=F){


  if(!is.null(case_coverage)){

    case_blacklist = case_coverage[over_10 < cov_threshold, protein_position]

    if(length(case_blacklist) > 0){
      case_coverage = case_coverage[!protein_position %in% case_blacklist]
      control_coverage = control_coverage[!protein_position %in% case_blacklist]

      message(paste0("Some regions ignored due to low coverage in cases (threshold = ",
                     cov_threshold, ", residues = ", findRanges(case_blacklist, "str"), ")"))

      if(any(case_residues %in% case_blacklist)){
        index = which(case_residues %in% case_blacklist)
        if(length(index) > 0){
          case_residues = case_residues[-index]
          message(paste0("Removed ", length(index), " case variants in region blacklisted by case coverage"))
        }
      }

      if(any(control_residues %in% case_blacklist)){
        index2 = which(control_residues %in% case_blacklist)
        if(length(index2) > 0){
          control_residues = control_residues[-index2]
          message(paste0("Removed ", length(index2), " control variants in region blacklisted by case coverage"))
        }
      }
    }

  }
  if(!is.null(control_coverage)){

    control_blacklist = control_coverage[over_10 < cov_threshold, protein_position]

    if(length(control_blacklist) > 0){

      case_coverage = case_coverage[!protein_position %in% control_blacklist]
      control_coverage = control_coverage[!protein_position %in% control_blacklist]


      message(paste0("Some regions ignored due to low coverage in controls (threshold = ",
                     cov_threshold, ", residues = ", findRanges(control_blacklist, "str"), ")"))

      if(any(case_residues %in% control_blacklist)){
        index = which(case_residues %in% control_blacklist)
        if(length(index) > 0){
          case_residues = case_residues[-index]
          message(paste0("Removed ", length(index), " case variants in region blacklisted by control coverage"))
        }
      }

      if(any(control_residues %in% control_blacklist)){
        index2 = which(control_residues %in% control_blacklist)
        if(length(index2) > 0){
          control_residues = control_residues[-index2]
          message(paste0("Removed ", length(index2), " control variants in region blacklisted by control coverage"))
        }
      }

    }


  }

  # check input and restructure
  if(length(case_residues) == 0){
    message("Cases have 0 variants")
    return(NA_real_)
  }
  if(length(control_residues) == 0){
    message("Controls have 0 variants")
    return(NA_real_)
  }

  # check method type
  if(!method %in% c("mann-wald", "nbins")){
    method = "mann-wald"
    message(paste0("'method' arguement must equal'mann-wald' or 'nbins' but is ", method, ". Defaulting to mann-wald."))
  }

  # check the sfs for obvious outliers
  inspect_residues = function(residues, name){

    if(!is.vector(residues) & !is.numeric(residues)){
      stop(paste0(name, "_residues must be a numeric vector but is a ", class(residues)))
    }

    if(any(is.na(residues))){
      warning(paste0("NA values present in ", name, "_residues vector. Removing ", sum(is.na(residues)), " NAs from analysis."))
      residues = residues[!is.na(residues)]
    }

    counts = table(residues)

    most_freq_ac = names(which.max(table(counts)))
    if(most_freq_ac != "1") warning("Expected most variants to be singletons but the most frequent allele count is ", most_freq_ac)

    outlier_boundary = mean(counts) + 3*sd(counts) + 2
    outliers = counts[counts > outlier_boundary]

    if(length(outliers) > 0) warning("There are some outliers in the frequency some ", name, "_residues are observed. (residues = ", paste(names(outliers), collapse=", "), "; ac = ", paste(outliers, collapse=", "), "). Be careful as these variants may be founders, from related samples or common ancestry-specific variants (that have evaded frequency filtering). Such variants may cause artefactual clusters.")

    return(residues)

  }
  case_residues = inspect_residues(case_residues, "cases")
  control_residues = inspect_residues(control_residues, "controls")

  # params
  residues = c(case_residues, control_residues)
  n = length(residues)
  nresidues = max(residues)
  if(method == "nbins"){
    if(nbins > nresidues){
      warning("nbins > nresidues: nbins = nresidues")
      nbins = nresidues - 1
    }
    bins = nbins + 1
  } else bins = round(n ^ (2/5) + 1)


  # count variants across each bin
  breaks = seq(1, nresidues, length.out=bins)
  case_tab = table(cut(case_residues, breaks, include.lowest = T))
  control_tab = table(cut(control_residues, breaks, include.lowest = T))

  # remove bins with zero counts in both cases and controls
  contig = rbind(case_tab, control_tab)
  index = which(apply(contig, 2, function(x) sum(x) == 0))
  if(length(index) > 0) contig = contig[, -index]

  # get mean coverage across bins
  if(!is.null(case_coverage)){

    case_coverage[,bin:=cut(protein_position, breaks, include.lowest = T)]
    case_weight = case_coverage[!is.na(bin),mean(over_10), by=bin]$V1
    contig[1,] = contig[1,] / case_weight

  }
  if(!is.null(control_coverage)){

    control_coverage[,bin:=cut(protein_position, breaks, include.lowest = T)]
    control_weight = control_coverage[!is.na(bin),mean(over_10), by=bin]$V1
    contig[2,] = contig[2,] / control_weight

  }


  # calculate chi-squared two-sample test on contingency table
  chisq = chisq.test(contig)


  # plot the resulting residuals per bin to identify bins that contributed the most to test statistic
  if(plot_resids == T) print(plot_residuals(chisq))


  if(pval) return(chisq$p.value) else return(chisq)

}


#
#  nresidues = 1000 # length of the protein
#  probs = rexp(nresidues)^2 # probability of a missense variant at each residue
#
#  case_residues = sample(1:nresidues, 100, rep=T, probs)
#  control_residues = sample(1:nresidues, 100, rep=T, probs)
#
#  BIN_test(case_residues, control_residues)
#
#  # Example 2: Simulated DISEASE data
#  # simulate case-control residue positions from different distributions
#
#  nresidues = 1000 # length of the protein
#  probs = rexp(nresidues)^2
#  case_probs = probs * rep(c(1, 3, 1), c(200, 200, 600))
#  control_probs = probs * rep(c(2, 1, 2), c(200, 200, 600))
#
#  case_residues = sample(1:nresidues, 100, rep=T, case_probs)
#  control_residues = sample(1:nresidues, 100, rep=T, control_probs)
#
#  plot_distribs(case_residues, control_residues)
#
#  BIN_test(case_residues, control_residues)
#
#  case_coverage = rep(c(0, 1), c(200, 800))
#  control_coverage = rep(c(1, 0.8), c(800, 200))
#  case_coverage = data.table(protein_position=1:nresidues, over_10=case_coverage)
#  control_coverage = data.table(protein_position=1:nresidues, over_10=control_coverage)
#
#  BIN_test(case_residues, control_residues, case_coverage, control_coverage)
