#' Exclude residue positions from a dataset based on coverage files
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#'
#' @param dataset case-control data.table with at least columns: symbol and protein_position
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, index, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, index, over_10)
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param messages print messages to the terminal (e.g. how many variants removed by coverage)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Variant residues are removed from both cases and controls if they have low coverage in either case or control coverage. These positions are also removed from the coverage files to make downstream esimation of regional coverage accurate in light of variants removed.
#' @export


remove_uncovered_residues = function(dataset, case_coverage, control_coverage, cov_threshold, messages){

  setkey(dataset, symbol, protein_position)
  N = nrow(dataset)

  case_cov_raw = case_coverage
  control_cov_raw = control_coverage

  null1 = !is.null(case_coverage)
  null2 = !is.null(control_coverage)

  if(null1) setkey(case_coverage, symbol, protein_position)
  if(null2) setkey(control_coverage, symbol, protein_position)

  if(null1 & null2){

    case_coverage = case_coverage[control_coverage, nomatch=0]
    control_coverage = control_coverage[case_coverage, nomatch=0]

  }

  # case residues not covered
  if(null1){

    case_include = unique(case_coverage[over_10 >= cov_threshold, .(symbol, protein_position)])

    if(nrow(case_include) == 0) stop("No case variants after filtering for coverage")

    case_coverage = case_coverage[case_include, nomatch=0]

    control_coverage = control_coverage[case_include, nomatch=0]

    dataset = dataset[case_include, nomatch=0]

  }

  # control residues not covered
  if(null2){

    control_include = unique(control_coverage[over_10 >= cov_threshold, .(symbol, protein_position)])

    if(nrow(control_include) == 0) stop("No case variants after filtering for coverage")

    control_coverage = control_coverage[control_include, nomatch=0]

    case_coverage = case_coverage[control_include, nomatch=0]

    dataset = dataset[control_include, nomatch=0]

  }

  if(messages & (null1 | null2)) message(paste0(N-nrow(dataset), " variants removed with no coverage or less than ",
                  cov_threshold, " 10X coverage in either case_coverage or control_coverage"))

  blacklist = NULL
  if(null1) blacklist = rbind(blacklist, case_cov_raw[!case_coverage])
  if(null2) blacklist = rbind(blacklist, control_cov_raw[!control_coverage])

  if(!is.null(blacklist)){
    blacklist=unique(blacklist[,.(symbol, protein_position)])
    blacklist = blacklist[symbol=="A1BG",.(ranges=findRanges(protein_position, "str")), by=symbol]
  }

  return(list(dataset, case_coverage, control_coverage, blacklist))

}



