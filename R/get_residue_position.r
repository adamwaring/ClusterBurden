#' Extract amino acid residue position from HGVS protein consequence
#'
#' This function makes most sense for mutations that impact a single amino acid such as the missense mutation p.Trp26Cys.
#' For indels and duplications, where many amino acid positions are impacted, the first amino acid, i.e. the lowest index, is taken.
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a numeric vector of amino acid positions
#'
#' @param hgvs vector of variant protein consequences in HGVS format e.g. p.Trp26Cys or p.A54Y.
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @export
#' @examples
#' get_residue_position("p.Trp26Cys")


get_residue_position = function(hgvs){

  matches = regmatches(hgvs, gregexpr("[[:digit:]]+", hgvs))

  no_matches = sapply(matches, function(x) length(x) == 0)

  if(sum(no_matches) > 0){
    no_matches_input = paste0(which(no_matches), " (", hgvs[no_matches], ")")
    warning("Could not extract protein position from input: ", paste(no_matches_input, collapse=", "), ". Input must be in HGVS format. E.g. p.Trp26Cys.")
    matches[no_matches] = NA
  }

  multi_matches = sapply(matches, function(x) length(x) > 1)

  if(sum(multi_matches) > 0){
    multi_matches_input = paste0(which(multi_matches), " (", hgvs[multi_matches], ")")
    first_matches = sapply(matches[multi_matches], '[[', 1)
    warning("Multiple possible positions found in ", paste(multi_matches_input, collapse=", "), ". Taking first matches: ", paste(first_matches, collapse=", "),".")
    matches[multi_matches] = first_matches
  }

  return(as.numeric(unlist(matches)))

}

