#' Format coverage files from genome positions (chrom, pos) to residue positions (symbol, index)
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#'
#' @param coverage a coverage file with 3 columnms: chrom, pos, tenX (proportion of samples with 10X coverage)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @export


format_coverage = function(coverage){

  if(ncol(coverage) != 3) stop("Coverage must be a three column data.table: (chrom, pos, over_10)")

  coverage = setNames(coverage, c("chrom", "pos", "over_10"))

  data("grmap_reduced")

  residue_mapped = merge(grmap_reduced, coverage, by.x=c("chr", "pos"), by.y = c("chrom", "pos"))

  setNames(residue_mapped[, .(mean(over_10)), by=.(gene, index)],
           c("symbol", "protein_position", "over_10"))

}

