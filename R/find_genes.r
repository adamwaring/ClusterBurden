#' Select a subset of genes to use in an WES analysis
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#'
#' @param n_genes how many random null genes?
#' @param clndn_regex regular expression to find disease genes from ClinVAR
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @export



find_genes = function(n_genes, clndn_regex){

  data("canon_txs")

  genes = canon_txs[,sample(gene, n_genes)]

  data("clinvar_cases")

  genes = c(genes, clinvar_cases[grepl("[Pp]ath", clnsig) & grepl(paste(clndn_regex, collapse="|"), clndn, ignore.case = T), unique(symbol)])

  if(any(genes=="TTN")) genes = genes[-which(genes=="TTN")]

  return(genes)
}
