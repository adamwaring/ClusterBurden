#' Make example cases based on ClinVAR and gnomAD data
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return A data.table of example cases
#'
#' @param matchcontrols control group generated by collect_gnomad_controls()
#' @param clndn_regex ClinVAR disease annotation regular expression (require some knowledge of ClinVAR:clndn)
#' @param inframes include inframes as missense variants?
#' @param max_inframe_size maximum inframe size to include to a maximum of 3 (default=3)
#' @param filtertype type of allele frequency to filter MAF; "global", "popmax" or "strict" where global is the total allele frequency in controls, popmax is the maximum allele frequency in any ancestry group in GnomAD excluding ASJ, FIN, OTH and strict is the maximum allele frequency in any ancestry group in gnomAD or globally across GnomAD genomes
#' @param maxmaf maf threshold using for frequency filtering default is 0.0001 (i.e. 0.1%)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details This semi-synthetic data for example purposes. The resulting dataset is a composite of GnomAD gene variants and ClinVAR gene variants.
#' Disease genes are selected from ClinVAR according to the argument "clndn_regex" and null genes are selected from whichever GnomAD dataset is not used in the "matchcontrols" argument.
#' Allele counts and sample sizes are not available for ClinVAR so these are simulated: allele counts from an exponential distribution with rate=1.2 (resulting in
#' mostly singletons) and sample sizes to produce gene-level ORs x1:2 against gnomAD exomes.
#'
#' @export
#' @examples
#' # Filtering for cases must match filtering for controls
#' # Controls must be generated by collect_gnomad_controls()
#' # Genes for each function must be calculated first using find_genes()
#'
#' n_random_genes = 200
#' disease = "cardiomyopathy"
#'
#' genes = find_genes(n_random_genes, disease)
#'
#' controls = collect_gnomad_controls(genenames=genes)
#' cases = collect_example_cases(controls, disease)


collect_example_cases = function(matchcontrols, clndn_regex, inframes=T, max_inframe_size=3, filtertype="strict", maxmaf=0.0001){

  data("clinvar_cases")

  genes = clinvar_cases[symbol != "TTN" & cases[[filtertype]] < maxmaf & n_res <= max_inframe_size & grepl("[Pp]ath", clnsig) & grepl(paste(clndn_regex, collapse="|"), clndn, ignore.case = T), unique(symbol)]

  cases = clinvar_cases[symbol %in% genes]

  cases = cases[,ac:=ceiling(rexp(.N, 1.2))]

  cases[,group:="cln"]

  ctrlgenes = attributes(matchcontrols)$coverage[,.(group=unique(group)), by=symbol]
  ctrlgenes = ctrlgenes[!symbol %in% cases$symbol]


  cases_ctrls = NULL
  coverage = NULL

  if(any(ctrlgenes$group == "e")){

    data("gcontrols")
    gcontrols = gcontrols[af_g2 <= maxmaf]
    cases_ctrls = rbind(cases_ctrls, cbind(gcontrols[symbol%in%ctrlgenes[group=="e", symbol]], group="g"))

    data("genome_cov")
    coverage = rbind(coverage, genome_cov[symbol%in%ctrlgenes[group=="e", symbol], .(symbol, protein_position, over_10, group="g")])

  }

  if(any(ctrlgenes$group == "g")){

    data("econtrols")
    econtrols = econtrols[global <= maxmaf]
    cases_ctrls = rbind(cases_ctrls, cbind(econtrols[symbol%in%ctrlgenes[group=="g", symbol]], group="e"), fill=T)

    data("exome_cov")
    coverage = rbind(coverage, exome_cov[symbol%in%ctrlgenes[group=="g", symbol], .(symbol, protein_position, over_10, group="e")])
  }


  cases = rbind(cases, cases_ctrls, fill=T)

  # need coverage for the case genes from clinvar also after this
  ctrl_cov = attributes(matchcontrols)$coverage
  case_cov = ctrl_cov[symbol%in%genes]
  case_cov[,c("over_10", "group") := .(1, "cln")]
  coverage = unique(rbind(coverage, case_cov))


  if(!inframes) max_inframe_size=0

  cases = cases[cases[[filtertype]] < maxmaf & n_res <= max_inframe_size]

  # need to get the pretend ORs for the case genes
  # do this by calculating the freq for controls, then choose ss to get OR ~5:10 for cases
  ctrl_freqs = matchcontrols[symbol %in% genes, .(count=sum(ac)), by=symbol]
  ctrl_freqs = merge(ctrl_freqs, ctrl_cov[,.(group=unique(group), over_10=mean(over_10)), by=symbol], by="symbol")
  ctrl_freqs[,freq:=count/(over_10 * ifelse(group=="e", 125748, 15708))]

  case_counts = cases[symbol %in% genes,.(count=sum(ac)), by=symbol]
  case_counts = merge(case_counts, ctrl_freqs[,.(symbol, freq)], by="symbol")
  case_counts[,ss:=round(count/(freq*runif(.N, 1, 2)))]

  cln_ss = unique(case_counts[,.(symbol, ss)])

  ss = coverage[!symbol%in%genes,.(ss=ifelse(unique(group)=="e", 125748, 15708)), by=symbol]
  ss = rbind(ss, cln_ss)

  cases = cases[,.(aff=1, symbol, protein_position, ac, group)]

  setattr(cases, "coverage", coverage)
  setattr(cases, "ss", ss)
  setattr(cases, "provenance", "auto")

  return(cases)

}

