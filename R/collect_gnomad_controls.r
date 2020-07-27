#' Automatically get controls for BIN-test and ClusterBurden
#' GnomAD exomes or genomes v2 (missense only)
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return a data.table of GnomAD controls
#'
#' @param genenames genenames to collect controls for (NULL = all genes).
#' @param dataset GnomAD exomes or genomes? ("exome", "genome")
#' @param switch_dataset_threshold coverage threshold below which to switch control group for that gene (see Details)
#' @param inframes include inframes as missense variants? (default=TRUE)
#' @param max_inframe_size maximum inframe size to include to a maximum of 3 (default=3)
#' @param filtertype type of allele frequency to filter MAF; "global", "popmax" or "strict" where global is the total allele frequency in controls, popmax is the maximum allele frequency in any ancestry group in GnomAD excluding ASJ, FIN, OTH and strict is the maximum allele frequency in any ancestry group in gnomAD or globally across GnomAD genomes
#' @param maxmaf maf threshold using for frequency filtering default is 0.0001 (i.e. 0.1%)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Retrieve missense variants for GnomAD exomes or genomes v2 to use in an association analysis.
#' Columns in resulting dataset include genename (symbol), protein position and allele count.
#'
#' For switch_dataset_threshold the coverage for each gene in gnomAD exomes and genomes are calculated as the mean 10X coverage across all
#' bases in the exonic regions for that gene. If this input is not 0 then for all genes with coverage below the inputs value will switch to
#' the other control group if it has better coverage. For example, if 0.9 (90%) is selected and the dataset chosen is exomes, for all genes with
#' less than 90% coverage in exomes (n=4743), GnomAD genomes will be used instead if it has better coverage (n=4550).
#'
#' The sample sizes for each gene are calculated as the mean number of samples with at least 10X coverage across each base in the exonic regions
#' for that gene multiplied by the complete cohort size of gnomAD exomes (125748) or GnomAD genomes (15708). These sample sizes are attached as an
#' attribute to the control group which can be accessed by: attributes(controls)$ss.
#'
#' @export
#' @examples
#' # Filtering for cases must match filtering for controls
#' # for the simplest scenario using all filtering defaults
#' # Then for the genes of interest e.g. MYH7 and TNNI3
#'
#' controls = collect_gnomad_controls(c("MYH7", "TNNI3"))

collect_gnomad_controls = function(genenames=NULL, dataset="exome", switch_dataset_threshold=0, inframes=T, max_inframe_size=3, filtertype="pci95_strict", maxmaf=0.0001, messages=T){

  data("x10_gnomad")

  if(!is.null(genenames)) x10_gnomad = x10_gnomad[symbol %in% genenames]

  # find the chosen control group based on coverage
  if(dataset=="exome"){

    if(switch_dataset_threshold > 0){

      x10_gnomad[,chosen:=ifelse(is.na(x10_e2) | (!is.na(x10_g2) & x10_e2 < switch_dataset_threshold & x10_e2 < x10_g2), "g", "e")]

      switched = sum(x10_gnomad$chosen=="g")
      if(switched > 0 & messages) message(paste0(switched, " genes switched to GnomAD genomes for better coverage."))

    } else x10_gnomad[,chosen:="e"]

  } else if(dataset=="genome"){


    if(switch_dataset_threshold > 0){

      x10_gnomad[,chosen:=ifelse(is.na(x10_g2) | (!is.na(x10_e2) & x10_g2 < switch_dataset_threshold & x10_g2 < x10_e2), "e", "g")]

      switched = x10_gnomad[,sum(chosen=="e")]
      if(switched > 0 & messages) message(paste0(switched, " genes switched to GnomAD exomes for better coverage."))

    } else x10_gnomad[,chosen:="g"]

  }

  # get the variant data and coverage file for the correct control group for each gene
  controls = NULL
  coverage = NULL

  if(any(x10_gnomad$chosen == "e")){

    data("econtrols")
    econtrols = econtrols[af_e2 <= maxmaf]
    controls = rbind(controls, cbind(econtrols[symbol%in%x10_gnomad[chosen=="e", symbol]], group="e"))

    data("exome_cov")
    coverage = rbind(coverage, exome_cov[symbol%in%x10_gnomad[chosen=="e", symbol], .(symbol, protein_position, over_10, group="e")])

  }

  if(any(x10_gnomad$chosen == "g")){

    data("gcontrols")
    gcontrols = gcontrols[af_g2 <= maxmaf]
    controls = rbind(controls, cbind(gcontrols[symbol%in%x10_gnomad[chosen=="g", symbol]], group="g"), fill=T)

    data("genome_cov")
    coverage = rbind(coverage, genome_cov[symbol%in%x10_gnomad[chosen=="g", symbol], .(symbol, protein_position, over_10, group="g")])
  }

  ss = coverage[,.(ss=ifelse(unique(group)=="e", 125748, 15708)), by=symbol]

  # filter variant data
  if(!inframes) max_inframe_size=0

  controls = controls[controls[[filtertype]] < maxmaf & n_res <= max_inframe_size]

  if(messages & nrow(controls) == 0){
    message("No variants left after filtering")
  }

  controls = controls[,.(aff=0, symbol, protein_position, ac, group)]

  setattr(controls, "coverage", coverage)
  setattr(controls, "ss", ss)
  setattr(controls, "provenance", "auto")

  return(controls)

}
