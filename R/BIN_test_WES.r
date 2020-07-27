#' Calculate a p-values from the BIN-test across the whole-exome
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a data.table: symbol, p-value, covstats (optional)
#'
#' @param cases case data in format: data.table(aff, symbol, protein_position, ac)
#' @param controls control data in format: data.table(aff, symbol, protein_position, ac)
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param covstats include additional coverage information as a column in results?
#' @param messages print messages to the terminal (e.g. how many variants removed by coverage)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details
#' This method splits each gene into k = n ^ (2/5) bins where n is the total number of observed variants in the gene. To control for uneven coverage cell counts are adjusted by the reciprocal of the mean 10X coverage across the bin.
#' @export

# adjust allele count need to document
# added continuity correction - need to document this and also add to bintest single function
BIN_test_WES = function(cases, controls, case_coverage=NULL, control_coverage=NULL, cov_threshold=0.5, covstats=F, messages=T, ac_downweight=0.8){


  case_prov = attributes(cases)$provenance
  control_prov = attributes(controls)$provenance

  case_auto = !is.null(case_prov) && case_prov == "auto"
  control_auto = !is.null(control_prov) && control_prov == "auto"

  if(case_auto) case_coverage = attributes(cases)$coverage
  if(control_auto) control_coverage = attributes(controls)$coverage


  null1 = !is.null(case_coverage)
  null2 = !is.null(control_coverage)

  if(covstats & !(null1 | null2)){
    if(messages) message("Cannot provide covstats when no coverage files supplied.")
    covstats = F
  }


  dataset = rbind(cases, controls, fill=T)

  if(nrow(dataset) == 0){
    stop("dataset has zero rows.")
  }

  gene_counts = dataset[,.(cases=sum(aff==1), controls=sum(aff==0)), by=.(symbol)]

  zero_cases = gene_counts[cases==0, symbol]
  if(length(zero_cases) > 0){
    if(messages) message("Removed ", length(zero_cases), " genes with zero case variants")
    dataset = dataset[!symbol %in% zero_cases]
  }
  zero_controls = gene_counts[controls==0, symbol]
  if(length(zero_controls) > 0){
    if(messages) message("Removed ", length(zero_controls), " genes with zero controls variants")
    dataset = dataset[!symbol %in% zero_controls]
  }


  newdata = remove_uncovered_residues(dataset, case_coverage, control_coverage, cov_threshold, messages)
  dataset = newdata[[1]]
  case_coverage = newdata[[2]]
  control_coverage = newdata[[3]]
  blacklist = newdata[[4]]


  # check input and restructure
  if(nrow(dataset) == 0){
    stop(paste0("Zero variants left after excluded for low coverage regions (threshold = ", cov_threshold, ")"))
  }
  dataset = dataset[rep(1:.N, ceiling(ac^(ac_downweight)))]

  # params
  make_breaks = function(minpos, maxpos, n){
    bins = round(n ^ (2/5) + 1)
    if(bins<4) bins = 4
    breaks = seq(minpos-1, maxpos+1, length.out=bins)
    breaks2 = rep(breaks, each=2)
    breaks3 = data.table(matrix(breaks2[-c(1, length(breaks2))], ncol=2, byrow=T))
    breaks3[,.(start=ceiling(V1), end=floor(V2))]
  }
  break_points = dataset[,make_breaks(min(protein_position), max(protein_position), .N), by=symbol]
  setkey(break_points, symbol, start, end)

  dataset[,c("start", "end"):=.(protein_position, protein_position)]
  setkey(dataset, symbol, start, end)


  break_points$case_counts = foverlaps(break_points, dataset[aff==1], which=T, nomatch=NA)[, ifelse(all(is.na(yid)), 0L, .N), by=xid]$V1
  break_points$control_counts = foverlaps(break_points, dataset[aff==0], which=T, nomatch=NA)[, ifelse(all(is.na(yid)), 0L, .N), by=xid]$V1

  break_points[,case_counts := ifelse(case_counts==0, 0.5, case_counts)]
  break_points[,control_counts := ifelse(control_counts==0, 0.5, control_counts)]

  if(null1){

    case_coverage[,c("start", "end") := .(protein_position, protein_position)]
    setkey(case_coverage, symbol, start, end)

    break_points$case_weight = foverlaps(break_points, case_coverage)[,mean(over_10, na.rm=T), by=.(symbol, i.start)]$V1
    break_points[,case_counts := round(case_counts/case_weight)]

  } else break_points[,case_weight:=1]


  if(null2){

    control_coverage[,c("start", "end") := .(protein_position, protein_position)]
    setkey(control_coverage, symbol, start, end)

    break_points$control_weight = foverlaps(break_points, control_coverage)[,mean(over_10, na.rm=T), by=.(symbol, i.start)]$V1
    break_points[,control_counts := round(control_counts/control_weight)]

  } else break_points[,control_weight:=1]



  calc_p = function(case_counts, control_counts){

    contig = rbind(case_counts, control_counts)
    index = which(apply(contig, 2, function(x) sum(x, na.rm=T) == 0 | any(is.na(x))))
    if(length(index) > 0) contig = contig[, -index, drop=F]
    if(ncol(contig) < 3) NA_real_ else suppressWarnings(chisq.test(contig)$p.value)

  }

  break_points[,BIN.test_pvalue:=calc_p(case_counts, control_counts), by=symbol]


  if(null1 | null2){
    break_points = break_points[,.(BIN.test_pvalue=unique(BIN.test_pvalue), case_mean_cov = mean(case_weight, na.rm=T), control_mean_cov = mean(control_weight, na.rm=T), nbins=.N,
                               NAbins=sum(!complete.cases(.SD)), lowcovbins_cases=sum(case_weight < 0.8, na.rm=T)/.N, lowcovbins_controls=sum(control_weight < 0.8, na.rm=T)/.N), by=symbol]

    break_points[,cov_flag1:=ifelse(case_mean_cov < 0.6 | control_mean_cov < 0.6, "vhigh",
                                ifelse(case_mean_cov < 0.7 | control_mean_cov < 0.7, "high",
                                       ifelse(case_mean_cov < 0.8 | control_mean_cov < 0.8, "moderate", "pass")))]
    break_points[,cov_flag2:=ifelse(lowcovbins_cases > 0.5 | lowcovbins_controls > 0.5, "vhigh",
                                ifelse(lowcovbins_cases > 0.3 | lowcovbins_controls > 0.3, "high",
                                       ifelse(lowcovbins_cases > 0.1 | lowcovbins_controls > 0.1, "moderate", "pass")))]
  } else break_points = break_points[,.(BIN.test_pvalue=unique(BIN.test_pvalue)), by=symbol]


  if(covstats){

    if(!is.null(blacklist)) break_points = merge(break_points, blacklist, all.x=T)

    break_points[is.na(ranges), ranges:=""]
    break_points[,c("case_mean_cov", "control_mean_cov"):=.(signif(case_mean_cov, 2), signif(control_mean_cov, 2))]

    break_points[,covstats:=paste0(",xcov=", case_mean_cov, ",ycov=", control_mean_cov, ",nbins=", nbins-NAbins, "/", nbins,
                                   ",%x_lc_bins=", signif(lowcovbins_cases, 2), ",%y_lc_bins=", signif(lowcovbins_controls, 2), ",rem=(", ranges, ")")]

  }

  data("nresidues")
  break_points = merge(break_points, nresidues, by="symbol", all.x=T)
  break_points[,pl_flag:=ifelse(nresidues > 3000, "vhigh", ifelse(nresidues > 2000, "high", ifelse(nresidues > 1000, "moderate", "pass")))]


  return_cols = c("symbol", "BIN.test_pvalue", "pl_flag")
  if(null1 | null2) return_cols = c(return_cols, "cov_flag1", "cov_flag2")
  if(covstats) return_cols = c(return_cols, "covstats")

  break_points[,return_cols, with=F]

}
